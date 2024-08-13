library(optparse)
library(dplyr)
library(grid)
library(gridExtra)
library(ggplot2)
library(cowplot)
library(ggpubr)
library(reshape2)
library(data.table)
library(bedtoolsr)
library(ComplexUpset)
library(viridis)

FRACTION.OVERLAP <- 0.5

option_list = list(
  make_option(c("-c", "--callfiles"), action = "store", type = "character", default = NA,
              help = "List of callfiles with CNV data with columns Sample_ID, Chr, Start, End and Type"),
  make_option(c("-o", "--output"), action = "store", type = "character", default = "results",
              help = "Prefix for outputfiles generated during the analysis")
)

# Get all command line options
opt = parse_args(OptionParser(option_list = option_list))

# create textbox backgrounds for title
element_textbox <- function(...) {
  el <- element_text(...)
  class(el) <- c("element_textbox", class(el))
  el
}

element_grob.element_textbox <- function(element, ...) {
  text_grob <- NextMethod()
  rect_grob <- element_grob(calc_element("strip.background", theme_bw()))
  
  ggplot2:::absoluteGrob(
    grid::gList(
      element_grob(calc_element("strip.background", theme_bw())),
      text_grob
    ),
    height = grid::grobHeight(text_grob), 
    width = grid::unit(1, "npc")
  )
}

# read cnv call files and add to single dataframe
callfiles <- opt$callfiles
caller.files <- unlist(strsplit(callfiles, "\\s+"))
caller.files <- unlist(strsplit(opt$callfiles, "\\s+"))
cnv.calls <- data.frame(matrix(ncol = 6, nrow=0, dimnames=list(NULL, c("Sample_ID","Chr", "Start", "End", "Type", "Caller"))))
for(file in caller.files){
  cnvs <- read.table(file, header=T, sep="\t")
  cnvs$Caller <- gsub(".mapped.id.txt", "", gsub("results.", "", basename(file)))
  cnvs$Caller <- gsub("rgada", "R-GADA", cnvs$Caller)
  cnvs$Caller <- gsub("penncnv", "PennCNV", cnvs$Caller)
  cnvs$Caller <- gsub("quantisnp", "QuantiSNP", cnvs$Caller)
  cnvs$Caller <- gsub("ipattern", "iPattern", cnvs$Caller)
  cnvs$Caller <- gsub("ensemblecnv", "EnsembleCNV", cnvs$Caller)
  cnvs$Chr <- gsub("chr", "", cnvs$Chr)
  cnvs <- cnvs[,c("Sample_ID","Chr", "Start", "End", "Type", "Caller")]
  cnv.calls <- rbind(cnv.calls, cnvs)
}

# remove all cnvs from non-autosomal chromosomes
cnvs <- subset(cnv.calls, Chr %in% 1:22)
cnvs <- cnvs[complete.cases(cnvs),]
unique.cnvs <- unique(cnvs[,c("Chr", "Start", "End", "Type", "Caller")])

# add platform ids to datasets
penncnv <- subset(unique.cnvs, Caller == "PennCNV")
quantisnp <- subset(unique.cnvs, Caller == "QuantiSNP")
ipattern <- subset(unique.cnvs, Caller ==  "iPattern")
rgada <- subset(unique.cnvs, Caller == "R-GADA")
ensemblecnv <- subset(unique.cnvs, Caller == "EnsembleCNV")

# segment cnvs
chrs = (1:22)
cnv.segments = data.frame(Chr=character(0), Start=numeric(0), End=numeric(0), Type=character(0))
for(cur.chr in chrs){
  print(cur.chr)
  segments = data.frame(Chr=character(0), Start=numeric(0), End=numeric(0))
  cnvs.chr = subset(unique.cnvs, Chr == cur.chr)
  breakpoints = unique(c(cnvs.chr$Start, cnvs.chr$End))
  # sort breakpoints from low to high
  breakpoints = sort(breakpoints)
  for(i in 1:(length(breakpoints)-1)){
    segment = list(Chr=cur.chr, Start=breakpoints[i], End=breakpoints[i+1])
    segments = rbind(segments, segment)
  }
  # duplicate segments for deletions and duplications
  n.segments = nrow(segments)
  all.segments = rbind(segments, segments, segments)
  all.segments$Type[1:n.segments] <- "DEL"
  all.segments$Type[(n.segments+1):(2*n.segments)] <- "DUP"
  # add to all segments
  cnv.segments <- rbind(cnv.segments, all.segments)
}

cnv.segments$PennCNV <- NA
cnv.segments$QuantiSNP <- NA
cnv.segments$iPattern <- NA
cnv.segments$`R-GADA` <- NA
cnv.segments$EnsembleCNV <- NA
cnv.segments$CNV_ID <- paste0(cnv.segments$Chr, ":", cnv.segments$Start, "-", cnv.segments$End, "_", cnv.segments$Type)
print(unique(cnv.segments$Chr))

# separate deletions and duplications
cnvrs.dels <- subset(cnv.segments, Type == "DEL")
cnvrs.dups <- subset(cnv.segments, Type == "DUP")
cnvrs.subtypes <- list(DEL = cnvrs.dels, DUP = cnvrs.dups)

# for each caller and each cnvr subtype, iterate over cnvs and determine overlap
for (caller in unique(unique.cnvs$Caller)){
  print(caller)
  caller.dels <- subset(unique.cnvs, Caller == caller & Type == "DEL")
  caller.dups <- subset(unique.cnvs, Caller == caller & Type == "DUP")
  caller.subtypes <- list(DEL = caller.dels, DUP = caller.dups)
  for (subtype in c("DEL", "DUP")){
    print(subtype)
    cnvrs.subtype <- cnvrs.subtypes[[subtype]]
    caller.cnvs <- caller.subtypes[[subtype]]
    for (cur.chr in 22:1){
      print(cur.chr)
      cnvrs.subtype.chr <- subset(cnvrs.subtype, Chr == cur.chr)
      caller.cnvs.chr <- subset(caller.cnvs, Chr == cur.chr)
      overlapping.pairs <- bt.intersect(a = cnvrs.subtype.chr[,c("Chr", "Start", "End")],
                                        b = caller.cnvs.chr[,c("Chr", "Start", "End")],
                                        wa = T,
                                        f = FRACTION.OVERLAP
      )
      colnames(overlapping.pairs) <- c("Chr", "Start", "End")
      overlapping.pairs <- overlapping.pairs[!duplicated(overlapping.pairs),]
      overlapping.pairs$Type <- subtype
      overlapping.pairs[[caller]] <- T
      overlapping.pairs$CNV_ID <- paste0(overlapping.pairs$Chr, ":", overlapping.pairs$Start, "-", overlapping.pairs$End, "_", overlapping.pairs$Type)
      cnv.segments[[caller]][cnv.segments$CNV_ID %in% overlapping.pairs$CNV_ID] <- T
    }
  }
}

cnv.segments$Length <- cnv.segments$End - cnv.segments$Start

print(unique(cnv.segments$Chr))

# set bin size
breaks <- c(0, 500, 1000, 10000, 100000, 1000000)
tags <- c("100-1,000", "10-100", "1-10", "0.5-1", "0-0.5")

cnv.segments$bin <- cut(cnv.segments$Length,
                            breaks=breaks)

cnv.segments$bin <- factor(cnv.segments$bin,
                levels = c("(1e+05,1e+06]", "(1e+04,1e+05]", "(1e+03,1e+04]", "(500,1e+03]", "(0,500]"))
print(colnames(cnv.segments))
penncnv.cnvrs <- pull(subset(cnv.segments, PennCNV == T), CNV_ID)
quantisnp.cnvrs <- pull(subset(cnv.segments, QuantiSNP == T), CNV_ID)
ipattern.cnvrs <- pull(subset(cnv.segments, iPattern == T), CNV_ID)
rgada.cnvrs <- pull(subset(cnv.segments, `R-GADA` == T), CNV_ID)
ensemblecnv.cnvrs <- pull(subset(cnv.segments, EnsembleCNV == T), CNV_ID)

penncnv.dels <- pull(subset(cnv.segments, PennCNV == T & Type == "DEL"), CNV_ID)
quantisnp.dels <- pull(subset(cnv.segments, QuantiSNP == T & Type == "DEL"), CNV_ID)
ipattern.dels <- pull(subset(cnv.segments, iPattern == T & Type == "DEL"), CNV_ID)
rgada.dels <- pull(subset(cnv.segments, `R-GADA` == T & Type == "DEL"), CNV_ID)
ensemblecnv.dels <- pull(subset(cnv.segments, EnsembleCNV == T & Type == "DEL"), CNV_ID)

penncnv.dups <- pull(subset(cnv.segments, PennCNV == T & Type == "DUP"), CNV_ID)
quantisnp.dups <- pull(subset(cnv.segments, QuantiSNP == T & Type == "DUP"), CNV_ID)
ipattern.dups <- pull(subset(cnv.segments, iPattern == T & Type == "DUP"), CNV_ID)
rgada.dups <- pull(subset(cnv.segments, `R-GADA` == T & Type == "DUP"), CNV_ID)
ensemblecnv.dups <- pull(subset(cnv.segments, EnsembleCNV == T & Type == "DUP"), CNV_ID)

set.cnvrs <- list(PennCNV = penncnv.cnvrs, QuantiSNP = quantisnp.cnvrs, iPattern = ipattern.cnvrs, "R-GADA" = rgada.cnvrs, EnsembleCNV = ensemblecnv.cnvrs)
set.dels <- list(PennCNV = penncnv.dels, QuantiSNP = quantisnp.dels, iPattern = ipattern.dels, "R-GADA" = rgada.dels, EnsembleCNV = ensemblecnv.dels)
set.dups <- list(PennCNV = penncnv.dups, QuantiSNP = quantisnp.dups, iPattern = ipattern.dups, "R-GADA" = rgada.dups, EnsembleCNV = ensemblecnv.dups)

# remove rows with NA for all callers
is.na(cnv.segments) <- F
cnv.segments <- cnv.segments[rowSums(is.na(cnv.segments[,c("PennCNV", "QuantiSNP", "iPattern", "EnsembleCNV", "R-GADA")])) != ncol(cnv.segments[c("PennCNV", "QuantiSNP", "iPattern", "EnsembleCNV", "R-GADA")]), ]
colnames(cnv.segments)[8] <- "R-GADA"

print(unique(cnv.segments$Chr))

show_hide_scale = scale_color_manual(values=c('show'='white', 'hide'='transparent'), guide='none')
viridis_colors = viridis(n=length(tags))
size_scale = scale_fill_manual(values=c(
  `(0,500]`=viridis_colors[1], `(500,1e+03]`=viridis_colors[2],
  `(1e+03,1e+04]`=viridis_colors[3], `(1e+04,1e+05]`=viridis_colors[4],
  `(1e+05,1e+06]`=viridis_colors[5]
),
labels = tags)
tiff("Alt_UpSetPlot_overlap.tiff", units="cm", width=40, height=24, res=600)
upset(cnv.segments, c("PennCNV", "QuantiSNP", "iPattern", "EnsembleCNV", "R-GADA"), name="Caller",
      width_ratio=0.2,
      keep_empty_groups=F,
      queries=list(
        upset_query(set='PennCNV', color='#8c9cd5', fill='#8c9cd5'),
        upset_query(set='QuantiSNP', color='#e68d8d', fill='#e68d8d'),
        upset_query(set='iPattern', color='#8bcf8c', fill='#8bcf8c'),
        upset_query(set='EnsembleCNV', color='#f5f4a0', fill='#f5f4a0'),
        upset_query(set='R-GADA', color='#c7b497', fill='#c7b497')
      ),
      base_annotations=list(
        'Number of CNV segments in set'=(
          intersection_size(
            bar_number_threshold=1,  # show all numbers on top of bars
            width=0.5,   # reduce width of the bars
          )
          # add some space on the top of the bars
          + scale_y_continuous(expand=expansion(mult=c(0, 0.05)))
          + theme(
            # hide grid lines
            panel.grid.major=element_blank(),
            panel.grid.minor=element_blank(),
            # show axis lines
            axis.line=element_line(colour='black'),
            axis.title.y=element_text(size=14)
          )
        )
      ),
      themes=upset_modify_themes(
        list(
          'intersections_matrix'=theme(axis.text.y=element_text(size=14), axis.title.x=element_text(size=14)),
          'overall_sizes'=theme(axis.text.x=element_text(size=14), axis.title.x=element_text(size=14)),
          'Intersection size'=theme(text=element_text(size=14)),
          'default'=theme(text=element_text(size=14))
        )),
      matrix=intersection_matrix(
        geom=geom_point(
          shape='circle filled',
          size=5,
          stroke=0.45
        )
      ),
      annotations =list(
        'Percentage of CNV segments in set'=list(
          aes=aes(x=intersection, fill=bin),
          geom=list(
            geom_bar(stat='count', position='fill', na.rm=TRUE, width=0.5),
            geom_text(
              aes(
                label=!!aes_percentage(relative_to='intersection'),
                color=ifelse(bin %in% c("(0,500]", "(500,1e+03]", "(1e+03,1e+04]"), 'show', 'hide'),
              ),
              size=4,
              stat='count',
              position=position_fill(vjust = .5)
            ),
            labs(fill = "Segment size (kb)"),
            scale_y_continuous(labels=scales::percent_format()),
            show_hide_scale,
            size_scale,
            theme(axis.title = element_text(size = 14),
                  axis.text.y = element_text(size = 14),
                  plot.title = element_text(size=14),
                  # hide grid lines
                  panel.grid.major=element_blank(),
                  panel.grid.minor=element_blank(),
                  # show axis lines
                  axis.line=element_line(colour='black'))
          )
        )
      ),
      set_sizes=(
        upset_set_size()
        + theme(axis.text.x=element_text(angle=90, size=14), axis.ticks.x=element_line())
        + labs(y="Number of CNV segments")
      )
)
dev.off()
