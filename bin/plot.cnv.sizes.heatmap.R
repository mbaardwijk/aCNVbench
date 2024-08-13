library(optparse)
library(ggplot2)
library(dplyr)
library(viridis)
library(ggpubr)
library(grid)
library(GenomicRanges)
library(bedr)
library(stringr)
library(reshape2)
library(data.table)
library(scales)

option_list = list(
  make_option(c("-c", "--callfiles"), action = "store", type = "character", default = NA,
              help = "List of callfiles with CNV data with columns Sample_ID, Chr, Start, End and Type"),
  make_option(c("-g", "--goldstandard"), action = "store", type = "character", default = NA,
              help = "VCF file with gold standard cnvs"),
  make_option(c("--population"), action = "store", type = "character", default = NA,
              help = "Population panel file"),
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

parseGoldStandardVCF <- function(goldstandard.file) {
  invisible(capture.output(goldstandard <- read.vcf(goldstandard.file, split.info = FALSE, split.sample = TRUE)))
  goldstandard.cnvs.end <- str_split(str_split(goldstandard$vcf$INFO, "END=", simplify = TRUE)[,2], ";", simplify = TRUE)[,1]
  goldstandard.cnvs <- cbind(cbind(goldstandard$vcf[, c("CHROM", "POS", "ALT")], goldstandard.cnvs.end), goldstandard$GT)
  colnames(goldstandard.cnvs)[4] <- "END"
  goldstandard.cnvs <- melt(as.data.table(goldstandard.cnvs), id.vars=c("CHROM", "POS", "END", "ALT"), variable.name="SAMPLE", value.name="GT")
  goldstandard.cnvs <- subset(goldstandard.cnvs, GT != "0/0" & GT != "0")
  goldstandard.cnvs <- subset(goldstandard.cnvs, select=-c(GT))
  colnames(goldstandard.cnvs) <- c("CHROM", "POS", "END", "Type", "SAMPLE")
  goldstandard.cnvs$Type <- gsub("<", "", gsub(">", "", goldstandard.cnvs$Type))
  return(goldstandard.cnvs)
}

getAllCalls <- function(file.list, goldstandard.data){
  # read cnv call files and add to single dataframe
  cnv.calls <- data.frame(matrix(ncol = 6, nrow=0, dimnames=list(NULL, c("Sample_ID","Chr", "Start", "End", "Type", "Caller"))))
  for(file in file.list){
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
  # add gold standard cnvs
  samples <- unique(cnv.calls$Sample_ID)
  goldstandard.calls <- goldstandard.data[,c("SAMPLE", "CHROM", "POS", "END", "Type")]
  colnames(goldstandard.calls) <- c("Sample_ID","Chr", "Start", "End", "Type")
  goldstandard.calls <- subset(goldstandard.calls, Sample_ID %in% samples)
  goldstandard.calls$Chr <- gsub("chr", "", goldstandard.calls$Chr)
  goldstandard.calls$Caller <- "Gold Standard"
  goldstandard.calls$End <- sapply(goldstandard.calls$End, as.numeric)
  goldstandard.calls$Start <- sapply(goldstandard.calls$Start, as.numeric)
  cnv.calls <- rbind(cnv.calls, goldstandard.calls)
  # remove all cnvs from non-autosomal chromosomes
  cnv.calls <- subset(cnv.calls, Chr %in% seq(1, 22, 1))
  cnv.calls$Length <- cnv.calls$End - cnv.calls$Start
  cnv.calls <- subset(cnv.calls, Length > 0)
  return(cnv.calls)
}
caller.files <- unlist(strsplit(opt$callfiles, "\\s+"))
goldstandard.filename <- opt$goldstandard

# read the gold standard and caller files
goldstandard.cnvs <- parseGoldStandardVCF(goldstandard.filename)
cnv.calls <- getAllCalls(caller.files, goldstandard.cnvs)

# Execute analysis for all samples
unique.cnv.calls <- unique(cnv.calls[c("Chr", "Start", "End", "Type", "Caller")])
unique.cnv.calls$Length <- unique.cnv.calls$End - unique.cnv.calls$Start

# set bin size
breaks <- c(0, 500, 1000, 10000, 100000, 1000000, 10000000, 100000000, 1000000000, Inf)
tags <- c("0-0.5", "0.5-1", "1-10" , "10-100", "100-1,000", "1,000-10,000", "10,000-100,000", "100,000-1,000,000", ">1,000,000")
unique.cnv.calls$bin <- cut(unique.cnv.calls$Length,
                     breaks=breaks)
unique.cnv.calls$Caller <- factor(unique.cnv.calls$Caller,
                           levels = c("PennCNV", "QuantiSNP", "iPattern", "EnsembleCNV", "R-GADA", "Gold Standard"),ordered = TRUE)
# separate deletions and duplicationss
dels <- subset(unique.cnv.calls, Type == "DEL")
dups <- subset(unique.cnv.calls, Type == "DUP")

print(data.frame(caller=unique.cnv.calls$Caller,bin=factor(unique.cnv.calls$bin)) %>% count(caller,bin,.drop=FALSE)) 

lseq <- function(from=1, to=100000, length.out=6) {
  # logarithmic spaced sequence
  # blatantly stolen from library("emdbook"), because need only this
  exp(seq(log(from), log(to), length.out = length.out))
}

print(lseq(1, 925365, length.out = 5))

p1 <- data.frame(caller=unique.cnv.calls$Caller,bin=factor(unique.cnv.calls$bin)) %>% 
  count(caller,bin,.drop=FALSE) %>% 
  ggplot(aes(x=caller,y=bin,fill=n)) + geom_tile() +
  labs(x="", y="Size range of called CNVs (kb)") +
  theme_classic() +
  theme(axis.text.x = element_text(angle=90, vjust=0.5, hjust=1), legend.position = "none") +
  scale_y_discrete(labels = tags) +
  scale_fill_viridis( trans = 'log', name ="Number of \nunique CNVs", na.value="transparent", breaks=round(lseq(1, 925365, length.out = 5)))

p2 <- data.frame(caller=dels$Caller,bin=factor(dels$bin)) %>% 
  count(caller,bin,.drop=FALSE) %>% 
  ggplot(aes(x=caller,y=bin,fill=n)) + geom_tile() +
  labs(x="", y="Size range of called deletions (kb)") +
  theme_classic() +
  theme(axis.text.x = element_text(angle=90, vjust=0.5, hjust=1), legend.position = "none") +
  scale_y_discrete(labels = tags) +
  scale_fill_viridis( trans = 'log', name ="Number of \nunique CNVs", na.value="transparent", breaks=round(lseq(1, 925365, length.out = 5)))

p3 <- data.frame(caller=dups$Caller,bin=factor(dups$bin)) %>% 
  count(caller,bin,.drop=FALSE) %>% 
  ggplot(aes(x=caller,y=bin,fill=n)) + geom_tile() +
  labs(x="", y="Size range of called duplications (kb)") +
  theme_classic() +
  theme(axis.text.x = element_text(angle=90, vjust=0.5, hjust=1), legend.position = "none") +
  scale_y_discrete(labels = tags) +
  scale_fill_viridis( trans = 'log', name ="Number of \nunique CNVs", na.value="transparent", breaks=round(lseq(1, 925365, length.out = 5)))

figure <- ggarrange(p1, p2, p3, ncol=3, common.legend = T, legend='right', labels=c("A", "B", "C"), label.y=0.05)
tiff(paste0(opt$output, "_CNV_sizes_unique_log.tiff"), units="cm", width=28, height=14, res=600)
grid.draw(figure)
dev.off()

# Execute analysis per super population
if(!is.na(opt$population)){
  for(population in unique(population.info$super_pop)){
    population.samples <- population.info[population.info$super_pop == population, 'sample']
    population.cnv.calls <- cnv.calls[cnv.calls$Sample_ID %in% population.samples, ]
    
    unique.population.cnv.calls <- unique(population.cnv.calls[c("Chr", "Start", "End", "Type", "Caller")])
    unique.population.cnv.calls$Length <- unique.population.cnv.calls$End - unique.population.cnv.calls$Start
    
    # set bin size
    breaks <- c(0, 500, 1000, 10000, 100000, 1000000, 10000000, 100000000, 1000000000, Inf)
    tags <- c("0-0.5", "0.5-1", "1-10" , "10-100", "100-1,000", "1,000-10,000", "10,000-100,000", "100,000-1,000,000", ">1,000,000")
    unique.population.cnv.calls$bin <- cut(unique.population.cnv.calls$Length,
                                breaks=breaks)
    unique.population.cnv.calls$Caller <- factor(unique.population.cnv.calls$Caller,
                                      levels = c("PennCNV", "QuantiSNP", "iPattern", "EnsembleCNV", "R-GADA", "Gold Standard"),ordered = TRUE)
    # separate deletions and duplicationss
    population.dels <- subset(unique.population.cnv.calls, Type == "DEL")
    population.dups <- subset(unique.population.cnv.calls, Type == "DUP")
    
    p1 <- data.frame(caller=unique.population.cnv.calls$Caller,bin=factor(unique.population.cnv.calls$bin)) %>% 
      count(caller,bin,.drop=FALSE) %>% 
      ggplot(aes(x=caller,y=bin,fill=n)) + geom_tile() +
      labs(x="", y="", title="All CNVs", fill="Count") +
      theme_bw() +
      theme(axis.text.x = element_text(angle=90, vjust=0.5, hjust=1), legend.position = "none", plot.title = element_textbox(hjust = 0.5, margin = margin(t = 5, b = 5))) +
      scale_y_discrete(labels = tags) +
      scale_fill_viridis( trans = 'log', breaks = c(0, 20, 403, 8103, 162754), name ="CNV Count\n(log transformed)", limit=range(1, 925365))
    
    p2 <- data.frame(caller=population.dels$Caller,bin=factor(population.dels$bin)) %>% 
      count(caller,bin,.drop=FALSE) %>% 
      ggplot(aes(x=caller,y=bin,fill=n)) + geom_tile() +
      labs(x="", y="", title="Deletions", fill="Count") +
      theme_bw() +
      theme(axis.text.x = element_text(angle=90, vjust=0.5, hjust=1), legend.position = "none", plot.title = element_textbox(hjust = 0.5, margin = margin(t = 5, b = 5))) +
      scale_y_discrete(labels = tags) +
      scale_fill_viridis( trans = 'log', breaks = c(0, 20, 403, 8103, 162754), name ="CNV Count\n(log transformed)", limit=range(1, 925365))
    
    p3 <- data.frame(caller=population.dups$Caller,bin=factor(population.dups$bin)) %>% 
      count(caller,bin,.drop=FALSE) %>% 
      ggplot(aes(x=caller,y=bin,fill=n)) + geom_tile() +
      labs(x="", y="", title="Duplications", fill="Count") +
      theme_bw() +
      theme(axis.text.x = element_text(angle=90, vjust=0.5, hjust=1), legend.position = "none", plot.title = element_textbox(hjust = 0.5, margin = margin(t = 5, b = 5))) +
      scale_y_discrete(labels = tags) +
      scale_fill_viridis( trans = 'log', breaks = c(0, 20, 403, 8103, 162754), name ="CNV Count\n(log transformed)", limit=range(1, 925365))
    
    figure <- ggarrange(p1, p2, p3, ncol=3, common.legend = TRUE, legend="right")
    annotated_figure <- annotate_figure(figure, left = textGrob("CNV size in kb", rot = 90, vjust = 1, gp = gpar(cex = 1.3)))
    tiff(paste0(opt$output, "_CNV_sizes_unique_log_", population, ".tiff"), units="cm", width=28, height=14, res=600)
    grid.draw(annotated_figure)
    dev.off()
  }
}
