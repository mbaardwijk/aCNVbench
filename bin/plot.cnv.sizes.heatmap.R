suppressPackageStartupMessages(library(optparse))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(viridis))
suppressPackageStartupMessages(library(ggpubr))
suppressPackageStartupMessages(library(grid))
suppressPackageStartupMessages(library(GenomicRanges))
suppressPackageStartupMessages(library(bedr))
suppressPackageStartupMessages(library(stringr))
suppressPackageStartupMessages(library(reshape2))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(scales))

COLUMN.NAMESS <- c('Sample_ID', 'Chr', 'Start', 'End', 'Type')

# define set colors for all callers
light_rhg_cols <- c("PennCNV"="#8c9cd5", "QuantiSNP"="#e68d8d", "iPattern"="#8bcf8c", "EnsembleCNV"="#f5f4a0" ,"R-GADA"="#c7b497", "Gold Standard"="#e9f3fc")

option_list = list(
  make_option(c("-c", "--callfiles"), action = "store", type = "character", default = NA,
              help = "List of callfiles with CNV data with columns Sample_ID, Chr, Start, End and Type"),
  make_option(c("-g", "--goldstandard"), action = "store", type = "character", default = NA,
              help = "csv file with gold standard cnvs"),
  make_option(c("-o", "--output"), action = "store", type = "character", default = "results",
              help = "Prefix for outputfiles generated during the analysis")
)

# Get all command line options
opt = parse_args(OptionParser(option_list = option_list))

plotCNVSize <- function(unique.cnv.calls, color.palette, filename){
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
}

getAllCalls <- function(file.list, goldstandard.data){
  # read cnv call files and add to single dataframe
  cnv.calls <- data.frame(matrix(ncol = 6, nrow=0, dimnames=list(NULL, c("Sample_ID","Chr", "Start", "End", "Type", "Caller"))))
  for(file in file.list){
    cnvs <- read.table(file, header=T, sep="\t")
    cnvs$Caller <- gsub(".mapped.id.txt", "", gsub(".results.txt", "", basename(file)))
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
  goldstandard.calls <- goldstandard.data[, COLUMN.NAMESS]
  colnames(goldstandard.calls) <- COLUMN.NAMESS
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
goldstandard.cnvs <- read.csv(goldstandard.filename, header=T)
cnv.calls <- getAllCalls(caller.files, goldstandard.cnvs)

# Execute analysis for all samples
unique.cnv.calls <- unique(cnv.calls[c("Chr", "Start", "End", "Type", "Caller")])
unique.cnv.calls$Length <- unique.cnv.calls$End - unique.cnv.calls$Start

# generate plots
plotCNVSize(unique.cnv.calls, light_rhg_cols, filename=paste0(opt$output,"_CNV_counts.tiff"))
