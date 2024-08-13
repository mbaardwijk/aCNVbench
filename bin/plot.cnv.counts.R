#!/usr/bin/env Rscript
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(ggpubr))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(purrr))
suppressPackageStartupMessages(library(grid))
#suppressPackageStartupMessages(library(reshape2))
suppressPackageStartupMessages(library(stringr))
suppressPackageStartupMessages(library(optparse))
suppressPackageStartupMessages(library(RColorBrewer))
suppressPackageStartupMessages(library(bedr))

VALID.SV.TYPES <- c('BND', 'CNV', 'DEL', 'DUP', 'INS', 'INV')
CALLER.COLUMNS <- c('Sample_ID', 'Chr', 'Start', 'End', 'Type', 'Length')
CALLER.GR.COLUMNS <- c('Sample_ID', 'seqnames', 'start', 'end', 'Type', 'Length')
POSITION.COLUMNS.GOLDSTANDARD <- c('CHROM', 'POS', 'END')
POSITION.COLUMNS.CALLER <- c('Chr', 'Start', 'End')
AUTOSOMAL.CHRS <- paste0("chr", seq(1:22))
FRACTION.OVERLAP <- 0.5
LENGTH.LARGE.CNVS <- 100000

# define set colors for all callers
light_rhg_cols <- c("PennCNV"="#8c9cd5", "QuantiSNP"="#e68d8d", "iPattern"="#8bcf8c", "EnsembleCNV"="#f5f4a0" ,"R-GADA"="#c7b497", "Gold Standard"="#e9f3fc")
dark_rhg_cols <- c("PennCNV"="#0e31a6", "QuantiSNP"="#cd1b1b", "iPattern"="#089908", "EnsembleCNV"="#eae83b", "R-GADA"="#8b6228", "Gold Standard"="#2fa9e2")

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

getCNVNumbers <- function(cnv.calls){
  small <- subset(cnv.calls, Length < LENGTH.LARGE.CNVS)
  large <- subset(cnv.calls, Length >= LENGTH.LARGE.CNVS)
  dels <- subset(cnv.calls, Type %in% c("Loss", "cxLoss", "DEL"))
  dups <- subset(cnv.calls, Type %in% c("Gain", "cxGain", "DUP"))
  # determine number of cnvs for each category and caller
  callers <- unique(cnv.calls$Caller)
  samples <- unique(cnv.calls$Sample_ID)
  cnv.call.numbers <- data.frame()
  for(sample in samples){
    for(caller in callers){
      n.cnvs <- nrow(subset(cnv.calls, Sample_ID == sample & Caller == caller))
      n.small <- nrow(subset(small, Sample_ID == sample & Caller == caller))
      n.large <- nrow(subset(large, Sample_ID == sample & Caller == caller))
      n.dels <- nrow(subset(dels, Sample_ID == sample & Caller == caller))
      n.dups <- nrow(subset(dups, Sample_ID == sample & Caller == caller))
      cnv.call.numbers[sample, paste0(caller, ".All.cnvs")] <- n.cnvs
      cnv.call.numbers[sample, paste0(caller, ".Small")] <- n.small
      cnv.call.numbers[sample, paste0(caller, ".Large")] <- n.large
      cnv.call.numbers[sample, paste0(caller, ".Dels")] <- n.dels
      cnv.call.numbers[sample, paste0(caller, ".Dups")] <- n.dups
    }
  }
  #cnv.call.numbers$Sample_ID <- rownames(cnv.call.numbers)
  return(cnv.call.numbers)
}

plotCNVCount <- function(cnv.calls, color.palette, filename){
  cnv.call.numbers <- getCNVNumbers(cnv.calls)
  # separate the different categories
  all.cnv.numbers <- select(cnv.call.numbers, contains(".All.cnvs"))
  small.cnv.numbers <- select(cnv.call.numbers, contains(".Small"))
  large.cnv.numbers <- select(cnv.call.numbers, contains(".Large"))
  del.cnv.numbers <- select(cnv.call.numbers, contains(".Dels"))
  dup.cnv.numbers <- select(cnv.call.numbers, contains(".Dups"))
  colnames(all.cnv.numbers) <- gsub(".All.cnvs", "", colnames(all.cnv.numbers))
  colnames(small.cnv.numbers) <- gsub(".Small", "", colnames(small.cnv.numbers))
  colnames(large.cnv.numbers) <- gsub(".Large", "", colnames(large.cnv.numbers))
  colnames(del.cnv.numbers) <- gsub(".Dels", "", colnames(del.cnv.numbers))
  colnames(dup.cnv.numbers) <- gsub(".Dups", "", colnames(dup.cnv.numbers))
  all.cnv.numbers <- melt(all.cnv.numbers, value.name="Count", variable.name="Caller")
  small.cnv.numbers <- melt(small.cnv.numbers, value.name="Count", variable.name="Caller")
  large.cnv.numbers <- melt(large.cnv.numbers, value.name="Count", variable.name="Caller")
  del.cnv.numbers <- melt(del.cnv.numbers, value.name="Count", variable.name="Caller")
  dup.cnv.numbers <- melt(dup.cnv.numbers, value.name="Count", variable.name="Caller")
  # prepare for plotting
  all.cnv.numbers$Caller <- factor(all.cnv.numbers$Caller,
                                   levels = c("PennCNV", "QuantiSNP", "iPattern", "EnsembleCNV", "R-GADA", "Gold Standard"),ordered = TRUE)
  small.cnv.numbers$Caller <- factor(small.cnv.numbers$Caller,
                                     levels = c("PennCNV", "QuantiSNP", "iPattern", "EnsembleCNV", "R-GADA", "Gold Standard"),ordered = TRUE)
  large.cnv.numbers$Caller <- factor(large.cnv.numbers$Caller,
                                     levels = c("PennCNV", "QuantiSNP", "iPattern", "EnsembleCNV", "R-GADA", "Gold Standard"),ordered = TRUE)
  dup.cnv.numbers$Caller <- factor(dup.cnv.numbers$Caller,
                                   levels = c("PennCNV", "QuantiSNP", "iPattern", "EnsembleCNV", "R-GADA", "Gold Standard"),ordered = TRUE)
  del.cnv.numbers$Caller <- factor(del.cnv.numbers$Caller,
                                   levels = c("PennCNV", "QuantiSNP", "iPattern", "EnsembleCNV", "R-GADA", "Gold Standard"),ordered = TRUE)
  write.table(all.cnv.numbers, file="all.cnv.numbers.txt", sep="\t", quote=F, row.names=F)
  write.table(del.cnv.numbers, file="del.cnv.numbers.txt", sep="\t", quote=F, row.names=F)
  write.table(dup.cnv.numbers, file="dup.cnv.numbers.txt", sep="\t", quote=F, row.names=F)
  # generate plots
  p1 <- ggplot(all.cnv.numbers, aes(x=Caller, y=Count, fill=Caller)) + 
    geom_boxplot() +
    stat_summary(fun=mean, geom="point", shape=20, size=3, color="red", fill="red") +
    labs(x="", y="", title = "All CNVs") +
    theme_bw() +
    theme(axis.text.x = element_text(angle=90, vjust=0.5, hjust=1), legend.position = "none", plot.title = element_textbox(hjust = 0.5, margin = margin(t = 5, b = 5))) +
    scale_fill_manual(values = color.palette) +
    scale_y_continuous(limits = c(0, 6500))

  p2 <- ggplot(del.cnv.numbers, aes(x=Caller, y=Count, fill=Caller)) + 
    geom_boxplot() +
    stat_summary(fun=mean, geom="point", shape=20, size=3, color="red", fill="red")  +
    labs(x="", y="", title = "Deletions") +
    theme_bw() +
    theme(axis.text.x = element_text(angle=90, vjust=0.5, hjust=1), legend.position = "none", plot.title = element_textbox(hjust = 0.5, margin = margin(t = 5, b = 5))) +
    scale_fill_manual(values = color.palette)+
    scale_y_continuous(limits = c(0, 3500))
  
  p3 <- ggplot(dup.cnv.numbers, aes(x=Caller, y=Count, fill=Caller)) + 
    geom_boxplot() +
    stat_summary(fun=mean, geom="point", shape=20, size=3, color="red", fill="red")  +
    labs(x="", y="", title = "Duplications") +
    theme_bw() +
    theme(axis.text.x = element_text(angle=90, vjust=0.5, hjust=1), legend.position = "none", plot.title = element_textbox(hjust = 0.5, margin = margin(t = 5, b = 5))) +
    scale_fill_manual(values = color.palette)+
    scale_y_continuous(limits = c(0, 3500))
  
  p4 <- ggplot(small.cnv.numbers, aes(x=Caller, y=Count, fill=Caller)) + 
    geom_boxplot() +
    stat_summary(fun=mean, geom="point", shape=20, size=3, color="red", fill="red")  +
    labs(x="", y="", title = "Small CNVs") +
    theme_bw() +
    theme(axis.text.x = element_text(angle=90, vjust=0.5, hjust=1), legend.position = "none", plot.title = element_textbox(hjust = 0.5, margin = margin(t = 5, b = 5))) +
    scale_fill_manual(values = color.palette)+
    scale_y_continuous(limits = c(0, 3500))
  
  p5 <- ggplot(large.cnv.numbers, aes(x=Caller, y=Count, fill=Caller)) + 
    geom_boxplot() +
    stat_summary(fun=mean, geom="point", shape=20, size=3, color="red", fill="red")  +
    labs(x="", y="", title = "Large CNVs") +
    theme_bw() +
    theme(axis.text.x = element_text(angle=90, vjust=0.5, hjust=1), legend.position = "none", plot.title = element_textbox(hjust = 0.5, margin = margin(t = 5, b = 5))) +
    scale_fill_manual(values = color.palette)+
    scale_y_continuous(limits = c(0, 3500))
   
  figure <- ggarrange(p1, p2, p3, p4, p5, ncol=3, nrow=2, common.legend = TRUE, legend="top")
  annotated_figure <- annotate_figure(figure, left = textGrob("CNV Count", rot = 90, vjust = 1, gp = gpar(cex = 1.3)))
  tiff(filename, units="cm", width=24, height=24, res=600)
  grid.draw(annotated_figure)
  dev.off()
}

# determine individual CNVs based on ALT and GT fields from VCF
parseCNVGenotype <- function(vcf.cnv.data, gt.column.name){
  vcf.cnv.data <- separate(data = vcf.cnv.data, col = gt.column.name, into = c("GT.1", "GT.2"))
  vcf.cnv.data$Type <- "REF"
  vcf.cnv.data$Type[vcf.cnv.data$GT.1 == 1 | vcf.cnv.data$GT.2 == 1] <- vcf.cnv.data$ALT[vcf.cnv.data$GT.1 == 1 | vcf.cnv.data$GT.2 == 1]
  vcf.cnv.data$Type <- gsub("<", "", vcf.cnv.data$Type)
  vcf.cnv.data$Type <- gsub(">", "", vcf.cnv.data$Type)
  return(vcf.cnv.data)
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
  
# generate plots
plotCNVCount(cnv.calls, light_rhg_cols, filename=paste0(opt$output,"_CNV_counts.tiff"))

# Execute analysis per super population
if(!is.na(opt$population)){
  population.info <- read.table(opt$population, sep="\t", header=T)
  for(population in unique(population.info$super_pop)){
    population.samples <- population.info[population.info$super_pop == population, 'sample']
    print(head(cnv.calls))
    
    population.cnv.calls <- cnv.calls[cnv.calls$Sample_ID %in% population.samples, ]
    #plotCNVSize(cnv.calls, dark_rhg_cols, filename=paste0(opt$output,"_CNV_sizes_", population, ".tiff"))
    plotCNVCount(population.cnv.calls, light_rhg_cols, filename=paste0(opt$output, "_CNV_counts_", population, ".tiff"))
    
  }
}
