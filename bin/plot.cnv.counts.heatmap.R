suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(viridis))
suppressPackageStartupMessages(library(ggpubr))
suppressPackageStartupMessages(library(grid))
suppressPackageStartupMessages(library(optparse))

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

getCNVNumbers <- function(cnv.calls){
  dels <- subset(cnv.calls, Type %in% c("Loss", "cxLoss", "DEL"))
  dups <- subset(cnv.calls, Type %in% c("Gain", "cxGain", "DUP"))
  # determine number of cnvs for each category and caller
  callers <- unique(cnv.calls$Caller)
  samples <- unique(cnv.calls$Sample_ID)
  cnv.call.numbers <- data.frame()
  for(sample in samples){
    for(caller in callers){
      n.cnvs <- nrow(subset(cnv.calls, Sample_ID == sample & Caller == caller))
      n.dels <- nrow(subset(dels, Sample_ID == sample & Caller == caller))
      n.dups <- nrow(subset(dups, Sample_ID == sample & Caller == caller))
      cnv.call.numbers[sample, paste0(caller, ".All.cnvs")] <- n.cnvs
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
  del.numbers <- select(cnv.call.numbers, contains(".Dels"))
  dup.numbers <- select(cnv.call.numbers, contains(".Dups"))
  colnames(all.cnv.numbers) <- gsub(".All.cnvs", "", colnames(all.cnv.numbers))
  colnames(del.numbers) <- gsub(".Dels", "", colnames(del.numbers))
  colnames(dup.numbers) <- gsub(".Dups", "", colnames(dup.numbers))
  all.cnv.numbers <- melt(all.cnv.numbers, value.name="Count", variable.name="Caller")
  del.numbers <- melt(del.numbers, value.name="Count", variable.name="Caller")
  dup.numbers <- melt(dup.numbers, value.name="Count", variable.name="Caller")
  # prepare for plotting
  all.cnv.numbers$Caller <- factor(all.cnv.numbers$Caller,
                                  levels = c("PennCNV", "QuantiSNP", "iPattern", "EnsembleCNV", "R-GADA", "Gold Standard"),ordered = TRUE)
  dup.numbers$Caller <- factor(dup.numbers$Caller,
                                  levels = c("PennCNV", "QuantiSNP", "iPattern", "EnsembleCNV", "R-GADA", "Gold Standard"),ordered = TRUE)
  del.numbers$Caller <- factor(del.numbers$Caller,
                                  levels = c("PennCNV", "QuantiSNP", "iPattern", "EnsembleCNV", "R-GADA", "Gold Standard"),ordered = TRUE)
  write.table(all.cnv.numbers, file="all.cnv.numbers.txt", sep="\t", quote=F, row.names=F)
  write.table(del.numbers, file="del.numbers.txt", sep="\t", quote=F, row.names=F)
  write.table(dup.numbers, file="dup.numbers.txt", sep="\t", quote=F, row.names=F)
  # set bin size
  breaks <- c(1, 10, 50, 100, 1000, Inf)
  tags <- c("0", "1-10", "11-50" , "51-100", "101-1,000", ">1,000")

  all.cnv.numbers$Count <- as.numeric(all.cnv.numbers$Count)
  del.numbers$Count <- as.numeric(del.numbers$Count)
  dup.numbers$Count <- as.numeric(dup.numbers$Count)

  all.cnv.numbers$bin <- cut(all.cnv.numbers$Count,
                              breaks=breaks)
  all.cnv.numbers$bin <- addNA(all.cnv.numbers$bin)
  levels(all.cnv.numbers$bin) <- c(levels(all.cnv.numbers$bin), 0)
  all.cnv.numbers$bin[is.na(all.cnv.numbers$bin)] <- 0
  all.cnv.numbers$bin <- factor(all.cnv.numbers$bin, levels = c("0", "(1,10]", "(10,50]", "(50,100]", "(100,1e+03]", "(1e+03,Inf]"))
  del.numbers$bin <- cut(del.numbers$Count,
                          breaks=breaks)
  del.numbers$bin <- addNA(del.numbers$bin)
  levels(del.numbers$bin) <- c(levels(del.numbers$bin), 0)
  del.numbers$bin[is.na(del.numbers$bin)] <- 0
  del.numbers$bin <- factor(del.numbers$bin, levels = c("0", "(1,10]", "(10,50]", "(50,100]", "(100,1e+03]", "(1e+03,Inf]"))
  dup.numbers$bin <- cut(dup.numbers$Count,
                          breaks=breaks)
  dup.numbers$bin <- addNA(dup.numbers$bin)
  levels(dup.numbers$bin) <- c(levels(dup.numbers$bin), 0)
  dup.numbers$bin[is.na(dup.numbers$bin)] <- 0
  dup.numbers$bin <- factor(dup.numbers$bin, levels = c("0", "(1,10]", "(10,50]", "(50,100]", "(100,1e+03]", "(1e+03,Inf]"))

  all.cnv.numbers$Caller <- factor(all.cnv.numbers$Caller,
                                  levels = c("PennCNV", "QuantiSNP", "iPattern", "EnsembleCNV", "R-GADA", "Gold Standard"),ordered = TRUE)
  del.numbers$Caller <- factor(del.numbers$Caller,
                                  levels = c("PennCNV", "QuantiSNP", "iPattern", "EnsembleCNV", "R-GADA", "Gold Standard"),ordered = TRUE)
  dup.numbers$Caller <- factor(dup.numbers$Caller,
                                  levels = c("PennCNV", "QuantiSNP", "iPattern", "EnsembleCNV", "R-GADA", "Gold Standard"),ordered = TRUE)

  # Generate figures
  p1 <- data.frame(caller=all.cnv.numbers$Caller,bin=factor(all.cnv.numbers$bin)) %>% 
  count(caller,bin,.drop=FALSE) %>% 
  mutate(n = ifelse(n == 0, NA, n)) %>% 
  ggplot(aes(x=caller,y=bin,fill=n)) + geom_tile() +
  labs(x="", y="Number of called CNVs by single sample") +
  theme_classic() +
  theme(axis.text.x = element_text(angle=90, vjust=0.5, hjust=1), legend.position = "none") +
  scale_y_discrete(labels = tags) +
  scale_fill_viridis(name ="Number of samples", na.value="transparent")

  p2 <- data.frame(caller=del.numbers$Caller,bin=factor(del.numbers$bin)) %>% 
  count(caller,bin,.drop=FALSE) %>%
  mutate(n = ifelse(n == 0, NA, n)) %>% 
  ggplot(aes(x=caller,y=bin,fill=n)) + geom_tile() +
  labs(x="", y="Number of called CNVs by single sample") +
  theme_classic() +
  theme(axis.text.x = element_text(angle=90, vjust=0.5, hjust=1), legend.position = "none") +
  scale_y_discrete(labels = tags) +
  scale_fill_viridis(name ="Number of samples", na.value="transparent")

  p3 <- data.frame(caller=dup.numbers$Caller,bin=factor(dup.numbers$bin)) %>% 
  count(caller,bin,.drop=FALSE) %>% 
  mutate(n = ifelse(n == 0, NA, n)) %>% 
  ggplot(aes(x=caller,y=bin,fill=n)) + geom_tile() +
  labs(x="", y="Number of called CNVs by single sample") +
  theme_classic() +
  theme(axis.text.x = element_text(angle=90, vjust=0.5, hjust=1), legend.position = "none") +
  scale_y_discrete(labels = tags) +
  scale_fill_viridis(name ="Number of samples", na.value="transparent")

  figure <- ggarrange(p1, p2, p3, ncol=3, common.legend = T, legend='right', labels=c("A", "B", "C"), label.y=0.05)
  tiff(filename, units="cm", width=24, height=24, res=600)
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

# generate plots
plotCNVCount(cnv.calls, light_rhg_cols, filename=paste0(opt$output,"_CNV_counts.tiff"))
