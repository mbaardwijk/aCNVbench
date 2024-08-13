suppressPackageStartupMessages(library(optparse))
suppressPackageStartupMessages(library(tidyr))

option_list = list(
  make_option(c("-p", "--penncnv"), action = "store", type = "character", default = NA,
              help = "PennCNV data with columns ID, CHR, start, end"),
  make_option(c("-q", "--quantisnp"), action = "store", type = "character", default = NA,
              help = "QuantiSNP data with columns ID, CHR, start, end"),
  make_option(c("-i", "--ipattern"), action = "store", type = "character", default = NA,
              help = "iPattern data with columns ID, CHR, start, end"),
  make_option(c("-r", "--rgada"), action = "store", type = "character", default = NA,
              help = "R-GADA data"),
  make_option(c("--ensemblecnv_samplecnv"), action = "store", type = "character", default = NA,
              help = "Matrix_CN_after_GQ.rds file from EnsembleCNV with CNV_IDs as rows and samples as columns"),
  make_option(c("--ensemblecnv_matrixcn"), action = "store", type = "character", default = NA,
              help = "cnvr-after_GQ file from EnsembleCNV with all properties of individual CNVs"),
  make_option(c("-o", "--outfiles"), action = "store", type = "character", default = "results",
              help = "Prefix for outputfiles generated during the analysis")
)

opt = parse_args(OptionParser(option_list = option_list))

convert.penncnv <- function(penncnv.file, outfile) {
  penncnv.data <- read.table(penncnv.file, sep="", header=F)
  colnames(penncnv.data) <- c("CNV.region","Numsnp", "Length", "CNV.type", "Sample_ID", "startsnp", "endsnp", "conf")
  penncnv.data <- penncnv.data %>% separate("CNV.region", c("Chr", "Start", "End"), sep = "(:|-)")
  penncnv.data <- penncnv.data %>% separate("CNV.type", c("State", "Type"), sep = ",")
  penncnv.data$Type <- gsub("cn=", "", penncnv.data$Type)
  penncnv.data$Sample_ID <- gsub(".baflrr", "", penncnv.data$Sample_ID)
  penncnv.data <- penncnv.data[c("Sample_ID", "Chr", "Start", "End", "Type")]
  penncnv.data$Type <- ifelse(penncnv.data$Type > 2, "DUP", ifelse(penncnv.data$Type == 2, "REF", "DEL"))
  # change chromosomes to 'chr' notation if not current
  penncnv.data$Chr <- ifelse(grepl("chr", penncnv.data$Chr), penncnv.data$Chr, paste0("chr", penncnv.data$Chr))
  # write to output file
  write.table(penncnv.data, paste0(outfile, ".penncnv.txt"), sep="\t", row.names = F, col.names = T, quote=F)
}

convert.quantisnp <- function(quantisnp.file, outfile) {
  quantisnp.data <- read.table(quantisnp.file, sep="\t", header=T)
  quantisnp.data <- quantisnp.data[c("Sample.Name", "Chromosome", "Start.Position..bp.", "End.Position..bp.", "Copy.Number")]
  colnames(quantisnp.data) <- c("Sample_ID", "Chr", "Start", "End", "Type")
  quantisnp.data[, c("Start", "End")] <- sapply(quantisnp.data[, c("Start", "End")], as.numeric)
  quantisnp.data$Sample_ID <- gsub(".baflrr", "", quantisnp.data$Sample_ID)
  quantisnp.data$Type <- ifelse(quantisnp.data$Type > 2, "DUP", ifelse(quantisnp.data$Type == 2, "REF", "DEL"))
  quantisnp.data$Length <- quantisnp.data$End - quantisnp.data$Start
  # change chromosomes to 'chr' notation if not current
  quantisnp.data$Chr <- ifelse(grepl("chr", quantisnp.data$Chr), quantisnp.data$Chr, paste0("chr", quantisnp.data$Chr))
  # write to output file  
  write.table(quantisnp.data, paste0(outfile, ".quantisnp.txt"), sep="\t", row.names = F, col.names = T, quote=F)
}

convert.ipattern <- function(ipattern.file, outfile) {
  ipattern.data <- read.table(ipattern.file, sep="\t", header=F, skip=17)
  colnames(ipattern.data) <- c("Type", "Chr", "Start", "End", "Probe.n",  "On.probe.n", "clusterIdx", "gain_/loss_score", "cluster_score", "gain_/loss_sample.n", "sample_score", "Sample_ID", "CNV_event_ID", "CNVR_ID")
  ipattern.data <- ipattern.data[c("Sample_ID", "Chr", "Start", "End", "Type")]
  ipattern.data[, c("Start", "End")] <- sapply(ipattern.data[, c("Start", "End")], as.numeric)
  ipattern.data$Length <- ipattern.data$End - ipattern.data$Start
  ipattern.data[ipattern.data$Type %in% c("cxGain", "Gain"),]$Type <- "DUP"
  ipattern.data[ipattern.data$Type %in% c("cxLoss", "Loss"),]$Type <- "DEL"
  # change chromosomes to 'chr' notation if not current
  ipattern.data$Chr <- ifelse(grepl("chr", ipattern.data$Chr), ipattern.data$Chr, paste0("chr", ipattern.data$Chr))
  # write to output file  
  write.table(ipattern.data, paste0(outfile, ".ipattern.txt"), sep="\t", row.names = F, col.names = T, quote=F)
}

convert.rgada <- function(rgada.file, outfile) {
  rgada.data <- read.table(rgada.file, sep="\t", header=T)
  rgada.data <- rgada.data[c("ID", "chromosome", "IniProbe", "EndProbe", "State")]
  colnames(rgada.data) <- c("Sample_ID", "Chr", "Start", "End", "Type")
  rgada.data[, c("Start", "End")] <- sapply(rgada.data[, c("Start", "End")], as.numeric)
  rgada.data$Length <- rgada.data$End - rgada.data$Start
  rgada.data[rgada.data$Type == "Gain",]$Type <- "DUP"
  rgada.data[rgada.data$Type == "Loss",]$Type <- "DEL"
  # change chromosomes to 'chr' notation if not current
  rgada.data$Chr <- ifelse(grepl("chr", rgada.data$Chr), rgada.data$Chr, paste0("chr", rgada.data$Chr))
  # write to output file  
  write.table(rgada.data, paste0(outfile, ".rgada.txt"), sep="\t", row.names = F, col.names = T, quote=F)
}

convert.ensemblecnv <- function(ensemblecnv.sample.file, ensemblecnv.cnv.file, outfile) {
  # load files resulting from EnsembleCNV analysis
  ensemblecnv.sample.data <- as.data.frame(readRDS(ensemblecnv.sample.file))
  ensemblecnv.cnv.data <- read.table(ensemblecnv.cnv.file, sep="\t", header=T)
  rownames(ensemblecnv.cnv.data) <- ensemblecnv.cnv.data$CNVR_ID
  
  # for each individual sample, determine cnvs and add to results dataframe
  samples <- colnames(ensemblecnv.sample.data)
  ensemble.results <- data.frame(Sample_ID=character(0), Chr=numeric(0), Start=numeric(0), End=numeric(0), Type=numeric(0))
  for(sample in samples){
    print(sample)
    cnvs <- ensemblecnv.sample.data[ensemblecnv.sample.data[,sample] %in% c(0, 1, 3), sample, drop=FALSE]
    cnv.data <- ensemblecnv.cnv.data[rownames(cnvs),]
    results <- data.frame(Sample_ID=character(nrow(cnvs)), Chr=numeric(nrow(cnvs)), Start=numeric(nrow(cnvs)), End=numeric(nrow(cnvs)), Type=numeric(nrow(cnvs)))
    results$Sample_ID <- sample
    results$Chr <- cnv.data$chr
    results$Start <- cnv.data$posStart
    results$End <- cnv.data$posEnd
    results$Type <- cnvs[,sample]
    ensemble.results <- rbind(ensemble.results, results)
  }
  ensemble.results[, c("Start", "End")] <- sapply(ensemble.results[, c("Start", "End")], as.numeric)
  ensemble.results$Length <- ensemble.results$End - ensemble.results$Start
  ensemble.results$Type <- ifelse(ensemble.results$Type > 2, "DUP", ifelse(ensemble.results$Type == 2, "REF", "DEL"))
  # change chromosomes to 'chr' notation if not current
  ensemble.results$Chr <- ifelse(grepl("chr", ensemble.results$Chr), ensemble.results$Chr, paste0("chr", ensemble.results$Chr))
  # write to output file
  write.table(ensemble.results, paste0(outfile, ".ensemblecnv.txt"), sep="\t", row.names = F, col.names = T, quote=F)
}

if (!is.na(opt$penncnv)){
  convert.penncnv(opt$penncnv, opt$outfile)
}

if (!is.na(opt$quantisnp)){
  convert.quantisnp(opt$quantisnp, opt$outfile)
}

if (!is.na(opt$ipattern)){
  convert.ipattern(opt$ipattern, opt$outfile)
}

if (!is.na(opt$rgada)){
  convert.rgada(opt$rgada, opt$outfile)
}
 
if (!is.na(opt$ensemblecnv_samplecnv)){
  convert.ensemblecnv(opt$ensemblecnv_samplecnv, opt$ensemblecnv_matrixcn, opt$outfile)
}
