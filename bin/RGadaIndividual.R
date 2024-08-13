library(gada)
library(dplyr)
library(optparse)

option_list = list(
  make_option(c("-i", "--input"), action = "store", type = "character", default = NA,
              help = "Input .baflrr file"),
  make_option(c("-o", "--output"), action = "store", type = "character", default = "results",
              help = "Name of the output file with cnv calls"),
  make_option(c("-t", "--t_statistic"), action = "store", type = "numeric", default = 4,
              help = "Threshold for the t-statistic"),
  make_option(c("-a", "--a_alpha"), action = "store", type = "numeric", default = 0.8,
              help = "A alpha"),
  make_option(c("-m", "--min_seg_length"), action = "store", type = "numeric", default = 100,
              help = "Minimum segment length")
)

opt = parse_args(OptionParser(option_list = option_list))

# get sample name from file name
sample <- gsub(".baflrr", "", opt$input)
sample <- gsub("/", "", opt$input)

cnv.call <- setupGADA(opt$input,
                      log2ratioCol = 4,
                      BAFcol=5)

cnv.call <- SBL(cnv.call,
                estim.sigma2=TRUE,
                aAlpha=opt$a_alpha,
                verbose=T)
cnv.call <- BackwardElimination(cnv.call,
                                T=opt$t_statistic,
                                MinSegLen=opt$min_seg_length)

cnvs <- summary(cnv.call, length.base=c(500,10e6))

cnvs <- subset(cnvs, State != 0)
cnvs$State[cnvs$State == 1] <- "Gain"
cnvs$State[cnvs$State == -1] <- "Loss"
cnvs$ID <- sample

write.table(cnvs, file=opt$output, sep="\t", quote=F, row.names=F)