suppressPackageStartupMessages(library(optparse))
suppressPackageStartupMessages(library(tidyr))

option_list = list(
  make_option(c("-i", "--input"), action = "store", type = "character", default = NA,
              help = "Input csv file with columns Sample_ID, Chr, Start, End, Type"),
  make_option(c("-s", "--separator"), action = "store", type = "character", default = ",",
              help = "Separator for input csv file"),
  make_option(c("-o", "--output"), action = "store", type = "character", default = "results",
              help = "Suffix for outputfiles generated during the splitting")
)

opt = parse_args(OptionParser(option_list = option_list))

callset <- read.csv(opt$input, header=T, sep=opt$separator)
samples <- unique(callset$Sample_ID)

for (sample in samples){
    sample.set <- subset(callset, Sample_ID == sample)
    write.table(sample.set, file=paste0(sample, opt$output), row.names=F, quote=F, sep="\t")
}