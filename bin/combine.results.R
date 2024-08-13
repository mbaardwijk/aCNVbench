suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(optparse))

option_list <- list(
  make_option(c("-i", "--input"), action = "store", type = "character", default = NA,
              help = "List of input files"),
  make_option(c("-o", "--output"), action = "store", type = "character", default = NA,
              help = "Prefix of output files")
)

opt <- parse_args(OptionParser(option_list = option_list))

input.files <- unlist(strsplit(opt$input, "\\s+"))

combined.input <- as.data.frame(matrix(ncol=9, nrow=0))
colnames(combined.input) <- c("Sample_ID", "TP", "FP", "FN", "Overlap.t", "Category", "Sensitivity", "Precision", "F1.score")
# combine individual sample files for different catagories
for(file in input.files){
  print(file)
  results <- read.table(file, sep="\t", header=T)
  #results <- subset(results, select=colnames(combined.input))
  combined.input <- rbind(combined.input, results)
}

# write to output
write.table(combined.input, file=paste0(opt$output, ".performances.txt"), quote=F, sep="\t", col.names=T, row.names = F)
