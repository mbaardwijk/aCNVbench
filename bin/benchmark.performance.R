suppressPackageStartupMessages(library(bedr))
suppressPackageStartupMessages(library(GenomicRanges))
suppressPackageStartupMessages(library(liftOver))
suppressPackageStartupMessages(library(optparse))
suppressPackageStartupMessages(library(tidyr))

CALLER.COLUMNS <- c('Sample_ID', 'Chr', 'Start', 'End', 'Type')
CALLER.GR.COLUMNS <- c('Sample_ID', 'seqnames', 'start', 'end', 'Type')
POSITION.COLUMNS.GOLDSTANDARD <- c('Chr', 'Start', 'End')
POSITION.COLUMNS.CALLER <- c('Chr', 'Start', 'End')
AUTOSOMAL.CHRS <- paste0("chr", seq(1:22))
OVERLAP.FRACTIONS <- c(0.5, 0.4, 0.3, 0.2, 0.1, 0.01, 0)
#OVERLAP.FRACTIONS <- c(0)
#LENGTH.LARGE.CNVS <- 100000

option_list <- list(
  make_option(c("-g", "--goldstandard"), action = "store", type = "character", default = NA,
              help = "Path to the vcf file with gold standard cnvs."),
  make_option(c("-c", "--caller"), action = "store", type = "character", default = NA,
              help = "Path to the file with the called cnvs."),
  make_option(c("-l", "--liftoverchain"), action = "store", type = "character", default = "NA",
              help = "Path to the liftOver chainfile."),
  make_option(c("-o", "--output"), action = "store", type = "character", default = NA,
              help = "Name of the output file.")
)

opt <- parse_args(OptionParser(option_list = option_list))

determineOverlap <- function(cnv.calls.gr, goldstandard.gr) {
  if(length(cnv.calls.gr) > 0){
    # Create columns for 
    cnv.calls.gr$Overlap.calls <- 0
    cnv.calls.gr$Overlap.gs <- 0
    cnv.calls.gr$Overlap.set.calls <- 0
    cnv.calls.gr$Overlap.set.gs <- 0
    goldstandard.gr$Overlap.calls <- 0
    goldstandard.gr$Overlap.gs <- 0
    goldstandard.gr$Overlap.set.calls <- 0
    goldstandard.gr$Overlap.set.gs <- 0
    hits <- findOverlaps(cnv.calls.gr, goldstandard.gr)
    overlaps <- suppressWarnings(pintersect(cnv.calls.gr[queryHits(hits)], goldstandard.gr[subjectHits(hits)]))
    # Determine individual reciprocal overlap
    cnv.calls.gr$Overlap.calls[queryHits(hits)] <- width(overlaps) / width(cnv.calls.gr[queryHits(hits)])
    goldstandard.gr$Overlap.calls[subjectHits(hits)] <- width(overlaps) / width(cnv.calls.gr[queryHits(hits)])
    cnv.calls.gr$Overlap.gs[queryHits(hits)] <- width(overlaps) / width(goldstandard.gr[subjectHits(hits)])
    goldstandard.gr$Overlap.gs[subjectHits(hits)] <- width(overlaps) / width(goldstandard.gr[subjectHits(hits)])
    # Determine if there are sets of CNV calls overlapping individual gold standard CNVs
    call.sets.i <- unique(subjectHits(hits)[duplicated(subjectHits(hits), fromLast=TRUE)])
    for(i in call.sets.i) {
      call.set.overlap <- suppressWarnings(intersect(goldstandard.gr[i], cnv.calls.gr))
      set <- subjectHits(suppressWarnings(findOverlaps(goldstandard.gr[i], cnv.calls.gr)))
      percent.overlap.calls <- sum(width(call.set.overlap) / width(goldstandard.gr[i]))
      cnv.calls.gr$Overlap.set.calls[set] <- percent.overlap.calls
      goldstandard.gr$Overlap.set.calls[i] <- percent.overlap.calls
    }
    # Determine if there are individual CNV calls overlapping sets of gold standard CNVs
    gs.sets.i <- unique(queryHits(hits)[duplicated(queryHits(hits), fromLast=TRUE)])
    for(i in gs.sets.i) {
      gs.set.overlap <- suppressWarnings(intersect(cnv.calls.gr[i], goldstandard.gr))
      set <- subjectHits(suppressWarnings(findOverlaps(cnv.calls.gr[i], goldstandard.gr)))
      percent.overlap.gs <- sum(width(gs.set.overlap) / width(cnv.calls.gr[i]))
      cnv.calls.gr$Overlap.set.gs[i] <- percent.overlap.gs
      goldstandard.gr$Overlap.set.gs[set] <- percent.overlap.gs
    }
  }
  return(list(cnv.calls.gr, goldstandard.gr))
}

classifyCNVs <- function(cnv.calls.gr, goldstandard.gr, min.overlap){
  # All goldstandard CNVs not classified as TP will be classified as FN
  print("hier")
  goldstandard.gr$Class <- "FN"
  print("hier all?")
  if(length(cnv.calls.gr) > 0){
    # All called CNVs not classified as TP will be classified as FP
    cnv.calls.gr$Class <- "FP"
    
    # If min.overlap is 0 we have to use > instead of >=
    if(min.overlap == 0){
      # If any individual called CNVs overlap enough with individual CNVs in the gold standard
      if(any(cnv.calls.gr$Overlap.calls > min.overlap & cnv.calls.gr$Overlap.gs > min.overlap)){
        cnv.calls.gr[cnv.calls.gr$Overlap.calls > min.overlap & cnv.calls.gr$Overlap.gs > min.overlap]$Class <- "TP"
      }
      if(any(goldstandard.gr$Overlap.calls > min.overlap & goldstandard.gr$Overlap.gs > min.overlap)){
        goldstandard.gr[goldstandard.gr$Overlap.calls > min.overlap & goldstandard.gr$Overlap.gs > min.overlap]$Class <- "TP"
      }
      # If any called CNV sets overlap enough with an individual CNV in the gold standard
      if(any(cnv.calls.gr$Overlap.calls > min.overlap & cnv.calls.gr$Overlap.set.calls > min.overlap)){
        cnv.calls.gr[cnv.calls.gr$Overlap.calls > min.overlap & cnv.calls.gr$Overlap.set.calls > min.overlap]$Class <- "TP"
      }
      if(any(goldstandard.gr$Overlap.calls > min.overlap & goldstandard.gr$Overlap.set.calls > min.overlap)){
        goldstandard.gr[goldstandard.gr$Overlap.calls > min.overlap & goldstandard.gr$Overlap.set.calls > min.overlap]$Class <- "TP"
      }
      # If any individual called CNVs overlaps enough with a set of CNVs in the gold standard
      if(any(cnv.calls.gr$Overlap.gs > min.overlap & cnv.calls.gr$Overlap.set.gs > min.overlap)){
        cnv.calls.gr[cnv.calls.gr$Overlap.gs > min.overlap & cnv.calls.gr$Overlap.set.gs > min.overlap]$Class <- "TP"
      }
      if(any(goldstandard.gr$Overlap.gs > min.overlap & goldstandard.gr$Overlap.set.gs > min.overlap)){
        goldstandard.gr[goldstandard.gr$Overlap.gs > min.overlap & goldstandard.gr$Overlap.set.gs > min.overlap]$Class <- "TP"
      }
    } else {
      # If any individual called CNVs overlap enough with individual CNVs in the gold standard
      if(any(cnv.calls.gr$Overlap.calls >= min.overlap & cnv.calls.gr$Overlap.gs >= min.overlap)){
        cnv.calls.gr[cnv.calls.gr$Overlap.calls >= min.overlap & cnv.calls.gr$Overlap.gs >= min.overlap]$Class <- "TP"
      }
      if(any(goldstandard.gr$Overlap.calls >= min.overlap & goldstandard.gr$Overlap.gs >= min.overlap)){
        goldstandard.gr[goldstandard.gr$Overlap.calls >= min.overlap & goldstandard.gr$Overlap.gs >= min.overlap]$Class <- "TP"
      }
      # If any called CNV sets overlap enough with an individual CNV in the gold standard
      if(any(cnv.calls.gr$Overlap.calls >= min.overlap & cnv.calls.gr$Overlap.set.calls >= min.overlap)){
        cnv.calls.gr[cnv.calls.gr$Overlap.calls >= min.overlap & cnv.calls.gr$Overlap.set.calls >= min.overlap]$Class <- "TP"
      }
      if(any(goldstandard.gr$Overlap.calls >= min.overlap & goldstandard.gr$Overlap.set.calls >= min.overlap)){
        goldstandard.gr[goldstandard.gr$Overlap.calls >= min.overlap & goldstandard.gr$Overlap.set.calls >= min.overlap]$Class <- "TP"
      }
      # If any individual called CNVs overlaps enough with a set of CNVs in the gold standard
      if(any(cnv.calls.gr$Overlap.gs >= min.overlap & cnv.calls.gr$Overlap.set.gs >= min.overlap)){
        cnv.calls.gr[cnv.calls.gr$Overlap.gs >= min.overlap & cnv.calls.gr$Overlap.set.gs >= min.overlap]$Class <- "TP"
      }
      if(any(goldstandard.gr$Overlap.gs >= min.overlap & goldstandard.gr$Overlap.set.gs >= min.overlap)){
        goldstandard.gr[goldstandard.gr$Overlap.gs >= min.overlap & goldstandard.gr$Overlap.set.gs >= min.overlap]$Class <- "TP"
      }
    }
  }
  # Combine TP, FP and FN calls
  fns <- goldstandard.gr[goldstandard.gr$Class == "FN"]
  cnv.results <- suppressWarnings(pc(c(cnv.calls.gr, fns)))
  return(cnv.results)
}

scoreSample <- function(cnv.results){
  tp <- sum(cnv.results$Class == "TP")
  fp <- sum(cnv.results$Class == "FP")
  fn <- sum(cnv.results$Class == "FN")
  return(list(Sample_ID = cnv.results$Sample_ID[1], TP = tp, FP = fp, FN = fn))
}

test.cnvs <- function(cnv.calls.gr, goldstandard.gr) {
  # Find overlaps
  cnv.overlaps <- determineOverlap(cnv.calls.gr, goldstandard.gr)
  cnv.calls.gr <- cnv.overlaps[[1]]
  goldstandard.gr <- cnv.overlaps[[2]]
  # Determine TPs FPs and FNs for each overlap fraction
  all.cnv.scores <- data.frame(matrix(ncol=5,nrow=0, dimnames=list(NULL, c("Sample_ID", "TP", "FP", "FN", "Overlap.t"))))
  for(min.overlap in OVERLAP.FRACTIONS){
    print(min.overlap)
    cnv.results <- classifyCNVs(cnv.calls.gr, goldstandard.gr, min.overlap)
    print(cnv.results)
    # Determine scores
    cnv.scores <- scoreSample(cnv.results)
    print(cnv.scores)
    cnv.scores$Overlap.t <- min.overlap
    all.cnv.scores <- rbind(all.cnv.scores, cnv.scores)
  }
  return(all.cnv.scores)
}

# determine individual CNVs based on ALT and GT fields from VCF
parseCNVGenotype <- function(vcf.cnv.data, gt.column.name){
  vcf.cnv.data$Sample_ID <- gt.column.name
  vcf.cnv.data <- separate(data = vcf.cnv.data, col = gt.column.name, into = c("GT.1", "GT.2"))
  vcf.cnv.data$Type <- "REF"
  vcf.cnv.data$Type[vcf.cnv.data$GT.1 == 1 | vcf.cnv.data$GT.2 == 1] <- vcf.cnv.data$ALT[vcf.cnv.data$GT.1 == 1 | vcf.cnv.data$GT.2 == 1]
  vcf.cnv.data$Type <- gsub("<", "", vcf.cnv.data$Type)
  vcf.cnv.data$Type <- gsub(">", "", vcf.cnv.data$Type)
  return(vcf.cnv.data)
}

# change genome build of called CNVs if genome differs from goldstandard
liftover.cnvs <- function(caller, liftover.chain) {
  chain <- import.chain(liftover.chain)
  ranges <- makeGRangesFromDataFrame(caller, keep.extra.columns = T, na.rm=T)
  seqlevelsStyle(ranges) <- "UCSC"
  rangesLiftOver <- liftOver(ranges, chain)
  dataLiftOver <- as.data.frame(rangesLiftOver)[, CALLER.GR.COLUMNS]
  colnames(dataLiftOver) <- CALLER.COLUMNS
  return(dataLiftOver)
}

caller.filename <- opt$caller
goldstandard.filename <- opt$goldstandard
liftover.chain.filename <- opt$liftoverchain

# read the gold standard and caller files
caller <- read.table(caller.filename, sep="\t", header=TRUE, stringsAsFactors = FALSE)
goldstandard <- read.table(goldstandard.filename, sep="\t", header=TRUE, stringsAsFactors = FALSE)
sample.id <- colnames(goldstandard)[ncol(goldstandard)]

# print(head(caller))
# print(head(goldstandard))

# convert to chr notation with chr if not currently used and select autosomal chrs only
if(! any(grepl("chr", goldstandard$Chr))){
  goldstandard$Chr <- paste0("chr", goldstandard$Chr)
}
if(! any(grepl("chr", caller$Chr))){
  caller$Chr <- paste0("chr", caller$Chr)
}
goldstandard <- subset(goldstandard, Chr %in% AUTOSOMAL.CHRS)
#print(head(goldstandard))

caller <- subset(caller, Chr %in% AUTOSOMAL.CHRS)
#print(head(caller))

# if any autosomal CNVs are found
if(nrow(caller != 0)){
  # liftover if liftover chain is provided
  if(liftover.chain.filename != "NA") {
    caller <- liftover.cnvs(caller, liftover.chain.filename)
  }
  # retrieve deletions and duplications from genotype
  #goldstandard <- parseCNVGenotype(goldstandard, sample.id)
  caller.gr <- makeGRangesFromDataFrame(caller, keep.extra.columns=T)
  goldstandard.gr <- makeGRangesFromDataFrame(goldstandard, keep.extra.columns=T)
  # separate deletions and duplications
  goldstandard.dels <- goldstandard.gr[goldstandard.gr$Type == "DEL",]
  goldstandard.dups <- goldstandard.gr[goldstandard.gr$Type == "DUP",]
  caller.dels <- caller.gr[caller.gr$Type == "Loss",]
  caller.dups <- caller.gr[caller.gr$Type == "Gain",]
  
  print("hier")
  print(head(caller.gr))
  print(head(caller.dels))
  print(head(goldstandard.dels))
  
  # determine tp, fp, and fn for all catagories
  result.dels <- test.cnvs(caller.dels, goldstandard.dels)
  result.dups <- test.cnvs(caller.dups, goldstandard.dups)
  result.all <- as.data.frame(Map("+", result.dels[, c("TP", "FP", "FN")], result.dups[,c("TP", "FP", "FN")]))
  result.all$Sample_ID <- result.dels$Sample_ID
  result.all$Overlap.t <- result.dels$Overlap.t

  # combine all catagories in a single dataframe and calculate sensitivity, precision and f1
  print(head(result.dels))
  print(head(result.dups))
  print(head(result.all))
  results <- rbind(result.dels, result.dups, result.all)

  if(nrow(results) > 0){
    results$Category <- rep(c("Deletions", "Duplications","All"), each=length(OVERLAP.FRACTIONS))
    results$Sensitivity = results$TP / (results$TP + results$FN)
    results$Precision = results$TP / (results$TP + results$FP)
    results$F1.score = 2*((results$Precision*results$Sensitivity)/(results$Precision+results$Sensitivity))
  }
  print(head(results))
 
} else { # if no autosomal cnvs are found
  results <- as.data.frame(matrix(ncol=9, nrow=0))
  colnames(results) <- c("Sample_ID", "TP", "FP", "FN", "Overlap.t", "Category", "Sensitivity", "Precision", "F1.score")
}

# write to output
write.table(results, file=opt$output, quote=F, sep="\t", col.names=T, row.names=F)
