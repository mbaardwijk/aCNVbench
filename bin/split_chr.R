suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(optparse))

# define all necessary variables and files
option_list = list(
  make_option(c("-m", "--samples"), action = "store", type = "character", default = NA,
              help = "Text file listing the names of all samples"),
  make_option(c("-p", "--probes"), action = "store", type = "character", default = NA,
              help = "Tab separated file with sorted probe information, consisting of columns Name, Chr, Postition. No header."),
  make_option(c("-s", "--source"), action = "store", type = "character", default = NA,
              help = "Source directory of the input files"),
  make_option(c("-d", "--destination"), action = "store", type = "character", default = NA,
              help = "Destination directory of the output files"),
  make_option(c("-x", "--suffix"), action = "store", type = "character", default = NA,
              help = "Suffix of the input files"),
  make_option(c("-g", "--gender"), action = "store", type = "character", default = NA,
              help = "Text file storing the gender information of the samples"),
  make_option(c("-q", "--pq"), action = "store", type = "character", default = NA,
              help = "Text file storing the location of the split between p and q arms of the chromosomes"),
  make_option(c("-z", "--experimentname"), action = "store", type = "character", default = NA,
              help = "Experimentname used to name output files")
)

opt = parse_args(OptionParser(option_list = option_list))

samplefile = readLines(con = opt$samples)
probefile = read.table(opt$probes, header=F, col.names=c("Probe_id", "Chr", "Pos"))
chromosomes = unique(probefile$Chr)
# currently only works for autosomal chromosomes
chromosomes = chromosomes[chromosomes != "X" | chromosomes != "Y"] 
chromosome.arms = c("p", "q")
pqfile = read.table(opt$pq, header=F, col.names=c("Chr", "Position"))

sourcedir = opt$source
# sourcedir = paste(paste(getwd(), sourcedir, sep="/"))
pattern = paste0(".", opt$suffix)

# open all individual signal files with given extension
filenames = paste0(sourcedir, "/", samplefile, ".vn")
myfiles = lapply(filenames, fread)
firstfile = myfiles[[1]]

for(chr in chromosomes){
  pq.position = pqfile$Position[pqfile$Chr == chr]
  p.probes = subset(probefile, Chr == chr & Pos < pq.position)
  q.probes = subset(probefile, Chr == chr & Pos > pq.position)
  p.probes.indices = as.numeric(factor(rownames(probefile)[probefile$Chr == chr & probefile$Pos < pq.position]))
  q.probes.indices = as.numeric(factor(rownames(probefile)[probefile$Chr == chr & probefile$Pos > pq.position]))
  for(i in 1:length(myfiles)){
    cur.p.probes = myfiles[[i]][p.probes.indices]
    names(cur.p.probes) = samplefile[i]
    p.probes = cbind(p.probes, cur.p.probes)
    cur.q.probes = myfiles[[i]][q.probes.indices]
    names(cur.q.probes) = samplefile[i]
    q.probes = cbind(q.probes, cur.q.probes)
  }
  if(nrow(p.probes) > 1){
    write.table(p.probes, file=paste0(opt$destination, "/chr", chr, ".p.", opt$experimentname, ".int"), row.names=F, col.names=T, quote=F, sep="\t")
  }
  if(nrow(q.probes) > 1){
    write.table(q.probes, file=paste0(opt$destination, "/chr", chr, ".q.", opt$experimentname, ".int"), row.names=F, col.names=T, quote=F, sep="\t")
  }
}



