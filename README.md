# aCNVbench

![Pipeline Overview](https://github.com/mbaardwijk/aCNVbench/blob/main/Figure_1_workflow_figure.tif)
A Nextflow workflow for detecting and benchmarking CNVs from SNP-array data. The workflow includes 5 different CNV calling methods - PennCNV, QuantiSNP, iPattern, EnsembleCNV and R-GADA. 

## Requirements
- Nextflow version >= 20.xx.xx (tested on v23.04.1.5866)
- Docker (tested on v19.03.15)

## Citation
This workflow was developed for our publication titled 'A systematic benchmark of copy number variation detection tools for high density SNP genotyping arrays'. Please use the following citation:

van Baardwijk M.N., Heijnen L.S.E.M., Zhao H., Baudis M., Stubss A.P. (2024). A systematic benchmark of copy number variation detection tools for high density SNP genotyping arrays. Genomics, (in press). https://doi.org/10.1016/j.ygeno.2024.110962.

## Input data format
The workflow requires an input sheet specifying the input SNP-array data. This input sheet is a tab-separated file that includes the following:
-	The Illumina final report file, containing the columns Sample ID, SNP Name, Chr, Position, Allele 1 – Forward, Allele 2 – Forward, X, Y, B allele frequency and Log R ratio.
-	A sample map, with at least the columns Name and Gender.
-	A probe map, specifying the chromosome and position for every probe on the array.

An example can be found below:

| reportFile | samplesheetFile | snpMap |
|------------|-----------------|--------|
| Data/1000Genomes_GenomeStudio_FinalReport.txt | Data/1000Genomes_Sample_Map.txt | Data/HumanOmni2.5-4v1_B_SNP_Map.txt |

When the workflow is run in benchmarking mode (= default), a gold standard file specifying the gold standard CNVs is also expected. This file can be a VCF file or tab-separated file with the columns Sample_ID, Chr, Start, End and Type.

### Illumina final report file format
An example showing the header of an Illumina final report file can be found below:

> [Header]
>
> GSGT Version &emsp; 2.0.5
>
> Processing Date 4/20/2023 5:27 PM
>
> Content &emsp; HumanOmni2.5-4v1-Multi_B.bpm
>
> Num SNPs &emsp; 2450000
>
> Total SNPs &emsp; 2450000
>
> Num Samples &emsp; 2141
>
> Total Samples &emsp; 2141
>
> [Data]
> Sample ID &emsp; SNP Name &emsp; Chr &emsp; Position &emsp; Allele1 - Forward &emsp; Allele2 - Forward &emsp; X &emsp; Y &emsp; Log R Ratio &emsp; B Allele Freq
>
> 5434203053_R01C01 &emsp; GA008510 &emsp; Y &emsp; 11771305 &emsp; G &emsp; G &emsp; 0.052 &emsp; 0.631 &emsp; 0.5219 &emsp; 1.0000
>
> 5434203053_R01C01 &emsp; GA008524 &emsp; Y &emsp; 19612089 &emsp; T &emsp; T &emsp; 0.658 &emsp; 0.030 &emsp; 0.2345 &emsp; 0.0000

### Sample map file format
An example showing a sample map file can be found below:

> Index &emsp; Name &emsp; ID &emsp; Gender &emsp; Plate &emsp; Well &emsp; Group &emsp; Parent1 &emsp; Parent2 &emsp; Replicate &emsp; SentrixPosition
>
> 1 &emsp; 359146 &emsp; 5434203053_R01C01 &emsp; Female &emsp; WG0077073-AMP2 &emsp; A09 &emsp; &emsp; &emsp; &emsp; 5434203053_R01C01
>
> 2 &emsp; 359130 &emsp; 5434203053_R02C01 &emsp; Female &emsp; WG0077073-AMP2 &emsp; C09 &emsp; &emsp; &emsp; &emsp; 5434203053_R02C01

### Probe map file format
An example showing a probe map file can be found below:

> Name  Chromosome &emsp; Position
>
> GA008510 &emsp; Y &emsp; 11771305
> 
> GA008524 &emsp; Y &emsp; 19612089
>
> GA008529 &emsp; Y &emsp; 19613392
>
> GA008532 &emsp; Y &emsp; 26872641

## Run the workflow
An example of how to run the workflow for detecting and benchmarking CNVs:
```
nextflow run main.nf --inputSheet {path to tab-separated SNP array input file} --goldstandardFile {path to gold standard file in VCF or tab-separated format} –-genome ‘hg18’
```
The benchmarking step kan be skipped using the '--skipBenchmark' parameter. Individual CNV calling methods can also be skipped, for example the '--skipPennCNV can be added to remove PennCNV from the pipeline execution.

## Docker images
All methods were installed in separate docker images. While running the workflow automatically pulls the Docker images from DockerHub if they are not found locally, they can also be retrieved as follows:
```
docker pull mbaardwijk/acnvbench:16-11-2023
docker pull mbaardwijk/penncnv:18-10-2023
docker pull mbaardwijk/quantisnp:19-10-2023
docker pull mbaardwijk/ipattern:16-11-2023
docker pull mbaardwijk/ensemblecnv:27-03-2024
docker pull mbaardwijk/rgada:19-10-2023
```

## Collaborators
- Myrthe van Baardwijk (Erasmus Medical Center Rotterdam)
- Andrew Stubbs (Erasmus Medical Center Rotterdam)
- Hangjia Zhao (University of Zurich)
- Michael Baudis (University of Zurich)
