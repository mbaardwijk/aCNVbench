# aCNVbench
A Nextflow workflow for detecting and benchmarking CNVs from SNP-array data. The workflow includes 5 different CNV calling methods - PennCNV, QuantiSNP, iPattern, EnsembleCNV and R-GADA. 

## Requirements
- Nextflow version >= 20.xx.xx (tested on v23.04.1.5866)
- Docker (tested on v19.03.15)

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
