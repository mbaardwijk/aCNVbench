nextflow.enable.dsl=2

// Include the CNV detection workflows
include { RunPennCNV } from './Workflows/penncnv'
include { RunQuantiSNP } from './Workflows/quantisnp'
include { RuniPattern } from './Workflows/ipattern'
include { RunRGADA } from './Workflows/rgada'
include { RunEnsembleCNV } from './Workflows/ensemblecnv'

// Retrieve input files and parameters
Channel .fromPath( params.gcModel ).set{ gcModel }
Channel .fromPath( params.hmmFile ).set{ hmmFile }
Channel .fromPath( params.experimentName ).set{ experimentName }
Channel .empty().set{ callFiles }

process SamplesheetToGenderFile {
    /*   Retrieve sample ID and gender columns from input samplesheet
    */
    label 'cpu_one'

    input:
    tuple file(reportFile), file(samplesheetFile), file(snpMap)

    output:
    tuple file(reportFile), file(samplesheetFile), file(snpMap), path('*_gender_file.txt'), emit: inputFiles

    script:
    """
    awk 'NR==1{for (i=1; i<=NF; i++) if (\$i == "ID" ){a=i} else if(\$i == "Gender"){b=i}} {print \$a "\t" \$b}' ${samplesheetFile} > ${params.experimentName}_gender_file.txt
    sed -i 's/ID/Sample_ID/g' ${params.experimentName}_gender_file.txt
    """
}

process SliceInputFiles {
    /*   If a subset of chromosomes is given by the user, this process filters the
    Illumina report file on those chromosomes only
    */
    label 'cpu_one'

    input:
    tuple file(reportFile), file(samplesheetFile), file(snpMap), file(genderFile)

    output:
    tuple file('Filtered_report.txt'), file(samplesheetFile), file('Filtered_snpmap.txt'), file(genderFile), emit: filteredInput

    script:
    """
    dataline=\$(head -n 1000 ${reportFile} | grep -Fn "[Data]" | cut -f1 -d:)

    head -n \$(expr \${dataline} + 1) ${reportFile} > Header_report.txt
    tail -n +\$(expr \${dataline} + 2) ${reportFile} | awk 'BEGIN{split("${params.chrNames}", CHR)} {for(i in CHR) {if(\$3 == CHR[i]) print \$0}}' > Filtered_data_report.txt
    cat Header_report.txt Filtered_data_report.txt > Filtered_report.txt

    head -n 1 ${snpMap} > Header_snpmap.txt
    awk 'BEGIN{split("${params.chrNames}", CHR)} {for(i in CHR) {if(\$2 == CHR[i]) print \$0}}' ${snpMap} > Filtered_data_snpmap.txt
    cat Header_snpmap.txt Filtered_data_snpmap.txt > Filtered_snpmap.txt
    """
}

process SplitIlluminaReport {
    label 'cpu_one'
    label 'penncnv'

    input:
    tuple file(reportFile), file(samplesheetFile), file(snpMap), file(genderFile)

    output:
    tuple file(reportFile), file(samplesheetFile), file(snpMap), file(genderFile), path('*.signal'), emit: signalFile

    script:
    """
    /PennCNV-1.0.5/split_illumina_report.pl -suffix .signal ${reportFile}
    """
}

process SignalToBaflrr {
    /*   Convert signal files to PennCNV input format with BAF and LRR values
    */
    label 'cpu_8x'
    label 'acnvbench'
    
    input:
    tuple file(reportFile), file(samplesheetFile), file(snpMap), file(genderFile), file(signalFiles)

    output:
    tuple file(reportFile), file(samplesheetFile), file(snpMap), file(genderFile), file('*.baflrr'), emit: baflrrFile

    script:
    """
    echo ${snpMap}
    python3 /opt/signal_to_baflrr.py ./ ${snpMap}
    """
}

// Start detection and benchmark workflow
workflow {
    Channel .fromPath( params. inputSheet ) \
        | splitCsv(header:true, sep:'\t') \
        | map { row -> tuple(file(row.reportFile), file(row.samplesheetFile), file(row.snpMap) )}
        | set { inputData }
    SamplesheetToGenderFile( inputData )
    if(params.chrNames){ /* if only a subset of chromosomes is to be analyzed */
        SliceInputFiles( SamplesheetToGenderFile.out.inputFiles )
        SplitIlluminaReport( SliceInputFiles.out.filteredInput )
    }
    else {
        SplitIlluminaReport( SamplesheetToGenderFile.out.inputFiles )
    }
    SignalToBaflrr(      SplitIlluminaReport.out.signalFile )
    if(!params.skipPenncnv) {
        RunPennCNV(      SignalToBaflrr.out.baflrrFile )
    }
    else {
        penncnvCalls = 'NO_FILE'
    }
    if(!params.skipQuantisnp) {
        RunQuantiSNP(    SignalToBaflrr.out.baflrrFile )
    }
    else {
        quantisnpCalls = Channel.of('NO_FILE')
    }
     if(!params.skipIpattern) {
        if(params.chrNames){
            RuniPattern(     SliceInputFiles.out.filteredInput )
                       }
        else {
            RuniPattern(     SamplesheetToGenderFile.out.inputFiles )
        }
    }
    else {
        ipatternCalls = 'NO_FILE'
    }
    if(!params.skipRgada) {
        RunRGADA(        SignalToBaflrr.out.baflrrFile )
    }
    else {
        rgadaCalls = 'NO_FILE'
    }

    if(!params.skipEnsemblecnv) {
        if(params.chrNames){
            inputfiles = SliceInputFiles.out.filteredInput
            penncnvOutput = RunPennCNV.out.resultsDir
            quantisnpOutput = RunQuantiSNP.out.resultsDir
            ipatternOutput = RuniPattern.out.resultsDir
            ensemblecnvInput = inputfiles.join(penncnvOutput).join(quantisnpOutput).join(ipatternOutput)
            RunEnsembleCNV(     SignalToBaflrr.out.baflrrFile,
                                RunPennCNV.out.resultsDir,
                                RunQuantiSNP.out.resultsDir,
                                RuniPattern.out.resultsDir )
                       }
        else {
            inputfiles = SamplesheetToGenderFile.out.inputFiles
            penncnvOutput = RunPennCNV.out.resultsDir
            quantisnpOutput = RunQuantiSNP.out.resultsDir
            ipatternOutput = RuniPattern.out.resultsDir
            ensemblecnvInput = inputfiles.join(penncnvOutput).join(quantisnpOutput).join(ipatternOutput)
            RunEnsembleCNV(     SignalToBaflrr.out.baflrrFile,
                                RunPennCNV.out.resultsDir,
                                RunQuantiSNP.out.resultsDir,
                                RuniPattern.out.resultsDir )
                       }
        }
    else {
        ensemblecnvCalls = 'NO_FILE'
    }
}
