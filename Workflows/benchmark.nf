process CompareCalls {
    /*
    */
    container = 'mbaardwijk/bedr:galaga'
    
    input:
    tuple val(sampleID), path(callFile), val(platformName), path(goldstandardFile)

    output:
    tuple val(sampleID), val(platformName), path("${sampleID}.${platformName}.results.txt")
    
    script:
    """
    Rscript ${PWD}/bin/benchmark.performance.R -g ${goldstandardFile} -c ${callFile} -l ${params.chainFile} -o ${sampleID}.${platformName}.results.txt
    """
}

process CombineResults {
    /*
    */
    container 'mbaardwijk/acnvbench:16-11-2023'

    input:
    path resultFiles
    val platformName

    output:
    path "${platformName}.performances.txt"
    
    """
    Rscript ${PWD}/bin/combine.results.R --input '${resultFiles}' --output '${platformName}'
    """
}

process PlotCNVCounts {
    /*
    */
    publishDir "${params.resultsDir}", mode:'copy'
    container 'mbaardwijk/acnvbench:16-11-2023'

    input:
    path callFiles
    path goldstandardVcf

    output:
    path "*_CNV_counts.tiff" 

    """
    Rscript ${PWD}/bin/plot.cnv.counts.R --callfiles '${callFiles}' --goldstandard '${goldstandardVcf}' --output ${params.output}
    """
}

process PlotCNVSize {
    /*
    */
    publishDir "${params.resultsDir}", mode:'copy'
    container 'mbaardwijk/acnvbench:16-11-2023'

    input:
    path callFiles
    path goldstandardVcf

    output:
    path "aCNVBench_CNV_sizes_*.tiff"

    """
    Rscript ${PWD}/bin/plot.cnv.sizes.heatmap.R --callfiles '${callFiles}' --goldstandard '${goldstandardVcf}' --output ${params.output}
    """
}

process PlotUpSet {
    /*
    */
    publishDir "${params.resultsDir}", mode:'copy'
    container 'mbaardwijk/acnvbench:16-11-2023'

    input:
    path callFiles

    output:
    path "Alt_UpSetPlot_overlap.tiff"

    """
    Rscript ${PWD}/bin/plot.upset.overlap.cnv.calls.R --callfiles '${callFiles}' --output ${params.output}
    """
}

process PlotPrecisionRecall {
    /*
    */
    publishDir "${params.resultsDir}", mode:'copy'
    container 'mbaardwijk/acnvbench:16-11-2023'

    input:
    path performanceFiles

    output:
    path "aCNVBench_overall_precision_recall.tiff"

    """
    Rscript ${PWD}/bin/plot.precision.recall.R --performances '${performanceFiles}' --output ${params.output}
    """    
}

workflow {
    Channel .fromPath( params. individualInputSheet ) \
        | splitCsv(header:true, sep:'\t') \
        | map { row -> tuple(row.inputFile, row.sampleID, row.platform)} \
        | set { inputData }
    Channel .fromPath( params. goldstandardSheet ) \
        | splitCsv(header:true, sep:'\t') \
        | map { row -> tuple(row.goldstandardFile, row.sampleID)} \
        | set { goldstandardData }
    combinedData = inputData.join(  goldstandardData, by:[1])
    Channel .fromPath( params. callsetInputSheet ).splitText().map{it -> it.trim()}.collect().set{callsetFiles}
    callsetFiles.view()
    Channel .fromPath( params. goldstandardVcf ) \
        | set { goldstandardVcf }
    CompareCalls(                   combinedData )
    results = CompareCalls.out.branch{  sampleID, platformName, resultsFile ->
                                        penncnv: platformName == "PennCNV"
                                            return resultsFile
                                        quantisnp: platformName == "QuantiSNP"
                                            return resultsFile
                                        ipattern: platformName == "iPattern"
                                            return resultsFile
                                        ensemblecnv: platformName == "EnsembleCNV"
                                            return resultsFile
                                        rgada: platformName == "R-GADA"
                                            return resultsFile
    }
    penncnvFiles =                      results.penncnv.collect(sort: true)
    quantisnpFiles =                    results.quantisnp.collect(sort: true)
    ipatternFiles =                     results.ipattern.collect(sort: true)
    ensemblecnvFiles =                  results.ensemblecnv.collect(sort: true)
    rgadaFiles =                        results.rgada.collect(sort: true)
    resultFiles = penncnvFiles.concat(  quantisnpFiles, ipatternFiles, ensemblecnvFiles, rgadaFiles )
    platformNames = Channel.from(       "PennCNV", "QuantiSNP", "iPattern", "EnsembleCNV", "R-GADA" )
    CombineResults(                     resultFiles,
                                        platformNames )
    // PlotPrecisionRecall(                CombineResults.out.collect() )
    PlotCNVCounts(                      callsetFiles,
                                        goldstandardVcf )
    // PlotCNVSize(                        callsetFiles,
    //                                    goldstandardVcf )
    PlotUpSet(                          callsetFiles )
}
