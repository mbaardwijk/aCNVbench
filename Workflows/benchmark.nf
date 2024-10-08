process SplitGoldstandardData {
    /*
    */
    container 'mbaardwijk/acnvbench:16-11-2023'

    input:
    path goldstandardFile

    output:
    path "*_goldstandard.txt", emit: goldstandardFiles
    env samples, emit: sampleIDs

    """
    Rscript ${PWD}/bin/split_callset.R --input '${goldstandardFile}' --output '_goldstandard.txt'
    samples=`(ls *_goldstandard.txt | sed "s/_goldstandard.txt//g")`
    echo \$samples
    """
}

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
    Rscript ${PWD}/bin/plot.cnv.counts.heatmap.R --callfiles '${callFiles}' --goldstandard '${goldstandardVcf}' --output ${params.output}
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
    path "*_CNV_sizes_unique_log.tiff"

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
    path "UpSetPlot_overlap.tiff"

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
    path "aCNVbench_overall_precision_recall.tiff"

    """
    Rscript ${PWD}/bin/plot.precision.recall.R --performances '${performanceFiles}' --output ${params.output}
    """    
}

workflow RunBenchmark {
    take:
        callFiles
        individualCallFiles

    main:
        Channel .fromPath( params. goldstandardFile ) \
            | set { goldstandardData }
        
        SplitGoldstandardData(              goldstandardData )
        sampleIDs = SplitGoldstandardData.out.sampleIDs.splitCsv( sep: " " ) \
            | flatten()
        SplitGoldstandardData.out.goldstandardFiles.flatten() \
            | merge( sampleIDs ) \
            | set { goldstandardFiles}

        combinedData = individualCallFiles.combine(  goldstandardFiles, by: 1 )
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
        // resultFiles = penncnvFiles.concat(  quantisnpFiles, rgadaFiles )

        platformNames = Channel.from(       "PennCNV", "QuantiSNP", "R-GADA" )
        CombineResults(                     resultFiles,
                                            platformNames )                                 
        PlotPrecisionRecall(                CombineResults.out.collect() )
        PlotCNVCounts(                      callFiles.collect(),
                                            goldstandardData )
        PlotCNVSize(                        callFiles.collect(),
                                            goldstandardData )
        callFiles.collect().view()
        PlotUpSet(                          callFiles.collect() )
}
