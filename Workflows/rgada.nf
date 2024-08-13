process GetSampleNamesRgada {
    /*   Get all individual sample names to allow for parallel processing of samples
    */
    input:
    tuple file(reportFile), file(samplesheetFile), file(snpMap), file(genderFile), file(baflrrFiles)

    output:
    path 'samples.txt', emit: sampleFile

    script:
    """
    awk 'NR==1{for (i=1; i<=NF; i++) if (\$i == "ID" ){a=i}} {print \$a}' ${samplesheetFile} > 'samples.txt'
    grep -v "ID" samples.txt > tmpfile && mv tmpfile samples.txt
    """
}

process RGada {
    /*   Execute RGada.R script
    */
    label 'rgada'
    publishDir "${params.resultsDir}/R-GADA", mode:'copy'

    input:
    tuple file(reportFile), file(samplesheetFile), file(snpMap), file(genderFile), file(baflrrFiles)
    each sample

    output:
    path '*_rgada.tsv', emit: cnvFile

    script:
    """
    sampleName=`(echo "${sample}"|tr -d '\n')`
    
    Rscript ${PWD}/bin/RGadaIndividual.R \
    --input "\${sampleName}.baflrr" \
    --output "\${sampleName}_rgada.tsv" \
    --t_statistic ${params.tStatistic} \
    --a_alpha ${params.aAlpha} \
    --min_seg_length ${params.minSegLength} >& Rout.txt
    """
}

process ConvertRGadaCNVs  {
    /*  Convert cnvs to general benchmark format
    */
    label 'acnvbench'
    publishDir "${params.resultsDir}/CNVs/R-GADA", mode:'copy'

    input:
    file cnvFiles

    output:
    path 'results.rgada.txt', emit: callFile

    script:
    """
    files=(`ls -d *_rgada.tsv`)
    head -n 1 "\${files[0]}" > all.merged.rgada.tsv

    for file in \${files[*]}; do
        tail -n +2 \${file} >> all.merged.rgada.tsv
    done

    Rscript ${PWD}/bin/convert_cnv_results.R -r all.merged.rgada.tsv
    """
}

workflow RunRGADA {
    take:
        baflrrFiles
    main:
        GetSampleNamesRgada( baflrrFiles )
        samples = GetSampleNamesRgada.out.sampleFile.splitText()
        RGada(               baflrrFiles,
                             samples )
        ConvertRGadaCNVs(    RGada.out.cnvFile.collect() )
    emit:
        callFile = ConvertRGadaCNVs.out.callFile
}
