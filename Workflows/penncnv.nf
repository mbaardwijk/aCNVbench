process CompilePfb {
    /*   Compile B-allele frequencies
    */
    label 'penncnv'
    label 'cpu_16x'

    input:
    tuple file(reportFile), file(samplesheetFile), file(snpMap), file(genderFile), file(baflrrFiles)
    output:
    tuple file(reportFile), file(samplesheetFile), file(snpMap), file(genderFile), file(baflrrFiles), file('*.pfb'), emit: pfbFile

    script:
    """
    perl /PennCNV-1.0.5/compile_pfb.pl \
    --output output.pfb \
    --snpposfile ${snpMap} \
    ${baflrrFiles}
    """
}

process AdjustGCModel {
    /*   Make GC model file for correcting waviness
    */
    label 'penncnv'
    label 'cpu_16x'

    input:
    tuple file(reportFile), file(samplesheetFile), file(snpMap), file(genderFile), file(baflrrFiles), file(pfbFile)

    output:
    tuple file(reportFile), file(samplesheetFile), file(snpMap), file(genderFile), file(baflrrFiles), file(pfbFile), file('*.gcmodel'), emit: adjustedGcModel

    script:
    """
    perl /PennCNV-1.0.5/cal_gc_snp.pl \
    ${params.gcModel} \
    ${pfbFile} \
    --output output.gcmodel
    """
}

process DetectCNVs {
    /*   Detect copy numbers using hmmm model and pfb and adjustes gc files
    */
    label 'penncnv'
    label 'cpu_16x'

    input:
    tuple file(reportFile), file(samplesheetFile), file(snpMap), file(genderFile), file(baflrrFiles), file(pfbFile), file(adjustedGcModel)

    output:
    path 'all_cnvs.raw', emit: rawCNVs
    path 'detect_cnv.log', emit: qcLog

    script:
    """
    perl /PennCNV-1.0.5/detect_cnv.pl \
    --test \
    --confidence \
    --logfile detect_cnv.log \
    --hmmfile ${params.hmmFile} \
    --pfbfile ${pfbFile} \
    --gcmodel ${adjustedGcModel} \
    --output all_cnvs.raw \
    *.baflrr
    """
}

process CombineCNVs {
    /*   Combine CNV segments
    */
    label 'penncnv'
    label 'cpu_16x'

    input:
    tuple file(reportFile), file(samplesheetFile), file(snpMap), file(genderFile), file(baflrrFiles), file(pfbFile), file(adjustedGcModel)
    path rawCNVs

    output:
    path 'merged_all.rawcnv', emit: mergedCNVs

    script:
    """
    perl /PennCNV-1.0.5/clean_cnv.pl \
    combineseg \
    ${rawCNVs} \
    --signalfile ${snpMap} \
    --output merged_all.rawcnv \
    --fraction ${params.gapFraction}
    """
}

process FilterCNVs {
    /*   Filter CNV segments
    */
    label 'penncnv'
    label 'cpu_16x'
    
    input:
    path mergedCNVs
    path qcLog

    output:
    path 'merged_filtered_all.cnv', emit: filteredCNVs
    path 'CNV.PennCNV_qc_new.txt', emit: qcSum
    path 'penncnv_results', emit: resultDir

    script:
    """
    perl /PennCNV-1.0.5/filter_cnv.pl \
    ${mergedCNVs} \
    --numsnp ${params.numSNP} \
    --output merged_filtered_all.cnv \
    --qclogfile ${qcLog}\
    --qcpassout sampleall.qcpass \
    --qcsumout CNV.PennCNV_qc_new.txt \
    --qclrrsd ${params.qcLrrsd} \
    --qcbafdrift ${params.qcBafdrift} \
    --qcwf ${params.qcWf}

    mkdir penncnv_results
    cp merged_filtered_all.cnv penncnv_results/merged_filtered_all.cnv
    cp CNV.PennCNV_qc_new.txt penncnv_results/CNV.PennCNV_qc_new.txt
    """
}

process ConvertPennCNVs {
    /*  Convert CNVs to format used by EnsembleCNV
    */
    label 'penncnv'
    label 'cpu_one'
    publishDir "${params.resultsDir}/PennCNV", mode:'copy'

    input:
    tuple file(reportFile), file(samplesheetFile), file(snpMap), file(genderFile), file(baflrrFiles), file(pfbFile), file(adjustedGcModel)
    path resultDir

    output:
    path "${resultDir}", emit: resultDir

    script:
    """
    perl /PennCNV-1.0.5/convert_cnv.pl --intype penncnv --outtype tab ${resultDir}/merged_filtered_all.cnv > ${resultDir}/penncnv_all_cnvs.txt
    sed -i 's/.baflrr//g' ${resultDir}/penncnv_all_cnvs.txt
    sed -i 's/.baflrr//g' ${resultDir}/CNV.PennCNV_qc_new.txt
    cp ${pfbFile} ${resultDir}/SNP.pfb
    """
}

process ConvertPennCNVResults {
    /*  Convert cnvs to general benchmark format
    */
    label 'acnvbench'
    publishDir "${params.resultsDir}/CNVs/PennCNV", mode:'copy'

    input:
    tuple file(reportFile), file(samplesheetFile), file(snpMap), file(genderFile), file(baflrrFiles), file(pfbFile), file(adjustedGcModel)

    output:
    path "results.penncnv.txt", emit: callFile

    script:
    """
    Rscript ${PWD}/bin/convert_cnv_results.R -p merged_filtered_all.cnv
    """
}


workflow RunPennCNV {
    take:
        inputFiles
    main:
        CompilePfb(      inputFiles )
        AdjustGCModel(   CompilePfb.out.pfbFile )
        DetectCNVs(      AdjustGCModel.out.adjustedGcModel )
        CombineCNVs(     AdjustGCModel.out.adjustedGcModel,
                         DetectCNVs.out.rawCNVs )
        FilterCNVs(      CombineCNVs.out.mergedCNVs,
                         DetectCNVs.out.qcLog )
        ConvertPennCNVs( AdjustGCModel.out.adjustedGcModel,
                         FilterCNVs.out.resultDir )
        ConvertPennCNVResults( FilterCNVs.out.filteredCNVs )
    emit:
        callFile = ConvertPennCNVResults.out.callFile
        resultsDir = ConvertPennCNVs.out.resultDir
}
