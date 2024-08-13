process SetiPatternParameters {
    /*   Replace existing parameters with parameters set by nextflow.config
    */
    label 'ipattern'
    label 'cpu_one'

    output:
    path 'ipattern.conf', emit: paramFile

    script:
    """
    sed 's/winSize = 9/winSize = ${params.winSize}/g' ${params.ipatternConf} > ipattern.conf
    sed -i 's/peakSeparation = 0.6/peakSeparation = ${params.peakSeparation}/g' ipattern.conf
    sed -i 's/bandWidth = 0.8/bandWidth = ${params.bandWidth}/g' ipattern.conf
    sed -i 's/maxProbeDistance = 400000/maxProbeDistance = ${params.maxProbeDistance}/g' ipattern.conf
    """
}

process CreateiPatternBatches {
    /*   Determine separation of samples into batches as iPattern is not able to work with large
    dataset sizes
    */
    label 'ensemblecnv'
    label 'cpu_one'
    label 'mem_32'

    input:
    tuple file(reportFile), file(samplesheetFile), file(snpMap), file(genderFile), file(paramFile)

    output:
    tuple file(reportFile), file(samplesheetFile), file(snpMap), file('*_genderfile_ipattern.txt'), file(paramFile), path("Reports"), emit: reportDir
    path 'Data_file_Batch_*.txt', emit: dataFiles

    script:
    """
    sed 's/emale//g' ${genderFile} > ${params.experimentName}_genderfile_ipattern.txt
    sed -i 's/ale//g' ${params.experimentName}_genderfile_ipattern.txt

    mkdir Reports
    /opt/ensembleCNV/01_initial_call/prepare_IPQ_input_file/finalreport_to_iPattern.pl -prefix Reports/${params.experimentName} -suffix _report.txt ${reportFile}
    ls Reports/*_report.txt > Data_file.txt
    split -l ${params.batchSize} -d --additional-suffix=".txt" Data_file.txt Data_file_Batch_
    rm Data_file.txt
    """
}

process SplitAllDataFiles {
    /*   Split report file into single sample files in iPattern format
    */
    label 'ipattern'
    label 'cpu_16x'
    
    input:
    tuple file(reportFile), file(samplesheetFile), file(snpMap), file(genderFile), file(paramFile), path(reportDir)
    each file(dataFile)

    output:
    tuple file(dataFile), file(genderFile), file(paramFile), path(reportDir), path('*_input'), path('*_step1'), emit: splitData

    script:
    """
    prefix="Data_file_"
    filename="${dataFile.baseName}"
    mkdir -p "\${filename#"\$prefix"}_input"
    mkdir -p "\${filename#"\$prefix"}_step1"

    split_all_data_files.py \
    -f ${dataFile} \
    -m "\${filename#"\$prefix"}_input/matchfile.txt" \
    --interdir "\${filename#"\$prefix"}_input" \
    --output-directory "\${filename#"\$prefix"}_step1" \
    --samplefilesuffix zzz \
    --samplefile "\${filename#"\$prefix"}_input/sample_list.txt" \
    --allsamplefile "\${filename#"\$prefix"}_input/all_samples.txt" \
    --noqsub
    """
}

process GenerateSortProbes {
    /*   Generate and sort probes
    */
    label 'ipattern'

    input:
    tuple file(dataFile), file(genderFile), file(paramFile), path(reportDir), path(inputDir), path(step1Dir)

    output:
    tuple file(dataFile), file(genderFile), file(paramFile), path(reportDir), path(inputDir), path(step1Dir), path('*_step2'), emit: sortedProbes

    script:
    """
    prefix="Data_file_"
    filename="${dataFile.baseName}"
    mkdir -p "\${filename#"\$prefix"}_step2"

    reports=(${reportDir}/*)

    split_sample_probe.py \
    -i \${reports[0]} \
    -b "\${filename#"\$prefix"}_step2/probes.txt"

    sort_probes.py \
    -p "\${filename#"\$prefix"}_step2/probes.txt" \
    -x "\${filename#"\$prefix"}_step2/index.txt" \
    -o "\${filename#"\$prefix"}_step2/sorted_probes.txt"
    """
}

process SortSampleIntensities {
    /*   Sort the intensity values
    */
    label 'ipattern'

    input:
    tuple file(dataFile), file(genderFile), file(paramFile), path(reportDir), path(inputDir), path(step1Dir), path(step2Dir)

    output:
    tuple file(dataFile), file(genderFile), file(paramFile), path(reportDir), path(inputDir), path(step1Dir), path(step2Dir), path('*_step3'), emit: sortedIntensities

    script:
    """
    prefix="Data_file_"
    filename="${dataFile.baseName}"
    mkdir -p "\${filename#"\$prefix"}_step3"

    cp -r ${step1Dir}/*.zzz \${filename#"\$prefix"}_step3/
    
    ilmn_sort_all_int.py \
    -d "\${filename#"\$prefix"}_step3" \
    -m ${inputDir}/all_samples.txt \
    -a zzz \
    -b sort \
    -x ${step2Dir}/sorted_probes.txt \
    --noqsub

    rm -r \${filename#"\$prefix"}_step3/*.zzz
    """
}

process BalanceSamples {
    /*   Channel balancing of the samples
    */
    label 'ipattern'

    input:
    tuple file(dataFile), file(genderFile), file(paramFile), path(reportDir), path(inputDir), path(step1Dir), path(step2Dir), path(step3Dir)

    output:
    tuple file(dataFile), file(genderFile), file(paramFile), path(reportDir), path(inputDir), path(step1Dir), path(step2Dir), path('*_step4'), emit: balancedSamples

    script:
    """
    prefix="Data_file_"
    filename="${dataFile.baseName}"
    mkdir -p "\${filename#"\$prefix"}_step4"

    cp -r ${step3Dir}/*.sort \${filename#"\$prefix"}_step4/

    illuminaBalance.py \
    -s "\${filename#"\$prefix"}_step4/" \
    -m ${inputDir}/all_samples.txt \
    -a sort \
    -b bln \
    --noqsub

    rm -r \${filename#"\$prefix"}_step4/*.sort
    """
}

process QuantileNormalization {
    /*   Applies quantile normalization to the intensity values
    */
    label 'ipattern'

    input:
    tuple file(dataFile), file(genderFile), file(paramFile), path(reportDir), path(inputDir), path(step1Dir), path(step2Dir), path(step4Dir)

    output:
    tuple file(dataFile), file(genderFile), file(paramFile), path(reportDir), path(inputDir), path(step1Dir), path(step2Dir), path('*_step5'), emit: quantileNormalized

    script:
    """
    prefix="Data_file_"
    filename="${dataFile.baseName}"
    mkdir -p "\${filename#"\$prefix"}_step5"

    cp -r ${step4Dir}/*.bln \${filename#"\$prefix"}_step5/

    python /opt/ipn_0.582/preprocess/ilmn/quantile.py \
    -d "\${filename#"\$prefix"}_step5" \
    -m ${inputDir}/all_samples.txt \
    -a bln \
    -b nml

    rm -r \${filename#"\$prefix"}_step5/*.bln
    """
}

process RescaleIntensities {
    /*   Rescaling of intensity values
    */
    label 'ipattern'

    input:
    tuple file(dataFile), file(genderFile), file(paramFile), path(reportDir), path(inputDir), path(step1Dir), path(step2Dir), path(step5Dir)

    output:
    tuple file(dataFile), file(genderFile), file(paramFile), path(reportDir), path(inputDir), path(step1Dir), path(step2Dir), path('*_step6'), emit: rescaledIntensities

    script:
    """
    prefix="Data_file_"
    filename="${dataFile.baseName}"
    mkdir -p "\${filename#"\$prefix"}_step6"

    cp -r ${step5Dir}/*.nml \${filename#"\$prefix"}_step6/

    python /opt/ipn_0.582/preprocess/ilmn/rescale.py \
    -d "\${filename#"\$prefix"}_step6" \
    -m ${inputDir}/all_samples.txt \
    -a nml \
    -b rescale

    rm -r \${filename#"\$prefix"}_step6/*.nml
    """
}

process VarianceNormalization {
    /*   Normalization of variance
    */
    label 'ipattern'

    input:
    tuple file(dataFile), file(genderFile), file(paramFile), path(reportDir), path(inputDir), path(step1Dir), path(step2Dir), path(step6Dir)

    output:
    tuple file(dataFile), file(genderFile), file(paramFile), path(reportDir), path(inputDir), path(step1Dir), path(step2Dir), path('*_step7'), emit: normalizedVariance
    path('sample.stats.txt'), emit: sampleStats

    script:
    """
    prefix="Data_file_"
    filename="${dataFile.baseName}"
    mkdir -p "\${filename#"\$prefix"}_step7"

    cp -r ${step6Dir}/*.rescale \${filename#"\$prefix"}_step7/

    python /opt/ipn_0.582/preprocess/ilmn/variance.py \
    -d "\${filename#"\$prefix"}_step7" \
    -m ${inputDir}/all_samples.txt \
    -a rescale \
    -b vn \
    -v sample.stats.txt
    
    rm -r \${filename#"\$prefix"}_step7/*.rescale
    """
}

process SplitByChromosome {
    /*   Split probes according to chromosome
    */
    label 'acnvbench'
    label 'mem_64'

    input:
    tuple file(dataFile), file(genderFile), file(paramFile), path(reportDir), path(inputDir), path(step1Dir), path(step2Dir), path(step7Dir)

    output:
    tuple file(dataFile), file(genderFile), file(paramFile), path(reportDir), path(inputDir), path(step1Dir), path(step2Dir), path('*_step8'), emit: splitChromosome

    script:
    """
    prefix="Data_file_"
    filename="${dataFile.baseName}"
    mkdir -p "\${filename#"\$prefix"}_step8"

    cp ${step2Dir}/sorted_probes.txt "\${filename#"\$prefix"}_step8/sorted_probes.txt"

    Rscript ${PWD}/bin/split_chr.R \
    -m ${inputDir}/all_samples.txt \
    -p "\${filename#"\$prefix"}_step8/sorted_probes.txt" \
    -s ${step7Dir} \
    -d "\${filename#"\$prefix"}_step8" \
    -x vn \
    -g ${genderFile}.txt \
    -q ${params.pqFile} \
    -z ${params.experimentName}
    """
}

process SplitByGender {
    /*   Split probes according to gender
    */
    label 'ipattern'
    label 'cpu_2x'

    input:
    tuple file(dataFile), file(genderFile), file(paramFile), path(reportDir), path(inputDir), path(step1Dir), path(step2Dir), path(step8Dir)

    output:
    tuple file(dataFile), file(genderFile), file(paramFile), path(reportDir), path(inputDir), path(step1Dir), path(step2Dir), path('*_step9'), emit: splitGender
    path '*_step9/*.int', emit: intFiles

    script:
    """
    prefix="Data_file_"
    filename="${dataFile.baseName}"
    mkdir -p "\${filename#"\$prefix"}_step9"

    cp -r ${step8Dir}/*.int \${filename#"\$prefix"}_step9/

    python /opt/ipn_0.582/preprocess/ilmn/split_by_gender.py \
    -s "\${filename#"\$prefix"}_step9" \
    -d "\${filename#"\$prefix"}_step9" \
	-g ${genderFile}
    """
}

process ExecuteIpattern {
    /*  
    */
    label 'ipattern'
    label 'cpu_16x'

    input:
    tuple file(dataFile), file(genderFile), file(paramFile), path(reportDir), path(inputDir), path(step1Dir), path(step2Dir), path(step9Dir), path(intFile)

    output:
    path "${intFile.baseName}_results", emit: resultsDir
    path "done.*", emit: doneFile

    script:
    """
    prefix="Data_file_"
    filename="${dataFile.baseName}"
    mkdir -p "${intFile.baseName}_results"

    echo ${intFile}

    conffile=`readlink -e ${paramFile}`
    echo \$conffile

    ${PWD}/bin/Write_ipattern_jobs.sh ${intFile} ${step9Dir} "${intFile.baseName}_results" \${conffile} ${dataFile.baseName} "${params.knownCnvrFile}"

    bash ${intFile}.qsub

    DONEFILE=`mktemp done.XXXXXXXXXXXXXXXXXXXX`
    touch \${DONEFILE}
    """
}

process MergeIpatternJobs {
    /*   
    */
    label 'ipattern'
    label 'cpu_1x'

    input:
    tuple file(dataFile), file(genderFile), file(paramFile), path(reportDir), path(inputDir), path(step1Dir), path(step2Dir), path(step9Dir), path(intFile)
    path(resultDirs, stageAs: "?/*")

    output:
    path "all_reports/*_all_calls.txt", emit: callFile
    path "all_bmps", emit: bmpFiles

    script:
    """
    mkdir -p all_reports
    mkdir -p all_bmps

    for dir in ${resultDirs}; do
        [ -f \${dir}/*.report ] && cp \${dir}/*.report all_reports
        [ -f \${dir}/*.bmp ] && cp \${dir}/*.bmp all_bmps
    done

    python /opt/ipn_0.582/ipn/merge_call_files_chg_names.py \
        -d all_reports \
        -o all_reports/${dataFile.baseName}_all_calls.txt
    """
}

process CombineIpatternCNVs {
    /*  Convert cnvs to general benchmark format
    */
    input:
    path '*_all_calls.txt'
    path(sampleStats, stageAs: "sample.stats.*.txt")
    publishDir "${params.resultsDir}/iPattern", mode:'copy'

    output:
    path 'ipattern_results/ipattern_all_calls.txt', emit: callFile
    path 'ipattern_results', emit: resultsDir

    script:
    """
    mkdir -p ipattern_results
    for file in ${sampleStats}; do
        cat \${file} >> "ipattern_results/ipattern_sample.stats.txt"
    done

    callfiles=(*_all_calls.txt)
    head -n 17 \${callfiles[0]} > "ipattern_results/ipattern_all_calls.txt"
    for file in *_all_calls.txt; do
        tail -n +18 \$file >> "ipattern_results/ipattern_all_calls.txt"
    done
    """
}

process ConvertIpatternCNVs {
    /*  Convert cnvs to general benchmark format
    */
    label 'acnvbench'
    publishDir "${params.resultsDir}/CNVs/iPattern", mode:'copy'

    input:
    path allCalls

    output:
    path 'results.ipattern.txt', emit: callFile

    script:
    """
    Rscript ${PWD}/bin/convert_cnv_results.R -i ${allCalls}
    """
}

workflow RuniPattern {
    take:
        inputFiles

    main:
        SetiPatternParameters()
        CreateiPatternBatches(     inputFiles.combine(SetiPatternParameters.out.paramFile ) )
        // CreateiPatternDirectories( CreateiPatternBatches.out.reportDir,
        //                            CreateiPatternBatches.out.dataFiles.flatten() )
        SplitAllDataFiles(         CreateiPatternBatches.out.reportDir,
                                    CreateiPatternBatches.out.dataFiles.flatten() )
        GenerateSortProbes(        SplitAllDataFiles.out.splitData )
        SortSampleIntensities(     GenerateSortProbes.out.sortedProbes )
        BalanceSamples(            SortSampleIntensities.out.sortedIntensities )
        QuantileNormalization(     BalanceSamples.out.balancedSamples )
        RescaleIntensities(        QuantileNormalization.out.quantileNormalized )
        VarianceNormalization(     RescaleIntensities.out.rescaledIntensities )
        SplitByChromosome(         VarianceNormalization.out.normalizedVariance )
        SplitByGender(             SplitByChromosome.out.splitChromosome )
        intFiles = SplitByGender.out.intFiles.flatten().unique()
        ipatternJobs = SplitByGender.out.splitGender.combine(intFiles)
        ExecuteIpattern(           ipatternJobs )
        MergeIpatternJobs(         ipatternJobs.first(),
                                   ExecuteIpattern.out.resultsDir.collect() )
        CombineIpatternCNVs(       MergeIpatternJobs.out.callFile.collect(),
                                   VarianceNormalization.out.sampleStats.collect() )
        ConvertIpatternCNVs(       CombineIpatternCNVs.out.callFile )
    emit:
        callFile = ConvertIpatternCNVs.out.callFile
        resultsDir =  CombineIpatternCNVs.out.resultsDir
}
