process SetQuantiSNPParameters {
    /*   Replace existing parameters with parameters set by nextflow.config ana
    generate output directory
    */
    label 'quantisnp'
    label 'cpu_one'

    input:
    tuple file(reportFile), file(samplesheetFile), file(snpMap), file(genderFile), file(baflrrFiles)

    output:
    path 'alt_params.dat', emit: paramFile

    script:
    """
    sed 's/nComp\t3/nComp\t${params.nComp}/g' ${params.quantisnpParams} > alt_params.dat
    sed -i 's/v\t10/v\t${params.degreesFreedom}/g' alt_params.dat
    sed -i 's/nu_alpha\t1/nu_alpha\t${params.nuAlpha}/g' alt_params.dat
    sed -i 's/nu_beta\t1/nu_beta\t${params.nuBeta}/g' alt_params.dat
    sed -i 's/w_alpha\t1e4/w_alpha\t${params.wAlpha}/g' alt_params.dat
    sed -i 's/q_alpha\t10/q_alpha\t${params.qAlpha}/g' alt_params.dat
    sed -i 's/tau\t1e4/tau\t${params.tau}/g' alt_params.dat
    sed -i 's/S_alpha\t3/S_alpha\t${params.sAlpha}/g' alt_params.dat
    sed -i 's/S_alpha_homdel\t3/S_alpha_homdel\t${params.sAlphaHomdel}/g' alt_params.dat
    sed -i 's/longChromosome\t2e6/longChromosome\t${params.longChromosome}/g' alt_params.dat
    """
}

process GetSampleNames {
    /*   Get all individual sample names to allow for parallel processing of samples
    */
    label 'cpu_one'
    
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

process QuantiSNP {
    /*   Execute run quantisnp script
    */
    label 'quantisnp'
    label 'cpu_one'
    label 'mem_16'

    input:
    tuple file(reportFile), file(samplesheetFile), file(snpMap), file(genderFile), file(baflrrFiles)
    file paramFile
    each sample

    output:
    path 'quantisnp_output*', emit: outDir

    script:
    """
    sampleName=`(echo "${sample}"|tr -d '\n')`
    mkdir -p "quantisnp_output_\${sampleName}"
    file="\${sampleName}.baflrr"
    gender=(`grep -w \${sampleName} ${genderFile} | awk '{print \$2}'`);

    if [ \${gender} == "Unknown" ]
    then
        /opt/quantisnp/linux64/run_quantisnp2.sh \\
        /opt/MATLAB/MATLAB_Compiler_Runtime/v79/ \\
        --input-files \${file} \\
        --outdir "quantisnp_output_\${sampleName}" \\
        --sampleid \${sampleName} \\
        --levels ${params.levelsFile} \\
        --config ${paramFile}
    else
        /opt/quantisnp/linux64/run_quantisnp2.sh \\
        /opt/MATLAB/MATLAB_Compiler_Runtime/v79/ \\
        --input-files \${file} \\
        --outdir "quantisnp_output_\${sampleName}" \\
        --sampleid \${sampleName} \\
        --gender \${gender} \\
        --levels ${params.levelsFile} \\
        --config ${paramFile}
    fi
    """
}

process CreateResultDir {
    /* Create input format for EnsembleCNV input
    */
    label 'ensemblecnv'
    label 'cpu_one'
    publishDir "${params.resultsDir}/QuantiSNP", mode:'copy'
    
    input:
    path resultDirs

    output:
    path('quantisnp_output'), emit: resultDir

    script:
    """
    mkdir 'quantisnp_output'
    mkdir 'quantisnp_output/res'
    
    cp -r $resultDirs quantisnp_output/res
    cd quantisnp_output/res

    for dir in *; do
        mv "\$dir" "\$(echo "\$dir" | cut -c18-)";
    done

    pwd

    cd ../..

    pwd

    perl /opt/ensembleCNV/01_initial_call/run_QuantiSNP/step.3.combine.QuantiSNP.pl \
    --in_dir quantisnp_output/res \
    --out_dir quantisnp_output
    """
}

// process MergeQuantiSNPCNVs {
//     /* Merge CNVs and QC files from quantisnp to single call file
//     */
//     label 'acnvbench'
//     label 'cpu_one'
//     publishDir "${params.resultsDir}/QuantiSNP", mode:'copy'

//     input:
//     path resultDir

//     output:
//     path 'merged_quantisnp_output/all.merged.quantisnp.qc', emit: mergedQC
//     path 'merged_quantisnp_output/all.merged.quantisnp.cnv', emit: mergedCNVs

//     script:
//     """
//     mkdir merged_quantisnp_output
//     cd ${resultDir}/res

//     files=(`ls -d */`)
//     head -n 1 "\${files[0]}\${files[0]%/}.cnv" > all.merged.quantisnp.cnv
//     head -n 1 "\${files[0]}\${files[0]%/}.qc" > all.merged.quantisnp.qc

//     for file in \${files[*]}; do
//         tail -n +2 \${file}\${file%/}.cnv >> all.merged.quantisnp.cnv
//         tail -n +2 \${file}\${file%/}.qc >> all.merged.quantisnp.qc
//     done

//     cd ../..
//     mv 'quantisnp_output/res/all.merged.quantisnp.qc' 'merged_quantisnp_output/all.merged.quantisnp.qc'
//     mv 'quantisnp_output/res/all.merged.quantisnp.cnv' 'merged_quantisnp_output/all.merged.quantisnp.cnv'
//     """
// }

process ConvertQuantiSNPCNVs  {
    /*  Convert cnvs to general benchmark format
    */
    label 'acnvbench'
    label 'cpu_one'
    publishDir "${params.resultsDir}/CNVs/QuantiSNP", mode:'copy'

    input:
    each cnvDir

    output:
    tuple path('*.quantisnp.txt'), env(sample), val("QuantiSNP"), emit: individualCallFiles

    script:
    """
    sample=`(echo ${cnvDir.baseName} | sed "s/quantisnp_output_//g")`
    echo \$sample
    Rscript ${PWD}/bin/convert_cnv_results.R -q ${cnvDir}/*.cnv -o \${sample}
    """
}

process CombineQuantiSNPCNVs {
    /*  
    */
    label 'acnvbench'
    publishDir "${params.resultsDir}/CNVs/QuantiSNP", mode:'copy'

    input:
    path callFile

    output:
    path "quantisnp.results.txt", emit: callFile

    script:
    """
    callfiles=(*.quantisnp.txt)
    head -n 1 "\${callfiles[0]}" > quantisnp.results.txt
    for file in \${callfiles[*]}; do
        tail -n +2 \${file} >>  quantisnp.results.txt
    done
    """
}

workflow RunQuantiSNP {
    take:
        inputFiles
    main:
        SetQuantiSNPParameters( inputFiles )
        baflrrFiles = inputFiles.map { it[4] }.flatten()
        GetSampleNames( inputFiles )
        samples = GetSampleNames.out.sampleFile.splitText()
        QuantiSNP             ( inputFiles,
                                SetQuantiSNPParameters.out.paramFile,
                                samples )
        CreateResultDir       ( QuantiSNP.out.outDir.collect() )
        // MergeQuantiSNPCNVs    ( CreateResultDir.out.resultDir )
        ConvertQuantiSNPCNVs  ( QuantiSNP.out.outDir.collect() )
        callFilePaths = ConvertQuantiSNPCNVs.out.individualCallFiles.map{ callFile, sampleID, platform -> callFile }.collect()
        CombineQuantiSNPCNVs  ( callFilePaths )
    emit:
        callFile = CombineQuantiSNPCNVs.out.callFile
        individualCallFiles = ConvertQuantiSNPCNVs.out.individualCallFiles
        resultsDir = CreateResultDir.out.resultDir
}
