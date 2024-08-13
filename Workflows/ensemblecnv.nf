process CreateInputFiles {
    /*   Creates directory for the first step, puts files resulting from PennCNV, QuantiSNP and iPattern in the right folders,
    combines all QuantiSNP results and creates lrr and baf matrices
    */
    label 'ensemblecnv'
    memory = { 64.GB * task.attempt}
    maxRetries = 4
    errorStrategy = 'retry'
    cache 'lenient'
    
    input:
    tuple file(reportFile), file(samplesheetFile), file(snpMap), file(genderFile), file(baflrrFiles)
    path penncnvResultsDir
    path quantisnpResultsDir
    path ipatternResultsDir

    output:
    path "data", emit: dataDir
    path "01_initial_call", emit: initialCallDir

    script:
    """
    cp -r /opt/ensembleCNV/01_initial_call ./01_initial_call
    mkdir 01_initial_call/run_QuantiSNP/results
    mkdir data

    cp ${reportFile} data/final_report.txt

    cp ${genderFile} data/Samples_Table.txt
    cp -a ${penncnvResultsDir}/. 01_initial_call/run_PennCNV/results/
    cp -a ${quantisnpResultsDir}/. 01_initial_call/run_QuantiSNP/results/
    cp -a ${ipatternResultsDir}/. 01_initial_call/run_iPattern/results/

    perl /opt/ensembleCNV/01_initial_call/finalreport_to_matrix_LRR_and_BAF/finalreport_matrix_LRR_BAF.pl \
    data/final_report.txt \
    01_initial_call/finalreport_to_matrix_LRR_and_BAF

    cp 01_initial_call/finalreport_to_matrix_LRR_and_BAF/SNP_pos.txt data/SNP_pos.txt
    """
}

process CreateInputMatrices {
    /* Transforms chromosome-wise lrr and baf matrices to rds files for each chromosome separately
    */
    label 'ensemblecnv'
    label 'mem_32'

    input:
    path initialCallDir
    each chr

    output:
    path "${initialCallDir.baseName}_chr${chr}", emit: chrDir

    script:
    """
    mkdir "${initialCallDir.baseName}_chr${chr}"
    cp -r ${initialCallDir} "${initialCallDir.baseName}_chr${chr}"

    Rscript /opt/ensembleCNV/01_initial_call/finalreport_to_matrix_LRR_and_BAF/transform_from_tab_to_rds.R \
    --input ${initialCallDir}/finalreport_to_matrix_LRR_and_BAF \
    --output "${initialCallDir.baseName}_chr${chr}" \
    --startChr ${chr} \
    --endChr ${chr}
    """
}

process CombineInputMatrices {
    /* Combines chromosome-wise lrr and baf rds files of all separate chromosomes into one BAF and LRR directory
    */
    input:
    path chrs

    output:
    path "RDS", emit: rds

    script:
    """
    mkdir -p RDS/
    mkdir -p RDS/BAF
    mkdir -p RDS/LRR
    for chr in ${chrs}; do
        cp \${chr}/BAF/*.rds RDS/BAF
        cp \${chr}/LRR/*.rds RDS/LRR
    done
    """
}

process PcaOnLRR {
    /*   Execute principal component analysis on log likelihood ratios to test
         for potential batch effects
    */
    label 'ensemblecnv'

    input:
    path dataDir

    output:
    path "02_batch_effect", emit: batchEffectDir

    script:
    """
    cp -r /opt/ensembleCNV/02_batch_effect ./02_batch_effect
    
    Rscript 02_batch_effect/PCA_on_LRR/step.1.down.sampling.R \
    ${dataDir}/SNP_pos.txt \
    02_batch_effect/PCA_on_LRR

    perl 02_batch_effect/PCA_on_LRR/step.2.LRR.matrix.pl \
    02_batch_effect/PCA_on_LRR/snps.down.sample.txt \
    ${dataDir}/final_report.txt \
    02_batch_effect/PCA_on_LRR/LRR_matrix_for_PCA.txt

    Rscript 02_batch_effect/PCA_on_LRR/step.3.LRR.pca.R \
    02_batch_effect/PCA_on_LRR/ \
    LRR_matrix_for_PCA.txt
    """
}

process PcaOnSummaryStats {
    /*   Execute principal component analysis on sample level statistics generated
         by PennCNV, QuantiSNP and iPattern to test for potential batch effects
    */
    label 'ensemblecnv'

    input:
    path initialCallDir

    output:
    path "02_batch_effect", emit: batchEffectDir

    script:
    """
    cp -r /opt/ensembleCNV/02_batch_effect ./02_batch_effect

    Rscript 02_batch_effect/PCA_on_summary_stats/step.1.prepare.stats.R \
    ${initialCallDir}/run_iPattern/results \
    ${initialCallDir}/run_PennCNV/results \
    ${initialCallDir}/run_QuantiSNP/results/res \
    02_batch_effect/PCA_on_summary_stats

    Rscript 02_batch_effect/PCA_on_summary_stats/step.2.stats.PCA.R \
    02_batch_effect/PCA_on_summary_stats
    """
}
process CreateCNVR {
    /*   Identification of CNVRs: regions in which CNVs called from different individuals
         from different callers substantially overlap
    */
    label 'ensemblecnv'

    input:
    path dataDir
    path initialCallDir

    output:
    path "03_create_CNVR", emit: cnvrDir

    script:
    """
    cp -r /opt/ensembleCNV/03_create_CNVR ./03_create_CNVR
    
    Rscript 03_create_CNVR/step.1.CNV.data.R \
    03_create_CNVR \
    ${initialCallDir}/run_iPattern/results/ipattern_all_calls.txt \
    ${initialCallDir}/run_PennCNV/results/penncnv_all_cnvs.txt \
    ${initialCallDir}/run_QuantiSNP/results/quantisnp.cnv \
    ${dataDir}/Samples_Table.txt

    Rscript /opt/ensembleCNV/03_create_CNVR/step.2.create.CNVR.R \
    --icnv 03_create_CNVR/cnv.ipattern.txt \
    --pcnv 03_create_CNVR/cnv.penncnv.txt \
    --qcnv 03_create_CNVR/cnv.quantisnp.txt \
    --snp ${initialCallDir}/finalreport_to_matrix_LRR_and_BAF/SNP_pos.txt \
    --centromere ${params.centromere} \
    --output 03_create_CNVR
    """
}

process GenotypeCNVs {
    /*   Re-genotyping of CN status per individual by fitting a local likelihood model
    */
    label 'ensemblecnv'

    input:
    path initialCallDir
    path cnvrDir
    path rds

    output:

    path "04_CNV_genotype", emit: toGenotypeDir
    path "04.genotype.step2.jobs.txt", emit: jobFile

    script:
    """
    cp -r /opt/ensembleCNV/04_CNV_genotype ./04_CNV_genotype
    mkdir -p 04_CNV_genotype/data/
    mkdir -p 04_CNV_genotype/results/
    cp -r ${rds} 04_CNV_genotype/RDS

    cp ${initialCallDir}/run_PennCNV/results/SNP.pfb 04_CNV_genotype/data/SNP.pfb
    cp ${cnvrDir}/cnvr_clean.txt 04_CNV_genotype/data/cnvr_clean.txt
    cp ${cnvrDir}/cnv_clean.txt 04_CNV_genotype/data/cnv_clean.txt
    cp ${initialCallDir}/run_PennCNV/results/CNV.PennCNV_qc_new.txt 04_CNV_genotype/data/samples_QC.txt
    sed -i 's/.baflrr//g' 04_CNV_genotype/data/samples_QC.txt

    Rscript 04_CNV_genotype/step.1.split.cnvrs.into.batches.R \
    -i 04_CNV_genotype/data/cnvr_clean.txt \
    -o 04_CNV_genotype/data/cnvr_batch.txt \
    -n ${params.maxCNVRs}

    Rscript 04_CNV_genotype/step.2.submit.jobs.nextflow.R \
    --type 0 \
    --script     04_CNV_genotype \
    --sourcefile 04_CNV_genotype/scripts \
    --datapath   04_CNV_genotype/data \
    --matrixpath ${rds} \
    --resultpath 04_CNV_genotype/results \
    --joblog     04_CNV_genotype/results
    """
}

process SubmitGenotypeJobs {
    /* Execute jobs generated by GenotypeCNVs process, necessary for parallelization
    */
    label 'ensemblecnv'
    label 'big_mem'

    input:
    each job
    path toGenotypeDir
    path rds

    output:
    path "04_CNV_genotype_*", emit: genotypedDir

    script:
    """
    tmp_dir="04_CNV_genotype_\$RANDOM"
    cp -r ${toGenotypeDir} \${tmp_dir}

    job=\$(echo "$job" | sed "s/04_CNV_genotype/\$tmp_dir/g")
    echo \$job
    \$job
    """
}

process ResubmitGenotypeCNVs {
    /*   Checking which jobs from SubmitGenotypeJobs process failed to finish and
    make resubmission jobs
    */
    label 'ensemblecnv'

    input:
    path genotypedDirs, stageAs: "04_CNV_genotype_*/*"
    path toGenotypeDir
    path rds

    output:
    path "04_CNV_regenotype", emit: toRegenotypeDir
    path "04.genotype.step3.jobs.txt", emit: jobFile

    script:
    """
    cp -r /opt/ensembleCNV/04_CNV_genotype ./04_CNV_regenotype
    cp -r ${toGenotypeDir}/data/. ./04_CNV_regenotype/data
    
    mkdir -p 04_CNV_regenotype/results
    for dir in ${genotypedDirs}; do
        cp -r \${dir}/. ./04_CNV_regenotype/
    done

    Rscript 04_CNV_regenotype/step.3.check.and.resubmit.jobs.nextflow.R \
    --flag 1 \
    --script     04_CNV_regenotype \
    --sourcefile 04_CNV_regenotype/scripts \
    --datapath   04_CNV_regenotype/data \
    --matrixpath ${rds} \
    --resultpath 04_CNV_regenotype/results \
    --joblog     04_CNV_regenotype/results
    """
}

process ResubmitGenotypeJobs {
    /* Execute jobs generated by ResubmitGenotypeCNVs process, necessary for parallelization
    */
    label 'ensemblecnv'
    label 'big_mem'

    input:
    each job
    path toRegenotypeDir
    path rds

    output:
    path "04_CNV_regenotype_*", emit: regenotypedDir

    script:
    """
    tmp_dir="04_CNV_regenotype_\$RANDOM"
    cp -r ${toRegenotypeDir} \${tmp_dir}
    
    if [ "${job}" != "" ];
    then
        job=\$(echo "$job" | sed "s/04_CNV_genotype/\$tmp_dir/g")
        \$job
    fi

    job=\$(echo "$job" | sed "s/04_CNV_genotype/\$tmp_dir/g")
    \$job
    """
}

process PredictGenotypeResults {
    /*   Last step of genotyping, prediction of the genotyping results
    */
    label 'ensemblecnv'

    input:
    path regenotypedDirs, stageAs: "04_CNV_regenotype_*/*"
    path toRegenotypeDir

    output:
    path "04_CNV_regenotype", emit: predictedGenotypeDir

    script:
    """
    for dir in ${regenotypedDirs}; do
        cp -R -u -n \${dir}/. ${toRegenotypeDir}/
    done

    Rscript 04_CNV_regenotype/step.4.prediction.results.R \
    --datapath ${toRegenotypeDir}/data \
    --resultpath ${toRegenotypeDir}/results
    """
}

process BoundaryRefinement {
    /*   Further refine CNVR boundaries using highly correlated LRR signals across individuals
    */
    label 'ensemblecnv'

    input:
    path dataDir
    path predictedRegenotypedDir
    path rdsDir

    output:
    path "05_boundary_refinement", emit: toRefineDir
    path "05.refinement.step2.jobs.txt", emit: jobFile

    script:
    """
    cp -r /opt/ensembleCNV/05_boundary_refinement ./05_boundary_refinement
    
    mkdir -p 05_boundary_refinement/data/
    mkdir -p 05_boundary_refinement/results/

    cp ${dataDir}/SNP_pos.txt 05_boundary_refinement/data/SNP_pos.txt
    cp ${predictedRegenotypedDir}/results/cnvr_genotype.txt 05_boundary_refinement/data/cnvr_genotype.txt
    cp ${predictedRegenotypedDir}/results/matrix_CN.rds 05_boundary_refinement/data/matrix_CN.rds
    cp ${predictedRegenotypedDir}/results/matrix_GQ.rds 05_boundary_refinement/data/matrix_GQ.rds

    Rscript 05_boundary_refinement/step.1.common.CNVR.to.refine.R \
    --datapath 05_boundary_refinement/data \
    --resultpath 05_boundary_refinement/results \
    --freq ${params.freqCutOff}

    Rscript 05_boundary_refinement/step.2.submit.jobs.nextflow.R \
    --refinescript 05_boundary_refinement/CNVR.boundary.refinement.R \
    --rcppfile   05_boundary_refinement/refine.cpp \
    --datapath  05_boundary_refinement/data \
    --matrixpath ${rdsDir} \
    --centromere ${params.centromere} \
    --resultpath 05_boundary_refinement/results
    """
}

process SubmitBoundaryRefinementJobs {
    /* Execute jobs generated by BoundaryRefinement process, necessary for parallelization
    */
    label 'ensemblecnv'

    input:
    each job
    path toRefineDir
    path rds

    output:
    path "05_boundary_refinement_*", emit: refinedDir

    script:
    """
    tmp_dir="05_boundary_refinement_\$RANDOM"
    cp -r ${toRefineDir} \${tmp_dir}
    cp /opt/ensembleCNV/05_boundary_refinement/CNVR.boundary.refinement.R \${tmp_dir}

    job=\$(echo "$job" | sed "s/05_boundary_refinement/\$tmp_dir/g")
    \$job
    echo \$job
    """
}

process CleanResults {
    /*   Further refine CNVR boundaries using highly correlated LRR signals across individuals
    */
    label 'ensemblecnv'

    input:
    path initialCallDir
    path cnvrDir
    path refinedDirs
    path rds

    output:
    path "05a_regenotype_after_refinement", emit: cleanedDir
    path "04.genotype.step2.jobs.txt", emit: jobFile

    script:
    """
    cp -r /opt/ensembleCNV/05_boundary_refinement ./05a_regenotype_after_refinement
    
    mkdir -p 05a_regenotype_after_refinement
    for dir in ${refinedDirs}; do
        cp -r \${dir}/. 05a_regenotype_after_refinement
    done

    Rscript 05a_regenotype_after_refinement/step.3.clean.results.R \
    --resultpath 05a_regenotype_after_refinement/results

    cp ${initialCallDir}/run_PennCNV/results/SNP.pfb 05a_regenotype_after_refinement/data/SNP.pfb
    cp ${cnvrDir}/cnv_clean.txt 05a_regenotype_after_refinement/data/cnv_clean.txt
    cp ${initialCallDir}/run_PennCNV/results/CNV.PennCNV_qc_new.txt 05a_regenotype_after_refinement/data/samples_QC.txt
    sed -i 's/.baflrr//g' 05a_regenotype_after_refinement/data/samples_QC.txt

    Rscript /opt/ensembleCNV/04_CNV_genotype/step.1.split.cnvrs.into.batches.R \
    -i 05a_regenotype_after_refinement/results/cnvr_regenotype_after_refine.txt \
    -o 05a_regenotype_after_refinement/data/cnvr_batch.txt \
    -n ${params.maxCNVRs}

    Rscript /opt/ensembleCNV/04_CNV_genotype/step.2.submit.jobs.nextflow.R \
    --type 0 \
    --script     /opt/ensembleCNV/04_CNV_genotype \
    --sourcefile /opt/ensembleCNV/04_CNV_genotype/scripts \
    --datapath   05a_regenotype_after_refinement/data \
    --matrixpath ${rds} \
    --resultpath 05a_regenotype_after_refinement/results \
    --joblog     05a_regenotype_after_refinement/results
    """
}

process SubmitGenotypeAfterRefinementJobs {
    /* Execute jobs generated by CleanResults process, necessary for parallelization
    */
    label 'ensemblecnv'
    label 'big_mem'

    input:
    each job
    path cleanDir
    path rds

    output:
    path "05a_regenotype_after_refinement_*", emit: regenotypedDir

    script:
    """
    tmp_dir="05a_regenotype_after_refinement_\$RANDOM"
    cp -r ${cleanDir} \${tmp_dir}

    if [ "${job}" != "" ];
    then
        job=\$(echo "$job" | sed "s/05a_regenotype_after_refinement_/\$tmp_dir/g")
        \$job
    fi
    """
}

process ResubmitGenotypeAfterRefinement {
    /*   Checking which jobs from SubmitGenotypeAfterRefinementJobs process
    failed to finish and make resubmission jobs
    */
    label 'ensemblecnv'

    input:
    path regenotypedDirs
    path rds

    output:
    path "05b_regenotype_after_refinement", emit: toRegenotypeDir
    path "04.genotype.step3.jobs.txt", emit: jobFile

    script:
    """
    cp -r /opt/ensembleCNV/05_boundary_refinement ./05b_regenotype_after_refinement
    
    mkdir -p 05b_regenotype_after_refinement
    for dir in ${regenotypedDirs}; do
        cp -r \${dir}/. 05b_regenotype_after_refinement/
    done

    Rscript /opt/ensembleCNV/04_CNV_genotype/step.3.check.and.resubmit.jobs.nextflow.R \
    --flag 1 \
    --script     /opt/ensembleCNV/04_CNV_genotype/04_CNV_genotype \
    --sourcefile /opt/ensembleCNV/04_CNV_genotype/04_CNV_genotype/scripts \
    --datapath   05b_regenotype_after_refinement/data \
    --matrixpath ${rds} \
    --resultpath 05b_regenotype_after_refinement/results \
    --joblog     05b_regenotype_after_refinement/results
    """
}

process ResubmitGenotypeAfterRefinementJobs {
    /* Execute jobs generated by ResubmitGenotypeAfterRefinement process,
    necessary for parallelization
    */
    label 'ensemblecnv'
    label 'big_mem'

    input:
    each job
    path toRegenotypeDir

    output:
    path "05b_regenotype_after_refinement_*", emit: regenotypedDir

    script:
    """
    tmp_dir="05b_regenotype_after_refinement_\$RANDOM"
    cp -r ${toRegenotypeDir} \${tmp_dir}

    if [ "${job}" != "" ];
    then
        job=\$(echo "$job" | sed "s/05b_regenotype_after_refinement_/\$tmp_dir/g")
        \$job
    fi
    """
}

process PredictRegenotypeResults {
    /*   Last step of genotyping after boundary refinement, prediction of the genotyping results
    */
    label 'ensemblecnv'

    input:
    path regenotypedDirs

    output:
    path "05c_predict_regenotype_after_refinement", emit: predictedRegenotyping

    script:
    """
    mkdir -p 05c_predict_regenotype_after_refinement
    for dir in ${regenotypedDirs}; do
        cp -R -u -n \${dir}/. 05c_predict_regenotype_after_refinement/
    done
    cp /opt/ensembleCNV/04_CNV_genotype/step.4.prediction.results.R 05c_predict_regenotype_after_refinement/step.4.prediction.results.R

    Rscript 05c_predict_regenotype_after_refinement/step.4.prediction.results.R \
    --datapath 05c_predict_regenotype_after_refinement/data \
    --resultpath 05c_predict_regenotype_after_refinement/results
    """
}

process UpdateGenotypes {
    /*   Further refine CNVR boundaries after boundary refinement, using highly correlated
    LRR signals across individuals
    */
    label 'ensemblecnv'

    input:
    path beforeRefineDir
    path afterRefinedDir

    output:
    path "results", emit: updatedGenotypes

    script:
    """
    mkdir -p results

    Rscript /opt/ensembleCNV/05_boundary_refinement/step.4.update.genotype.matrix.R \
    --matrixbeforerefine ${beforeRefineDir}/results \
    --matrixrefine ${afterRefinedDir}/results \
    --refinepath ${afterRefinedDir}/results \
    --output results
    """
}

process FilterOnGQ {
    /*   Filter results on GQ score threshold to get final results
    */
    label 'ensemblecnv'
    publishDir "${params.resultsDir}/EnsembleCNV", mode:'copy'


    input:
    path updatedGenotypes

    output:
    path "results", emit: filteredCNVs

    script:
    """
    cp -r /opt/ensembleCNV/06_performance_assessment ./06_performance_assessment
    cp -r ${updatedGenotypes}/. 06_performance_assessment/

    Rscript 06_performance_assessment/step.2.set.GQ.generate.results.R \
    --matrixCN 06_performance_assessment/matrix_CN_final.rds \
    --matrixGQ 06_performance_assessment/matrix_GQ_final.rds \
    --cnvrfile 06_performance_assessment/cnvr_final.txt \
    --resultpath results \
    --gqscore ${params.gqScore}
    """
}

process ConvertEnsembleCNVs  {
    /*  Convert cnvs to general benchmark format
    */
    label 'acnvbench'
    publishDir "${params.resultsDir}/CNVs/EnsembleCNV", mode:'copy'

    input:
    path filteredCNVs

    output:
    path 'results.ensemblecnv.txt', emit: callFile

    script:
    """
    Rscript ${PWD}/bin/convert_cnv_results.R \
    --ensemblecnv_matrixcn ${filteredCNVs}/cnvr_after_GQ.txt \
    --ensemblecnv_samplecnv ${filteredCNVs}/matrix_CN_after_GQ.rds
    """
}

workflow RunEnsembleCNV {
    take:
        inputFiles
        penncnvResultsDir
        quantisnpResultsDir
        ipatternResultsDir
    main:
        if ( params.chrNames == "") {
            chrNames =                          Channel.of ( 1..22)
        }
        else {
            chrNames =                          Channel .from( params.chrNames.tokenize(",") )
        }
        CreateInputFiles(                       inputFiles,
                                                penncnvResultsDir,
                                                quantisnpResultsDir,
                                                ipatternResultsDir )
        CreateInputMatrices(                    CreateInputFiles.out.initialCallDir,
                                                chrNames )
        CombineInputMatrices(                   CreateInputMatrices.out.chrDir.collect() )
        // PcaOnLRR(                               CreateInputFiles.out.dataDir )
        // PcaOnSummaryStats(                      CreateInputFiles.out.initialCallDir )
        CreateCNVR(                             CreateInputFiles.out.dataDir,
                                                CreateInputFiles.out.initialCallDir )
        GenotypeCNVs(                           CreateInputFiles.out.initialCallDir,
                                                CreateCNVR.out.cnvrDir,
                                                CombineInputMatrices.out.rds )
        genotypeJobs =                          GenotypeCNVs.out.jobFile.splitText()
        SubmitGenotypeJobs(                     genotypeJobs,
                                                GenotypeCNVs.out.toGenotypeDir,
                                                CombineInputMatrices.out.rds)
        ResubmitGenotypeCNVs(                   SubmitGenotypeJobs.out.genotypedDir.collect(),
                                                GenotypeCNVs.out.toGenotypeDir,
                                                CombineInputMatrices.out.rds )
        regenotypeJobs =                        ResubmitGenotypeCNVs.out.jobFile.splitText()
        ResubmitGenotypeJobs(                   regenotypeJobs,
                                                ResubmitGenotypeCNVs.out.toRegenotypeDir,
                                                CombineInputMatrices.out.rds )
        PredictGenotypeResults(                 ResubmitGenotypeJobs.out.regenotypedDir.collect(),
                                                ResubmitGenotypeCNVs.out.toRegenotypeDir )
        BoundaryRefinement(                     CreateInputFiles.out.dataDir,
                                                PredictGenotypeResults.out.predictedGenotypeDir,
                                                CombineInputMatrices.out.rds )
        refinementJobs =                        BoundaryRefinement.out.jobFile.splitText()
        SubmitBoundaryRefinementJobs(           refinementJobs,
                                                BoundaryRefinement.out.toRefineDir,
                                                CombineInputMatrices.out.rds )
        CleanResults(                           CreateInputFiles.out.initialCallDir,
                                                CreateCNVR.out.cnvrDir,
                                                SubmitBoundaryRefinementJobs.out.refinedDir.collect(),
                                                CombineInputMatrices.out.rds )
        regenotypeAfterRefinementJobs =         CleanResults.out.jobFile.splitText()
        SubmitGenotypeAfterRefinementJobs(      regenotypeAfterRefinementJobs,
                                                CleanResults.out.cleanedDir,
                                                CombineInputMatrices.out.rds )
        ResubmitGenotypeAfterRefinement(        SubmitGenotypeAfterRefinementJobs.out.regenotypedDir.collect(),
                                                CombineInputMatrices.out.rds )
        resubmitRegenotypeJobs =                ResubmitGenotypeAfterRefinement.out.jobFile.splitText()
        ResubmitGenotypeAfterRefinementJobs(    resubmitRegenotypeJobs,
                                                ResubmitGenotypeAfterRefinement.out.toRegenotypeDir )
        PredictRegenotypeResults(               ResubmitGenotypeAfterRefinementJobs.out.regenotypedDir.collect() )
        UpdateGenotypes(                        PredictGenotypeResults.out.predictedGenotypeDir,
                                                PredictRegenotypeResults.out.predictedRegenotyping )
        FilterOnGQ(                             UpdateGenotypes.out.updatedGenotypes )
        ConvertEnsembleCNVs(                    FilterOnGQ.out.filteredCNVs )
}
