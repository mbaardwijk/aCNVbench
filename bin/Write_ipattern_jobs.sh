#!/bin/bash

infile=$1
srcdir=$2
tmpdir=$3
conffile=$4
batchname=$5
knownCnvrFile=$6

echo '#!/bin/bash' > ${infile}.qsub
echo "infile='${infile}'" >> ${infile}.qsub
echo "srcdir='${srcdir}'" >> ${infile}.qsub
echo "tmpdir='${tmpdir}'" >> ${infile}.qsub
echo "conffile='${conffile}'" >> ${infile}.qsub
#echo 'conffile=`readlink -e ${conffile}`' >> ${infile}.qsub

echo 'srcdir=`pwd`/${srcdir}' >> ${infile}.qsub
echo 'tmpdir=`pwd`/${tmpdir}' >> ${infile}.qsub

echo 'job_id=${JOB_ID:-$random}' >> ${infile}.qsub
echo 'thedir=$tmpdir/$job_id' >> ${infile}.qsub

echo 'mkdir -p $thedir' >> ${infile}.qsub
echo 'ln -sf ${conffile} ${thedir}/ipattern.conf' >> ${infile}.qsub
echo 'cd $thedir' >> ${infile}.qsub
echo 'ln -sf "${srcdir}/${infile}" "${thedir}/${infile}"' >> ${infile}.qsub
echo 'ln -sf "/opt/ipn_0.582/ipn/iPattern.Runner.R" ${thedir}/iPattern.Runner.R' >> ${infile}.qsub
echo 'ln -sf "/opt/ipn_0.582/ipn/iPattern.Func.R" ${thedir}/iPattern.Func.R' >> ${infile}.qsub
echo 'ln -sf "${knownCnvrFile}" ${thedir}/known.cnvr.txt' >> ${infile}.qsub

echo '/usr/bin/R --no-restore --no-save --no-readline ipattern.conf ${infile} < iPattern.Runner.R > /dev/null 2> ${tmpdir}/${infile}.$job_id.Rerr' >> ${infile}.qsub
#echo 'rm -f ipattern.conf ${infile} *.R known.cnvr.txt' >> ${infile}.qsub
#echo 'rm -f ${infile}' >> ${infile}.qsub
#echo 'rm -f *.R' >> ${infile}.qsub
#echo 'rm -f known.cnvr.txt' >> ${infile}.qsub

echo 'cp ${infile}.ipttn.txt ${tmpdir}' >> ${infile}.qsub
echo 'cd ..' >> ${infile}.qsub
echo 'tar -zcf ${infile}.ipttn.tar.gz $job_id' >> ${infile}.qsub

#echo 'rm -rf ${thedir}' >> ${infile}.qsub
echo 'cd ${tmpdir}' >> ${infile}.qsub

echo '/opt/ipn_0.582/ipn/iPatternRawReportParserAndFilter.pl < ./${infile}.ipttn.txt > ./${infile}.${batchname}.ipttn.txt.report' >> ${infile}.qsub
