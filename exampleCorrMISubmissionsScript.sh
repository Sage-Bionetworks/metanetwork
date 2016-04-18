#!/bin/bash

#number of cores to reserve for job
nthreads=2

#full s3 path where networks will go
s3="s3://metanetworks/testNetwork/"

#location of data file
dataFile="/shared/testNetwork/testData.csv"

#location of metanetwork synapse scripts
pathv="/shared/metanetworkSynapse/"

#output path for temporary result file prior to pushing to s3/synapse
outputpath="/shared/testNetwork/"

#path within s3
s3b="testNetwork"

#id of folder on Synapse that files will go to
parentId="syn5646308"

#path to csv file with annotations to add to file on Synapse
annotationFile="/shared/testNetwork/annoFile.txt"

#path to csv file with provenance to add to file on synapse
provenanceFile="/shared/testNetwork/provenanceFile.txt"

#path to error output
errorOutput="/shared/testNetwork/error.txt"

#path to out output
outOutput="/shared/testNetwork/out.txt"


echo "qsub -v s3=$s3,dataFile=$dataFile,pathv=$pathv,c3net=1,mrnet=1,aracne=1,wgcnaST=1,wgcnaTOM=1,sparrowZ=0,lassoCV1se=0,ridgeCV1se=0,genie3=0,tigress=0,numberCore=$nthreads,outputpath=$outputpath,s3b=$s3b,parentId=$parentId,annotationFile=$annotationFile,provenanceFile=$provenanceFile -pe orte $nthreads -S /bin/bash -V -cwd -N testingNetworkFunctionality -e $errorOutput -o $outOutput $pathv/buildNet.sh"

qsub -v s3=$s3,dataFile=$dataFile,pathv=$pathv,c3net=1,mrnet=1,aracne=1,wgcnaST=1,wgcnaTOM=1,sparrowZ=0,lassoCV1se=0,ridgeCV1se=0,genie3=0,tigress=0,numberCore=$nthreads,outputpath=$outputpath,s3b=$s3b,parentId=$parentId,annotationFile=$annotationFile,provenanceFile=$provenanceFile -pe orte $nthreads -S /bin/bash -V -cwd -N testingNetworkFunctionality -e $errorOutput -o $outOutput $pathv/buildNet.sh
