#!/bin/sh

executable=$1
inputFileTag=$2
outputFileTag=$3
#commitHash=$4
datasetName=$4
currDir=$(pwd)
######################################
# SETUP CMSSW STUFF...
######################################
source /cvmfs/cms.cern.ch/cmsset_default.sh
#export SCRAM_ARCH=slc6_amd64_gcc700
#scram project CMSSW 10_2_11_patch1
#cd 10_2_11_patch1/src
scram project CMSSW_10_3_1
cd CMSSW_10_3_1/src/
eval `scramv1 runtime -sh`
pwd

######################################
# SETUP PRIVATE STUFF...
######################################
echo "ls"
pwd
cd -
#cd $currDir
echo "RUNNING ANALYSIS"
pwd
./$executable $inputFileTag $outputFileTag $datasetName
echo "processed. ls"
ls
echo "COPYING OUTPUT"

xrdcp -f $outputFileTag root://cmseos.fnal.gov//store/user/vhegde/RA2bTrees_V17/EW_SUSY/skims/ForFullyBoosted/
rm $outputFileTag
