#!/bin/sh
root -l -q 'findFailedJobs.C("TTJets_DiLept_MC_MC2018")'
root -l -q 'findFailedJobs.C("TTJets_SingleLeptFromT_MC2018")'
root -l -q 'findFailedJobs.C("ZJetsToNuNu_HT_MC2018")'
root -l -q 'findFailedJobs.C("WJetsToLNu_HT_MC2018")'
root -l -q 'findFailedJobs.C("QCD_HT_MC2018")'
root -l -q 'findFailedJobs.C("ST__MC2018")'
#root -l -q 'findFailedJobs.C("Rare_MC2018")'
#root -l -q 'findFailedJobs.C("Rare_missingMC2018")'
root -l -q 'findFailedJobs.C("Rare_MC2018MC2017MC2016")'

root -l -q 'findFailedJobs.C("TChiWZ_1000_100_MC2018")'
root -l -q 'findFailedJobs.C("TChiWZ_800_100_MC2018")'
root -l -q 'findFailedJobs.C("TChiWZ_800_400_MC2018")'
root -l -q 'findFailedJobs.C("TChiWZ_800_700_MC2018")'
root -l -q 'findFailedJobs.C("TChiWZ_600_100_MC2018")'
root -l -q 'findFailedJobs.C("TChiWZ_1000_1_MC2018")'
root -l -q 'findFailedJobs.C("TChiWZ_800_1_MC2018")'
root -l -q 'findFailedJobs.C("TChiWZ_600_1_MC2018")'

root -l -q 'findFailedJobs.C("TChiWW_800_100_MC2018")'

# root -l -q 'findFailedJobs.C("TChiWH_800_100_MC2018")'
# root -l -q 'findFailedJobs.C("TChiWH_800_200_MC2018")'
# root -l -q 'findFailedJobs.C("TChiWH_800_400_MC2018")'
# root -l -q 'findFailedJobs.C("TChiWH_800_600_MC2018")'
# root -l -q 'findFailedJobs.C("TChiWH_600_100_MC2018")'
# root -l -q 'findFailedJobs.C("TChiWH_600_200_MC2018")'
# root -l -q 'findFailedJobs.C("TChiWH_600_400_MC2018")'
# root -l -q 'findFailedJobs.C("TChiWH_300_1_MC2018")'

hadd -f TTJets_MC2018.root TTJets_DiLept_MC_MC2018.root TTJets_SingleLeptFromT_MC2018.root
