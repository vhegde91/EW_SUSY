#!/bin/sh
root -l -q 'findFailedJobs.C("TTJets_DiLept_MC_MC2018")'
root -l -q 'findFailedJobs.C("TTJets_SingleLeptFromT_MC2018")'
root -l -q 'findFailedJobs.C("ZJetsToNuNu_HT_MC2018")'
root -l -q 'findFailedJobs.C("WJetsToLNu_HT_MC2018")'
#root -l -q 'findFailedJobs.C("QCD_HT_MC2018")'
root -l -q 'findFailedJobs.C("ST__MC2018")'
#root -l -q 'findFailedJobs.C("Rare_MC2018")'
#root -l -q 'findFailedJobs.C("Rare_missingMC2018")'
root -l -q 'findFailedJobs.C("Rare_MC2018MC2017MC2016")'

root -l -q 'findFailedJobs.C("TChiWZ_1000_1_MC2018")'
root -l -q 'findFailedJobs.C("TChiWZ_1000_100_MC2018")'
root -l -q 'findFailedJobs.C("TChiWZ_800_1_MC2018")'
root -l -q 'findFailedJobs.C("TChiWZ_800_100_MC2018")'
root -l -q 'findFailedJobs.C("TChiWZ_800_400_MC2018")'
root -l -q 'findFailedJobs.C("TChiWZ_600_1_MC2018")'
root -l -q 'findFailedJobs.C("TChiWZ_600_100_MC2018")'
root -l -q 'findFailedJobs.C("TChiWZ_300_1_MC2018")'

root -l -q 'findFailedJobs.C("TChiWW_800_100_MC2018")'

root -l -q 'findFailedJobs.C("TChiWH_1000_100_MC2018")'
root -l -q 'findFailedJobs.C("TChiWH_800_100_MC2018")'
root -l -q 'findFailedJobs.C("TChiWH_800_200_MC2018")'
root -l -q 'findFailedJobs.C("TChiWH_800_400_MC2018")'
root -l -q 'findFailedJobs.C("TChiWH_800_600_MC2018")'
root -l -q 'findFailedJobs.C("TChiWH_600_100_MC2018")'
root -l -q 'findFailedJobs.C("TChiWH_600_200_MC2018")'
root -l -q 'findFailedJobs.C("TChiWH_600_400_MC2018")'
root -l -q 'findFailedJobs.C("TChiWH_300_1_MC2018")'

hadd -f TTJets_MC2018.root TTJets_DiLept_MC_MC2018.root TTJets_SingleLeptFromT_MC2018.root
hadd -f VJets_MC2018.root WJetsToLNu_HT_MC2018.root ZJetsToNuNu_HT_MC2018.root
hadd -f TopVJets_MC2018.root VJets_MC2018.root TTJets_MC2018.root

# single lepton CR
root -l -q 'findFailedJobs.C("Ele_TTJets_DiLept_MC2018")'
root -l -q 'findFailedJobs.C("Ele_TTJets_SingleLeptFromT_MC2018")'
root -l -q 'findFailedJobs.C("Ele_WJetsToLNu_HT_MC2018")'
root -l -q 'findFailedJobs.C("Ele_QCD_HT_MC2018")'
root -l -q 'findFailedJobs.C("Ele_ST__MC2018")'
root -l -q 'findFailedJobs.C("Ele_Rare_MCMC2018")'
root -l -q 'findFailedJobs.C("Ele_data2018")'

root -l -q 'findFailedJobs.C("Muon_TTJets_DiLept_MC2018")'
root -l -q 'findFailedJobs.C("Muon_TTJets_SingleLeptFromT_MC2018")'
root -l -q 'findFailedJobs.C("Muon_WJetsToLNu_HT_MC2018")'
root -l -q 'findFailedJobs.C("Muon_QCD_HT_MC2018")'
root -l -q 'findFailedJobs.C("Muon_ST__MC2018")'
root -l -q 'findFailedJobs.C("Muon_Rare_MCMC2018")'
root -l -q 'findFailedJobs.C("Muon_data2018")'

hadd -f Ele_TTJets_MC2018.root Ele_TTJets_SingleLeptFromT_MC2018.root Ele_TTJets_DiLept_MC2018.root
hadd -f Muon_TTJets_MC2018.root Muon_TTJets_SingleLeptFromT_MC2018.root Muon_TTJets_DiLept_MC2018.root
