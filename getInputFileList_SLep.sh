#!/bin/sh
#for sample in DYJetsToLL_M-50_HT GJets_DR-0p4_HT QCD_HT TTGJets_ TTJets_DiLept TTJets_HT TTJets_SingleLeptFromT WGJets_MonoPhoton_PtG WJetsToLNu_HT ZJetsToNuNu_HT
for sample in TTJets_DiLept_MC TTJets_SingleLeptFromT ZJetsToNuNu_HT WJetsToLNu_HT QCD_HT ST_
do
#    grep $sample skimFilesSLep.txt | grep tree_SLe | grep MC2016 > Ele_${sample}_MC2016.txt
#    echo "root -l -q 'splitRunList.C(\"Ele_${sample}_MC2016.txt\",1)'"

#    grep $sample skimFilesSLep.txt | grep tree_SLe | grep MC2017 > Ele_${sample}_MC2017.txt
#    echo "root -l -q 'splitRunList.C(\"Ele_${sample}_MC2017.txt\",1)'"

    grep $sample skimFilesSLep.txt | grep tree_SLe | grep MC2018 > Ele_${sample}_MC2018.txt
    echo "root -l -q 'splitRunList.C(\"Ele_${sample}_MC2018.txt\",1)'"
#    grep $sample skimFilesSLep.txt | grep tree_SLm | grep MC2016 > Muon_${sample}_MC2016.txt
#    echo "root -l -q 'splitRunList.C(\"Muon_${sample}_MC2016.txt\",1)'"

#    grep $sample skimFilesSLep.txt | grep tree_SLm | grep MC2017 > Muon_${sample}_MC2017.txt
#    echo "root -l -q 'splitRunList.C(\"Muon_${sample}_MC2017.txt\",1)'"

    grep $sample skimFilesSLep.txt | grep tree_SLm | grep MC2018 > Muon_${sample}_MC2018.txt
    echo "root -l -q 'splitRunList.C(\"Muon_${sample}_MC2018.txt\",1)'"

done

for yEAR in MC2018
do
    grep TTTT skimFilesSLep.txt | grep tree_SLe | grep $yEAR > Ele_Rare_MC$yEAR.txt
    grep TTWJets skimFilesSLep.txt | grep tree_SLe | grep $yEAR >> Ele_Rare_MC$yEAR.txt
    grep TTZTo skimFilesSLep.txt | grep tree_SLe | grep $yEAR >> Ele_Rare_MC$yEAR.txt
    grep WWZ skimFilesSLep.txt | grep tree_SLe | grep $yEAR >> Ele_Rare_MC$yEAR.txt
    grep WWTo skimFilesSLep.txt | grep tree_SLe | grep $yEAR >> Ele_Rare_MC$yEAR.txt
    grep WZTo skimFilesSLep.txt | grep tree_SLe | grep $yEAR >> Ele_Rare_MC$yEAR.txt
    grep WZZ skimFilesSLep.txt | grep tree_SLe | grep $yEAR >> Ele_Rare_MC$yEAR.txt
    grep ZZTo skimFilesSLep.txt | grep tree_SLe | grep $yEAR >> Ele_Rare_MC$yEAR.txt
    echo "root -l -q 'splitRunList.C(\"Ele_Rare_MC${yEAR}.txt\",1)'"

    grep TTTT skimFilesSLep.txt | grep tree_SLm | grep $yEAR > Muon_Rare_MC$yEAR.txt
    grep TTWJets skimFilesSLep.txt | grep tree_SLm | grep $yEAR >> Muon_Rare_MC$yEAR.txt
    grep TTZTo skimFilesSLep.txt | grep tree_SLm | grep $yEAR >> Muon_Rare_MC$yEAR.txt
    grep WWZ skimFilesSLep.txt | grep tree_SLm | grep $yEAR >> Muon_Rare_MC$yEAR.txt
    grep WWTo skimFilesSLep.txt | grep tree_SLm | grep $yEAR >> Muon_Rare_MC$yEAR.txt
    grep WZTo skimFilesSLep.txt | grep tree_SLm | grep $yEAR >> Muon_Rare_MC$yEAR.txt
    grep WZZ skimFilesSLep.txt | grep tree_SLm | grep $yEAR >> Muon_Rare_MC$yEAR.txt
    grep ZZTo skimFilesSLep.txt | grep tree_SLm | grep $yEAR >> Muon_Rare_MC$yEAR.txt
    echo "root -l -q 'splitRunList.C(\"Muon_Rare_MC${yEAR}.txt\",1)'"
done

for sample in 2018
do
    grep tree_MET_$sample skimFilesSLep.txt | grep tree_SLe > Ele_data${sample}.txt
    echo "root -l -q 'splitRunList.C(\"Ele_data${sample}.txt\",1)'"
    grep tree_MET_$sample skimFilesSLep.txt | grep tree_SLm > Muon_data${sample}.txt
    echo "root -l -q 'splitRunList.C(\"Muon_data${sample}.txt\",1)'"
done
