#!/bin/sh
#for sample in DYJetsToLL_M-50_HT GJets_DR-0p4_HT QCD_HT TTGJets_ TTJets_DiLept TTJets_HT TTJets_SingleLeptFromT WGJets_MonoPhoton_PtG WJetsToLNu_HT ZJetsToNuNu_HT
for sample in TTJets_DiLept_MC TTJets_SingleLeptFromT ZJetsToNuNu_HT WJetsToLNu_HT QCD_HT ST_
do
#    grep $sample skimFilesSR.txt | grep 2016 > ${sample}_MC2016.txt
#    echo "root -l -q 'splitRunList.C(\"${sample}_MC2016.txt\",1)'"

#    grep $sample skimFilesSR.txt | grep 2017 > ${sample}_MC2017.txt
#    echo "root -l -q 'splitRunList.C(\"${sample}_MC2017.txt\",1)'"

    grep $sample skimFilesSR.txt | grep 2018 > ${sample}_MC2018.txt
    echo "root -l -q 'splitRunList.C(\"${sample}_MC2018.txt\",1)'"

done

for yEAR in 2018
do
    grep TTTT skimFilesSR.txt | grep $yEAR > Rare_MC$yEAR.txt
    grep TTWJets skimFilesSR.txt |grep $yEAR >> Rare_MC$yEAR.txt
    grep TTZTo skimFilesSR.txt |grep $yEAR >> Rare_MC$yEAR.txt
    grep WWZ skimFilesSR.txt |grep $yEAR >> Rare_MC$yEAR.txt
    grep WWTo skimFilesSR.txt |grep $yEAR >> Rare_MC$yEAR.txt
    grep WZTo skimFilesSR.txt |grep $yEAR >> Rare_MC$yEAR.txt
    grep WZZ skimFilesSR.txt |grep $yEAR >> Rare_MC$yEAR.txt
    grep ZZTo skimFilesSR.txt |grep $yEAR >> Rare_MC$yEAR.txt
done

#for sample in 2016 2017 2018
#do
#    grep MET$sample skimFilesSR_data.txt > MET_Run$sample.txt
#done