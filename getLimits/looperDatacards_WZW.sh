#!/bin/sh
echo "Exp limits" > a.txt
for i in TChiWW_800_100_MC2018 TChiWZ_800_1_MC2018 TChiWZ_800_100_MC2018 TChiWZ_800_400_MC2018 TChiWZ_1000_1_MC2018 TChiWZ_1000_100_MC2018 TChiWZ_600_1_MC2018 TChiWZ_600_100_MC2018 TChiWZ_300_1_MC2018
do
#    for h in METvBin_2T2M METvBin_1T2M METvBin_1T1M
    for h in METvBin_Wmc_WWZ METvBin_Wmd_WWZ METvBin_Tau21_WWZ
#    for h in METvBin_Wmd_WWZ
    do
    	root -l 'makeDatacard_Simple.C("'$i'","","'$h'")' -q
    	combineCards.py dataCards/$i*$h*.txt > dataCards/combined_$i$h.txt
    	rm dataCards/$i*.txt
        echo $i |tee -a a.txt
    	combine -M AsymptoticLimits dataCards/combined_${i}${h}.txt -n ${i}_${h}.txt -t -1 |tee -a a.txt
    done 
    echo "=================================================="    
    # echo $i |tee -a a.txt
    # combineCards.py dataCards/combined_${i}*METvBin_2T2M*.txt dataCards/combined_${i}*METvBin_1T2M*.txt dataCards/combined_${i}*METvBin_1T1M*.txt > dataCards/combined_${i}_2T2M_1T2M_1T1M.txt
    # combine -M AsymptoticLimits dataCards/combined_${i}_2T2M_1T2M_1T1M.txt -t -1 |tee -a a.txt
#    combineCards.py dataCards/combined_${i}*METvBin_2T2M.txt dataCards/combined_${i}*METvBin_1T2M.txt > dataCards/combined_${i}_2T2M_1T2M.txt
#    combine -M AsymptoticLimits dataCards/combined_${i}_2T2M_1T2M.txt -t -1 |tee -a a.txt
done
grep 'TChi\|50.0' a.txt > opLimitCalc.txt
rm a.txt