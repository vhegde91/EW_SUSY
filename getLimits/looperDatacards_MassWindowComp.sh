#!/bin/sh
selType="_CompMass_"
rm dataCards/*.txt
#for i in TChiWH_300_1_MC2018 TChiWH_600_100_MC2018 TChiWH_600_200_MC2018 TChiWH_600_400_MC2018 TChiWH_800_100_MC2018 TChiWH_800_200_MC2018 TChiWH_800_400_MC2018 TChiWH_800_600_MC2018 TChiWH_1000_100_MC2018 TChiWW_800_100_MC2018 TChiWZ_1000_100_MC2018 TChiWZ_1000_1_MC2018 TChiWZ_300_1_MC2018 TChiWZ_600_100_MC2018 TChiWZ_600_1_MC2018 TChiWZ_800_100_MC2018 TChiWZ_800_1_MC2018 TChiWZ_800_400_MC2018
for i in TChiWZ_800_100_MC2018 TChiWH_800_100_MC2018
do
    for h in METvBin_FBWZ METvBin_FBWH METvBin_FBZ METvBin_FBH
    do
    	root -l 'makeDatacard_Simple.C("'$i'","'$h'")' -q
    	combineCards.py dataCards/$i*$h*.txt > dataCards/combined_${i}${selType}$h.txt
    	rm dataCards/$i*$h*.txt
	nlspM=$(echo $i | grep -o -E '[0-9]+' |sed -n '1p')
	lspM=$(echo $i | grep -o -E '[0-9]+' |sed -n '2p')
	mass=$(echo "$nlspM+$lspM/10000" |bc -l)
	echo $mass
    	combine -M AsymptoticLimits dataCards/combined_${i}${selType}${h}.txt -n $selType${i}_${h} -m $mass -t -1 
    done 
    echo "=================================================="
    # >=1 b-tagged
    echo $i "bTag"
    combineCards.py dataCards/combined_${i}${selType}*METvBin_FBWH.txt dataCards/combined_${i}${selType}*METvBin_FBH.txt > dataCards/combined_${i}${selType}_75To135.txt
    combineCards.py dataCards/combined_${i}${selType}*METvBin_FBWZ.txt dataCards/combined_${i}${selType}*METvBin_FBZ.txt > dataCards/combined_${i}${selType}_65To135.txt
    echo "-------- 75To135 -----------"${i}
    combine -M AsymptoticLimits dataCards/combined_${i}${selType}_75To135.txt -n $selType${i}_75To135 -m $mass -t -1
    echo "-------- 65To135 -----------"${i}
    combine -M AsymptoticLimits dataCards/combined_${i}${selType}_65To135.txt -n $selType${i}_65To135 -m $mass -t -1
done

for i in TChiWZ TChiWH
do
    ls -1 higgsCombine${selType}${i}_*_MC2018_75To135*.root > runList_${i}_75To135_${selType}.txt
    ./plotlimit runList_${i}_75To135_${selType}.txt 75To135${selType}${i}.root a

    ls -1 higgsCombine${selType}${i}_*_MC2018_65To135*.root > runList_${i}_65To135_${selType}.txt
    ./plotlimit runList_${i}_65To135_${selType}.txt 65To135${selType}${i}.root a
done

