#!/bin/sh
echo "Exp limits" > a.txt
for i in TChiWH_800_100_MC2018 TChiWH_800_200_MC2018 TChiWH_800_400_MC2018 TChiWH_800_600_MC2018 TChiWH_600_100_MC2018 TChiWH_600_200_MC2018 TChiWH_600_400_MC2018 TChiWH_300_1_MC2018
do
    for h in _FB_WH _FB_W _FB_H _leftOv _FB_H_ak84W _FB_W_ak84H
    do
	root -l 'makeDatacard_Simple.C("'$i'","'$h'")' -q
	combineCards.py dataCards/$i*$h*.txt > dataCards/combined_$i$h.txt
	rm dataCards/$i*.txt
    
    done 
    echo "=================================================="

    echo $i " FB_WH" |tee -a a.txt
    combine -M AsymptoticLimits dataCards/combined_${i}_FB_WH.txt -n ${i}_FB_WH -t -1 |tee -a a.txt
    echo $i " FB_W" |tee -a a.txt
    combine -M AsymptoticLimits dataCards/combined_${i}_FB_W.txt -n ${i}_FB_W -t -1 |tee -a a.txt
    echo $i " FB_H" |tee -a a.txt
    combine -M AsymptoticLimits dataCards/combined_${i}_FB_H.txt -n ${i}_FB_H -t -1 |tee -a a.txt
    echo $i " FB_WH_FB_W_FB_H" |tee -a a.txt
    combineCards.py dataCards/combined_${i}_FB_WH.txt dataCards/combined_${i}_FB_W.txt dataCards/combined_${i}_FB_H.txt > dataCards/combined_FB_WH_FB_W_FB_H$i.txt
    combine -M AsymptoticLimits dataCards/combined_FB_WH_FB_W_FB_H$i.txt -n ${i}_FB_WH_FB_W_FB_H -t -1 |tee -a a.txt

    echo $i " leftOv" |tee -a a.txt
    combine -M AsymptoticLimits dataCards/combined_${i}_leftOv.txt -n ${i}_leftOv -t -1 |tee -a a.txt
    combineCards.py dataCards/combined_${i}_leftOv.txt dataCards/combined_FB_WH_FB_W_FB_H$i.txt > dataCards/combined_FB_WH_FB_W_FB_H_leftOv$i.txt
    echo $i " FB_WH_FB_W_FB_H_leftOv" |tee -a a.txt
    combine -M AsymptoticLimits dataCards/combined_FB_WH_FB_W_FB_H_leftOv$i.txt -n ${i}_FB_WH_FB_W_FB_H_leftOv -t -1 |tee -a a.txt

    echo $i " _FB_H_ak84W" |tee -a a.txt
    combine -M AsymptoticLimits dataCards/combined_${i}_FB_H_ak84W.txt -n ${i}_FB_H_ak84W -t -1 |tee -a a.txt
    echo $i " _FB_W_ak84H" |tee -a a.txt
    combine -M AsymptoticLimits dataCards/combined_${i}_FB_W_ak84H.txt -n ${i}_FB_W_ak84H -t -1 |tee -a a.txt

    echo $i " FB_WH_FB_W_ak84H_FB_H_ak84W" |tee -a a.txt
    combineCards.py dataCards/combined_${i}_FB_WH.txt dataCards/combined_${i}_FB_W_ak84H.txt dataCards/combined_${i}_FB_H_ak84W.txt > dataCards/combined_FB_WH_FB_W_ak84H_FB_H_ak84W$i.txt
    combine -M AsymptoticLimits dataCards/combined_FB_WH_FB_W_ak84H_FB_H_ak84W$i.txt -n ${i}_FB_WH_FB_W_ak84H_FB_H_ak84W -t -1 |tee -a a.txt
    
done
grep 'TChi\|50.0' a.txt > opLimitCalc.txt
rm a.txt