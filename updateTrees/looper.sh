#!/bin/sh
for j in 800
do
    for i in 500 600
    do
	if [ "$i" -gt "$j" ] || [ "$i" -eq "$j" ]
	then
	    break
	fi
#	echo "root://cmseos.fnal.gov//store/user/vhegde/RA2bTrees_V17/EW_SUSY/Autumn18Fast.TChiWZ_NLSP${j}_LSP${i}_Autumn18_Private_0_RA2AnalysisTree.root" > smallrunList.txt
	echo "root://cmseos.fnal.gov//store/user/vhegde/RA2bTrees_V17/EW_SUSY/Autumn18Fast.TChiWH_${j}_${i}_Autumn18_Private_0_RA2AnalysisTree.root" > smallrunList.txt
##	echo "root://kodiak-se.baylor.edu//store/user/akanugan/TChipmWW/TChiWW_${j}_${i}_Autumn18_RA2AnalysisTree.root" > smallrunList.txt
	./signalReg smallrunList.txt Autumn18Fast.TChiWH_NLSP${j}_LSP${i}_Autumn18_Private_0_RA2AnalysisTree.root TChiWH_${j}
	root -l -q 'checkTree.C("Autumn18Fast.TChiWH_NLSP'${j}'_LSP'${i}'_Autumn18_Private_0_RA2AnalysisTree.root")'
	mv Autumn18Fast.TChiWH_NLSP${j}_LSP${i}_Autumn18_Private_0_RA2AnalysisTree.root /eos/uscms/store/user/vhegde/RA2bTrees_V17/EW_SUSY/
    done
done