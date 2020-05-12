#include<iostream>
#include<iomanip>
#include"TH1.h"
#include"TROOT.h"
#include"TH2.h"
#include"TFile.h"
#include"TDirectory.h"
#include"TTree.h"
#include"TBrowser.h"
#include"TF1.h"
#include<string>
#include<vector>
#include"TGraphErrors.h"
#include"TGraph.h"
#include"TLegend.h"
#include"TLatex.h"
#include"TCanvas.h"
#include"THStack.h"
#include"TStyle.h"


TString name;
vector<TString> fname={"TChiWH_800_100_MC2018.root",
		       "TChiWH_800_200_MC2018.root",
		       "TChiWH_800_400_MC2018.root",
		       "TChiWH_800_600_MC2018.root",
		       "TChiWH_600_400_MC2018.root",
		       "TChiWH_600_200_MC2018.root",
		       "TChiWH_600_100_MC2018.root",
		       "TChiWH_300_1_MC2018.root",
		       "ST__MC2018.root",
		       "Rare_MC2018MC2017MC2016.root",
		       "TTJets_MC2018.root",
		       "WJetsToLNu_HT_MC2018.root",
		       "ZJetsToNuNu_HT_MC2018.root"};
const int nfiles=13;    //Specify no. of files
TFile *f[nfiles];
TFile *fout[nfiles];
  
TString nameOfOutHist = "METvBin";
void merger(vector<TString>,TString);
void merger(vector<TString> hlist,TString name){
  TH1D *h0,*h1;
  TString name1;
  for(int i=0;i<nfiles;i++){
    f[i] = TFile::Open(fname[i]);
    name1 = fname[i];
    name1.ReplaceAll(".root","_"+name+".root");
    fout[i] = new TFile(name1,"recreate");
    for(int j=0;j<hlist.size();j++){
      if(j==0) f[i]->GetObject(hlist[j],h0);
      else{
	f[i]->GetObject(hlist[j],h1);
	h0->Add(h1);
      }
    }
    fout[i]->cd();
    h0->Write(nameOfOutHist);
    f[i]->Close();
    fout[i]->Close();
  }
}
void regionsCombineHists(){
  vector<TString> name1;
  //----- Fully boosted W&H
  name1.push_back("WH/METvBin_1Wt1Wm1Ht1Hm");//WH
  merger(name1,"FB_WH");
  
  //------ Fully W-tagged only
  name1.resize(0);
  name1.push_back("WH/METvBin_1Wt1Wm1Ht0Hm");//W
  name1.push_back("WH/METvBin_1Wt1Wm0Ht1Hm");//W
  name1.push_back("WH/METvBin_1Wt1Wm0Ht0Hm");//W
  merger(name1,"FB_W");

  //------ Fully H-tagged only
  name1.resize(0);
  name1.push_back("WH/METvBin_1Wt0Wm1Ht1Hm");//H
  name1.push_back("WH/METvBin_0Wt1Wm1Ht1Hm");//H
  name1.push_back("WH/METvBin_0Wt0Wm1Ht1Hm");//H
  merger(name1,"FB_H");
  
  //------- No AK8 passing W/H tag AND nass. Use AK4 pair for W and H.
  name1.resize(0);
  name1.push_back("WH/METvBin_1Wt0wm1Ht0hm");
  name1.push_back("WH/METvBin_1Wt0wm0Ht0hm");
  name1.push_back("WH/METvBin_0Wt1wm0Ht0hm");
  name1.push_back("WH/METvBin_0Wt0wm1Ht0hm");

  name1.push_back("WH/METvBin_0Wt1wm1Ht0hm");
  name1.push_back("WH/METvBin_0Wt0wm0Ht0hm");
  
  name1.push_back("WH/METvBin_0Wt1Wm0Ht1Hm");
  name1.push_back("WH/METvBin_1Wt0wm0Ht1hm");
  name1.push_back("WH/METvBin_0Wt0wm0Ht1hm");
  merger(name1,"leftOv");

  //------ Fully W-tagged, AK8/AK4 H-mass
  name1.resize(0);
  name1.push_back("WH/METvBin_1Wt1wm1Ht0hm");//Wh
  name1.push_back("WH/METvBin_1Wt1Wm0Ht1Hm");//Wh
  name1.push_back("WH/METvBin_1Wt1wm0Ht0hm");//Wh
  merger(name1,"FB_W_ak84H");

  //------ Fully H-tagged, AK8/AK4 W-mass
  name1.resize(0);
  name1.push_back("WH/METvBin_1Wt0wm1Ht1hm");//Hw
  name1.push_back("WH/METvBin_0Wt1Wm1Ht1Hm");//Hw
  name1.push_back("WH/METvBin_0Wt0wm1Ht1hm");//Hw
  merger(name1,"FB_H_ak84W");
  
  //----- using AK4 jets only. No AK8 info is used.
  //name1.push_back("WH/METvBin_0-1Wm1Hm");//hm
  //name1.push_back("WH/METvBin_0-1Wm0Hm");//w
  //name1.push_back("WH/METvBin_0-0Wm1Hm");//h
}
