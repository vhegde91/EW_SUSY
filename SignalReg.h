#ifndef SignalReg_H
#define SignalReg_H

#include <iostream>
#include <fstream>
#include <cmath>
#include <vector>
#include "NtupleVariables.h"
#include "TH1F.h"
#include "TH2.h"
#include <TProfile.h>
#include "TFile.h"
#include "TLorentzVector.h"
#include "TDirectory.h"

class SignalReg : public NtupleVariables{

 public:
  SignalReg(const TString &inputFileList="foo.txt", const char *outFileName="histo.root",const char *dataset="data");
  ~SignalReg();
  Bool_t   FillChain(TChain *chain, const TString &inputFileList);
  Long64_t LoadTree(Long64_t entry);
  void     EventLoop(const char *,const char *);
  void     BookHistogram(const char *);
  void print(Long64_t);

  //Variables defined
  bool isMC=true;
  double wt=0,lumiInfb=35.815165;

  TH1D *h_filters;
  TH1D *h_MET;
  TH1D *h_MHT;
  TH1D *h_HT;
  TH1D *h_NJets;
  TH1D *h_BTags;

  TH1D *h_AK8J1Pt, *h_AK8J1Mass, *h_AK8J1Eta, *h_AK8J1Tau21;

  TH1D *h_dPhi1;
  TH1D *h_dPhi2;
  TH1D *h_dPhi3;
  TH1D *h_dPhi4;

  TH1F *h_cutflow;
  TFile *oFile;
  
};
#endif

#ifdef SignalReg_cxx

void SignalReg::BookHistogram(const char *outFileName) {

  //  char hname[200], htit[200];
  //  double xlow = 0.0,  xhigh = 2000.0;
  //  int nbins = 2000;
  TString name,title;
 
  oFile = new TFile(outFileName, "recreate");
  TH1::SetDefaultSumw2(1);

  h_cutflow = new TH1F("CutFlow","cut flow",25,0,25);
  h_filters = new TH1D("Filters","Filters: Bin1 : all nEvnts, other bins: filter pass/fail",10,0,10);
  
  h_MET = new TH1D("MET","MET",200,0,2000);
  h_MHT = new TH1D("MHT","MHT",200,0,2000);
  h_HT = new TH1D("HT","HT",100,0,5000);
  h_NJets = new TH1D("NJets","NJets with pT > 30, |eta| < 20.4",20,0,20);  
  h_BTags = new TH1D("BTags","BTags with DeepCSV MedWP",10,0,10);  

  h_AK8J1Pt = new TH1D("AK8Pt","Leading AK8 jets Pt",200,0,2000);
  h_AK8J1Eta = new TH1D("AK8Eta","AK8Eta",120,-6,6);
  h_AK8J1Mass = new TH1D("AK8Mass","AK8Mass",60,0,300);
  h_AK8J1Tau21 = new TH1D("AK8Tau21","Tau21",20,0,1);
  
  h_dPhi1 = new TH1D("DeltaPhi1","DeltaPhi1",40,0,4);
  h_dPhi2 = new TH1D("DeltaPhi2","DeltaPhi2",40,0,4);
  h_dPhi3 = new TH1D("DeltaPhi3","DeltaPhi3",40,0,4);
  h_dPhi4 = new TH1D("DeltaPhi4","DeltaPhi4",40,0,4);

}

SignalReg::SignalReg(const TString &inputFileList, const char *outFileName, const char* dataset) {
  string nameData=dataset;
  TString nameData2 = nameData;
  TChain *tree = new TChain("tree");
  if(nameData2.Contains("TChiWZ")) tree = new TChain("TreeMaker2/PreSelection");
  if( ! FillChain(tree, inputFileList) ) {
    std::cerr << "Cannot get the tree " << std::endl;
  } else {
    std::cout << "Initiating analysis of dataset " << dataset << std::endl;
  }

  cout<<"Treating the input files as "<<nameData<<" for setting tree branches"<<endl;
  NtupleVariables::Init(tree,nameData);

  BookHistogram(outFileName);
  
}

Bool_t SignalReg::FillChain(TChain *chain, const TString &inputFileList) {
  int itr=0;
  TFile *filePointer;
  ifstream infile(inputFileList, ifstream::in);
  std::string buffer;

  if(!infile.is_open()) {
    std::cerr << "** ERROR: Can't open '" << inputFileList << "' for input" << std::endl;
    return kFALSE;
  }

  std::cout << "TreeUtilities : FillChain " << std::endl;
  while(1) {
    infile >> buffer;
    if(!infile.good()) break;
    //    std::cout << "Adding tree from " << buffer.c_str() << std::endl;
    chain->Add(buffer.c_str());
  }
  std::cout << "No. of Entries in this tree : " << chain->GetEntries() << std::endl;
  return kTRUE;
}

Long64_t SignalReg::LoadTree(Long64_t entry) {
  // Set the environment to read one entry                                                                                          
  if (!fChain) return -5;
  Long64_t centry = fChain->LoadTree(entry);
  if (centry < 0) return centry;
  if (!fChain->InheritsFrom(TChain::Class()))  return centry;
  TChain *chain = (TChain*)fChain;
  if (chain->GetTreeNumber() != fCurrent) {
    fCurrent = chain->GetTreeNumber();
    //    Notify();
  }
  return centry;
}

SignalReg::~SignalReg() { 

  if (!fChain) return;
  delete fChain->GetCurrentFile();
  oFile->cd();
  oFile->Write();
  oFile->Close();

}

#endif

