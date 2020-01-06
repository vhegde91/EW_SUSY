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

  vector<double> METvbinsZZMET={300,450,600,800,1000,1200,2000};
  vector<double> METvbins={200,250,300,350,400,450,550,650,800,1200};
  vector<double> mTvbins={400,500,600,700,800,900,1050,1300,2000};
  vector<double> mT2Jvbins={0,100,200,300,400,500,600,800,1000,1400};
  vector<double> mTSumvbins={500,700,900,1100,1400,1700,2000,2500};

  TH1D *h_filters;
  TH1D *h_MET;
  TH1D *h_METvBin, *h_METvBinZZMET;
  TH1D *h_MHT;
  TH1D *h_HT;
  TH1D *h_NJets;
  TH1D *h_BTags;

  TH1D *h_MT, *h_MTvBin;
  TH1D *h_MT2J, *h_MT2JvBin;
  TH1D *h_mTRatio, *h_mTSum, *h_mTSumvBin;
  TH1D *h_dPhiMETAK8;
  TH1D *h_dPhiAK8J1J2;
  TH1D *h_RA2bBins;
  TH1D *h_dRbosons;
  TH1D *h_nAk4jNotAK8;
  TH1D *h_AK4BosonCandMass;
  TH1D *h_dRaK4BosonGenBoson;
  TH1D *h_dRbestAK8AK4Cand;
  TH1D *h_dPhiBestAK8AK4Cand;

  TH1D *h_dPhibJetAK4BosonCand;
  TH1D *h_dRbJetAK4BosonCand;

  TH1D *h_nAK8;
  TH1D *h_AK8J1Pt, *h_AK8J1Mass, *h_AK8J1Eta, *h_AK8J1Tau21;
  TH1D *h_AK8J2Pt, *h_AK8J2Mass, *h_AK8J2Eta, *h_AK8J2Tau21;

  TH1D *h_InvMassAK8Jets,*h_dPhibJetMET,*h_dPhibJetAK8,*h_dPhibJetAK8J2;
  TH1D *h_dRbJetAK8, *h_dRbJetAK8J2, *h_AK8J1J2MassRatio;

  TH2D *h2_AK8J1J2Tau21;
  TH2D *h2_AK8J1J2Mass;
  TH2D *h2_SusyPDGMass;
  TH2D *h2_LSPPtBosonPt;
  TH2D *h2_bJetPos;

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
  h_METvBinZZMET = new TH1D("METvBinZZMET","MET variable bins",METvbinsZZMET.size()-1,&(METvbinsZZMET[0]));
  h_METvBin = new TH1D("METvBin","MET variable bins",METvbins.size()-1,&(METvbins[0]));
  h_HT = new TH1D("HT","HT",100,0,5000);
  h_NJets = new TH1D("NJets","NJets with pT > 30, |eta| < 20.4",20,0,20);  
  h_BTags = new TH1D("BTags","BTags with DeepCSV MedWP",10,0,10);  
  h_dRbosons = new TH1D("dRbosons","#DeltaR b/w Gen W/Z/H",100,0,5);
  h_nAk4jNotAK8 = new TH1D("nAk4jNotAK8","No. of AK4 jets not within 0.8 of best AK8 jet or second best, if exists, AK8 jet",8,0,8);
  h_AK4BosonCandMass = new TH1D("AK4BosonCandMass","Boson->qq Candidate mass, calc from AK4",60,0,300);
  h_dRaK4BosonGenBoson = new TH1D("dRaK4BosonGenBoson","Min dR(aK4BosonCandidate, GenBoson)",100,0,5);
  h_dRbestAK8AK4Cand = new TH1D("dRbestAK8AK4Cand","dR(Best AK8, AK4 candidate for W/Z) for 1 boosted AK8 events",100,0,5);
  h_dPhiBestAK8AK4Cand = new TH1D("dPhiBestAK8AK4Cand","dPhi(Best AK8, AK4 candidate for W/Z) for 1 boosted AK8 events",40,0,4);

  h_dPhibJetAK4BosonCand = new TH1D("dPhibJetAK4BosonCand","dPhi(leading b, AK4 pair candidate for W/Z) for 1 boosted AK8 events",40,0,4);
  h_dRbJetAK4BosonCand = new TH1D("dRJetAK4BosonCand","dR(leading b, AK4 pair candidate for W/Z) for 1 boosted AK8 events",100,0,5);

  h_MT = new TH1D("mT","mT(MET,AK8J)",200,0,2000);
  h_MT2J = new TH1D("mT2J","mT(MET,AK8J2)",200,0,2000);  
  h_mTSum = new TH1D("mTSum","mT+mT2J",200,0,4000);
  h_MTvBin = new TH1D("mTvBin","mT(MET,AK8J)",mTvbins.size()-1,&(mTvbins[0]));
  h_MT2JvBin = new TH1D("mT2JvBin","mT(MET,AK8J2)",mT2Jvbins.size()-1,&(mT2Jvbins[0]));
  h_mTSumvBin = new TH1D("mTSumvBin","mT+mT2J",mTSumvbins.size()-1,&(mTSumvbins[0]));

  h_mTRatio = new TH1D("mTRatio","mT/mT2J: ratio of mT1 and mT2J",500,0,5);

  h_dPhiMETAK8 = new TH1D("dPhiMETAK8","dPhi(MET,AK8J)",40,0,4);
  h_dPhiAK8J1J2 = new TH1D("dPhiAK8J1J2","dPhi(AK8J1,AK8J2)",40,0,4);

  h_nAK8 = new TH1D("nAK8","no. of AK8 jets",5,0,5);
  h_AK8J1Pt = new TH1D("AK8Pt","Leading AK8 jets Pt",200,0,2000);
  h_AK8J1Eta = new TH1D("AK8Eta","AK8 Eta",120,-6,6);
  h_AK8J1Mass = new TH1D("AK8Mass","AK8 Mass",60,0,300);
  h_AK8J1Tau21 = new TH1D("AK8Tau21","AK8 Tau21",20,0,1);

  h_AK8J2Pt = new TH1D("AK8J2Pt","2nd leading AK8 jets Pt",200,0,2000);
  h_AK8J2Eta = new TH1D("AK8J2Eta","AK8J2 Eta",120,-6,6);
  h_AK8J2Mass = new TH1D("AK8J2Mass","AK8J2 Mass",60,0,300);
  h_AK8J2Tau21 = new TH1D("AK8J2Tau21","AK8J2 Tau21",20,0,1);

  h_InvMassAK8Jets = new TH1D("InvMassAK8Jets","Invariant mass of leading 2 AK8 jets",200,0,2000);
  h_AK8J1J2MassRatio = new TH1D("AK8J1J2MassRatio","AK8J2 mass/AK8J1 mass",100,0,5);
  h_dRbJetAK8 = new TH1D("dRbJetAK8","#DeltaR(AK8,b-jet)",200,0,2);
  h_dRbJetAK8J2 = new TH1D("dRbJetAK8J2","#DeltaR(AK8J2,b-jet)",200,0,2);

  h_dPhibJetMET = new TH1D("dPhibJetMET","#Delta#Phi(leading b, MET)",40,0,4);
  h_dPhibJetAK8 = new TH1D("dPhibJetAK8","#Delta#Phi(leading b, AK8)",40,0,4);
  h_dPhibJetAK8J2 = new TH1D("dPhibJetAK8J2","#Delta#Phi(leading b, AK8J2)",40,0,4);

  h2_AK8J1J2Tau21 = new TH2D("AK8J1J2Tau21","x:AK8J1 #tau21 vs y:AK8J2 #tau21",20,0,1,20,0,1);
  h2_AK8J1J2Mass = new TH2D("AK8J1J2Mass","x:AK8J1 Mass vs y:AK8J2 Mass",60,0,300,60,0,300);
  h2_SusyPDGMass = new TH2D("SusyPDGMass","x:PDG ID of SUSY particle, y:Mass of particle",50,1000000,1000050,300,0,3000);
  h2_LSPPtBosonPt = new TH2D("LSPPtBosonPt","x:LSP Pt, y:Boson Pt",200,0,2000,200,0,2000);
  h2_bJetPos = new TH2D("bJetPos","x:dR(leading bJet, AK8 jets), y: dR(leading bJet, AK8/AK4 Cand)",100,0,5,100,0,5);

  h_dPhi1 = new TH1D("DeltaPhi1","DeltaPhi1",40,0,4);
  h_dPhi2 = new TH1D("DeltaPhi2","DeltaPhi2",40,0,4);
  h_dPhi3 = new TH1D("DeltaPhi3","DeltaPhi3",40,0,4);
  h_dPhi4 = new TH1D("DeltaPhi4","DeltaPhi4",40,0,4);

  h_RA2bBins = new TH1D("RA2bBins","RA2b bins",175,0,175);
}

SignalReg::SignalReg(const TString &inputFileList, const char *outFileName, const char* dataset) {
  string nameData=dataset;
  TString nameData2 = nameData;
  TChain *tree = new TChain("tree");
  //  if(nameData2.Contains("TChiWZ")) tree = new TChain("TreeMaker2/PreSelection");
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

