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
#include "MT2_ROOT.h"

class SignalReg : public NtupleVariables{

 public:
  SignalReg(const TString &inputFileList="foo.txt", const char *outFileName="histo.root",const char *dataset="data");
  ~SignalReg();
  Bool_t   FillChain(TChain *chain, const TString &inputFileList);
  Long64_t LoadTree(Long64_t entry);
  void     EventLoop(const char *,const char *);
  void     BookHistogram(const char *);
  int  getEventType();
  void getEventTypeFine();
  TString getEventTypeWH();
  void print(Long64_t);
  //
  TLorentzVector visa = TLorentzVector( -18.1222 , -14.4356 , 0 , 158.653);
  TLorentzVector visb = TLorentzVector( 48.2681 , 38.449 , 0 , 62.513);
  TLorentzVector met4 = TLorentzVector( -30.1459 , -24.0134 , 0 , 74.8509);
  //  ComputeMT2 mycalc;
  //
  //Variables defined
  bool isMC=true;
  double wt=0,lumiInfb=35.815165;
  double massLow = 65., massHigh = 90.; //65-90, 55-100
  double massLowH = 85., massHighH = 135.; //85-135, 75-145
  double deepDoubleBDiscriminatorValue = 0.45; //0.3 for DoubleBDiscriminatorValue
  double deepCSVvalue = 0;
  vector<TLorentzVector> bjets;
  vector<TLorentzVector> nonbjets;
  //    double massLow = 0., massHigh = 1000.;

  vector<float> tau21;
  float i_tau21;
  TLorentzVector bestAK8J1, bestAK8J2;
  int bestAK8J1IdxInMainColl = -1, bestAK8J2IdxInMainColl = -1;
  vector<TLorentzVector> goodAk8;
  vector<double> goodAk8Mass;
  vector<int> ak8IdxInMainColl;

  vector<double> METvbinsZZMET={300,450,600,800,1000,1200,2000};
  vector<double> METvbins={200,250,300,350,400,450,550,650,800,1200};
  vector<double> mTvbins={400,500,600,700,800,900,1050,1300,2000};
  vector<double> mT2Jvbins={0,100,200,300,400,500,600,800,1000,1400};
  vector<double> mTSumvbins={500,700,900,1100,1400,1700,2000,2500};
  vector<double> MT2vbins={0,50,100,150,200,250,300,350,400,450,550,650,800,1200};

  TH1D *h_filters;
  TH1D *h_EvtType;
  TH1D *h_EvtTypeFine;
  TH1D *h_EvtTypeWH;
  TH1D *h_EvtTypeWH_ak84;
  TH1D *h_MET;
  TH1D *h_METvBin, *h_METvBinZZMET;
  TH1D *h_MHT;
  TH1D *h_HT;
  TH1D *h_NJets;
  TH1D *h_BTags;
  TH1D *h_LeadbPairMass;
  TH1D *h_LeadNonbPairMass;
  TH1D *h_DeepWdiscr;
  TH1D *h_DeepDoubleBdiscr;
  TH1D *h_LeadAK8Mass;
  TH1D *h_AK8MassNearGenH;
  TH1D *h_AK8MassNearGenW;
  TH1D *h_AK8MassNearGenZ;

  TH1D *h_MT, *h_MTvBin, *h_MT2, *h_MT2vBin;
  TH1D *h_MT2J, *h_MT2JvBin;
  TH1D *h_mTRatio, *h_mTSum, *h_mTSumvBin;
  TH1D *h_mEfft;
  TH1D *h_mTbMin, *h_mCT;
  TH1D *h_dPhiMETAK8;
  TH1D *h_dPhiAK8J1J2;
  TH1D *h_RA2bBins;
  TH1D *h_dRbosons;
  TH1D *h_nAk4jNotAK8;
  TH1D *h_AK4BosonCandMass;
  TH1D *h_dRaK4BosonGenBoson;
  TH1D *h_dRbestAK8AK4Cand;
  TH1D *h_dPhiBestAK8AK4Cand;
  TH1D *h_deepDoubleBdiscrForHcand;
  TH1D *h_TaggedHMassFailM;
  TH1D *h_TaggedWMassFailM;
  TH1D *h_dRTaggedHFailM;
  TH1D *h_dRTaggedWFailM;

  TH1D *h_dPhibJetAK4BosonCand;
  TH1D *h_dRbJetAK4BosonCand;

  TH1D *h_nAK8;
  TH1D *h_AK8J1Pt, *h_AK8J1Mass, *h_AK8J1Eta, *h_AK8J1Tau21;
  TH1D *h_AK8J2Pt, *h_AK8J2Mass, *h_AK8J2Eta, *h_AK8J2Tau21;

  TH1D *h_InvMassAK8Jets,*h_dPhibJetMET,*h_dPhibJetAK8,*h_dPhibJetAK8J2;
  TH1D *h_dRbJetAK8, *h_dRbJetAK8J2, *h_AK8J1J2MassRatio;

  TH2D *h2_mTbMin_mCT;
  TH2D *h2_AK8J1J2Tau21;
  TH2D *h2_AK8J1J2Mass;
  TH2D *h2_SusyPDGMass;
  TH2D *h2_LSPPtBosonPt;
  TH2D *h2_bJetPos;

  TH1D *h_dPhi1;
  TH1D *h_dPhi2;
  TH1D *h_dPhi3;
  TH1D *h_dPhi4;

  TH1D *h_EvtTypeWH_0AK8M;
  TH1D *h_MET_catWH[36];
  TH1D *h_METvBin_catWH[36];
  TH1D *h_mT_catWH[36];
  TH1D *h_AK8HPt_catWH[36];
  TH1D *h_AK8HEta_catWH[36];
  TH1D *h_AK8HM_catWH[36];
  TH1D *h_dPhiMETAK8H_catWH[36];
  TH1D *h_AK8WPt_catWH[36];
  TH1D *h_AK8WEta_catWH[36];
  TH1D *h_AK8WM_catWH[36];
  TH1D *h_mT2J_catWH[36];

  TH1D *h_LeadbPairMass_catWH[36];
  TH1D *h_LeadNonbPairMass_catWH[36];

  TH1D *h_mTLead_catWH[36];
  TH1D *h_mT2Lead_catWH[36];

  TH1D *h_mEfft_catWH[36];
  TH1D *h_mTbMin_catWH[36];
  TH1D *h_mCT_catWH[36];
  TH1D *h_mTqMin_catWH[36];
  TH1D *h_mCTq_catWH[36];

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
  h_EvtType = new TH1D("EvtType","Event type",10,0,10);
  h_EvtTypeFine = new TH1D("EvtTypeFine","Event type fine event category",20,0,20);
  h_EvtTypeWH = new TH1D("EvtTypeWH","Event type for WH event category",16,0,16);
  h_EvtTypeWH_ak84 = new TH1D("EvtTypeWH_ak84","Event type for WH event category, if W/H AK8 mass fail, use AK4 pair mass",12,0,12);
  h_EvtTypeWH_0AK8M = new TH1D("EvtTypeWH_0AK8M","Event type for WH event category, no AK8 jet within W/H mass window",10,0,10);

  h_MET = new TH1D("MET","MET",200,0,2000);
  h_MHT = new TH1D("MHT","MHT",200,0,2000);
  h_METvBinZZMET = new TH1D("METvBinZZMET","MET variable bins",METvbinsZZMET.size()-1,&(METvbinsZZMET[0]));
  h_METvBin = new TH1D("METvBin","MET variable bins",METvbins.size()-1,&(METvbins[0]));
  h_HT = new TH1D("HT","HT",100,0,5000);
  h_NJets = new TH1D("NJets","NJets with pT > 30, |eta| < 20.4",20,0,20);  
  h_BTags = new TH1D("BTags","BTags with DeepCSV MedWP",10,0,10);  

  h_LeadNonbPairMass = new TH1D("LeadNonbPairMass","M(Lead 2 Non-bjets)",60,0,300);
  h_LeadbPairMass = new TH1D("LeadbPairMass","M(Lead 2 bjets)",60,0,300);

  h_DeepWdiscr = new TH1D("DeepWdiscr","Deep W dicriminator for leading AK8",100,0,1);
  h_DeepDoubleBdiscr = new TH1D("DeepDoubleBdiscr","deepDoubleBdicsr for leading AK8",100,0,1);
  h_LeadAK8Mass = new TH1D("LeadAK8Mass","Leading AK8 SD Mass",60,0,300);
  h_AK8MassNearGenH = new TH1D("AK8MassNearGenH","AK8 SD mass near GenH (dR<0.3)",60,0,300);;
  h_AK8MassNearGenW = new TH1D("AK8MassNearGenW","AK8 SD mass near GenW (dR<0.3)",60,0,300);
  h_AK8MassNearGenZ = new TH1D("AK8MassNearGenZ","AK8 SD mass near GenZ (dR<0.3)",60,0,300);

  h_dRbosons = new TH1D("dRbosons","#DeltaR b/w Gen W/Z/H",100,0,5);
  h_nAk4jNotAK8 = new TH1D("nAk4jNotAK8","No. of AK4 jets not within 0.8 of best AK8 jet or second best, if exists, AK8 jet",8,0,8);
  h_AK4BosonCandMass = new TH1D("AK4BosonCandMass","Boson->qq Candidate mass, calc from AK4",60,0,300);
  h_dRaK4BosonGenBoson = new TH1D("dRaK4BosonGenBoson","Min dR(aK4BosonCandidate, GenBoson)",100,0,5);
  h_dRbestAK8AK4Cand = new TH1D("dRbestAK8AK4Cand","dR(Best AK8, AK4 candidate for W/Z) for 1 boosted AK8 events",100,0,5);
  h_dPhiBestAK8AK4Cand = new TH1D("dPhiBestAK8AK4Cand","dPhi(Best AK8, AK4 candidate for W/Z) for 1 boosted AK8 events",40,0,4);

  h_dPhibJetAK4BosonCand = new TH1D("dPhibJetAK4BosonCand","dPhi(leading b, AK4 pair candidate for W/Z) for 1 boosted AK8 events",40,0,4);
  h_dRbJetAK4BosonCand = new TH1D("dRJetAK4BosonCand","dR(leading b, AK4 pair candidate for W/Z) for 1 boosted AK8 events",100,0,5);
  h_TaggedHMassFailM = new TH1D("TaggedHMassFailM","SD mass of AK8 tagged as H, but not within H-mass",200,0,200);
  h_dRTaggedHFailM = new TH1D("dRTaggedHFailM","dR(H tagged AK8 and failing H mass,Gen H)",60,0,3);
  h_TaggedWMassFailM = new TH1D("TaggedWMassFailM","SD mass of AK8 tagged as W, but not within W-mass",200,0,200);
  h_dRTaggedWFailM = new TH1D("dRTaggedWFailM","dR(W tagged AK8 and failing W mass,Gen W)",60,0,3);

  h_MT2 = new TH1D("MT2","MT2(AKJ1, AK8J2)",200,0,2000);
  h_MT = new TH1D("mT","mT(MET,AK8J)",200,0,2000);
  h_MT2J = new TH1D("mT2J","mT(MET,AK8J2)",200,0,2000);  
  h_mTSum = new TH1D("mTSum","mT+mT2J",200,0,4000);
  h_MTvBin = new TH1D("mTvBin","mT(MET,AK8J)",mTvbins.size()-1,&(mTvbins[0]));
  h_MT2JvBin = new TH1D("mT2JvBin","mT(MET,AK8J2)",mT2Jvbins.size()-1,&(mT2Jvbins[0]));
  h_mTSumvBin = new TH1D("mTSumvBin","mT+mT2J",mTSumvbins.size()-1,&(mTSumvbins[0]));
  h_MT2vBin = new TH1D("MT2vBin","MT2(AKJ1, AK8J2) variable bins",MT2vbins.size()-1,&(MT2vbins[0]));

  h_mTbMin = new TH1D("mTbMin","min(mTb1,mTb2)",200,0,2000);
  h_mCT = new TH1D("mCT","mCT(b1,b2)",200,0,2000);

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
  h_deepDoubleBdiscrForHcand = new TH1D("deepDoubleBdiscrForHcand","deepDoubleBdiscriminator value for AK8 jet matched (dR < 0.3) to Gen H",100,-1,1);

  h2_mTbMin_mCT = new TH2D("mTbMin_mCT","x:min(mTb1,mTb2), mCT(b1,b2)",200,0,2000,200,0,2000);
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

  h_EvtType->Fill("2 Boosted",0);
  h_EvtType->Fill("1 Boosted",0);
  h_EvtType->Fill("1 Good AK8",0);
  h_EvtType->Fill("2 Good AK8",0);
  h_EvtType->Fill("0 Good AK8",0);

  TString catName;
  for(int t=2;t>=0;t--)
    for(int M=2;M>=0;M--){
      catName = to_string(t)+"T"+to_string(M)+"M";
      h_EvtTypeFine->Fill(catName,0);
    }
  for(int t=1;t>=0;t--)
    for(int M=1;M>=0;M--)
      for(int m=1;m>=0;m--){
	catName = to_string(t)+"T"+to_string(M)+"M"+to_string(m)+"m";
	h_EvtTypeFine->Fill(catName,0);
      }
  //----
  for(int t=1;t>=0;t--)
    for(int M=1;M>=0;M--)
      for(int th=1;th>=0;th--)
	for(int Mh=1;Mh>=0;Mh--){
	  catName = to_string(t)+"Wt"+to_string(M)+"Wm"+to_string(th)+"Ht"+to_string(Mh)+"Hm";
	  h_EvtTypeWH->Fill(catName,0);
	}
  //-----
  for(int j=1;j>=0;j--)
    for(int W=1;W>=0;W--)
      for(int H=1;H>=0;H--){
	catName = to_string(j)+"-"+to_string(W)+"Wm"+to_string(H)+"Hm";
	h_EvtTypeWH_0AK8M->Fill(catName,0);
      }
  //-----
  for(int t=1;t>=0;t--)
    for(int M=1;M>=0;M--)
      for(int th=1;th>=0;th--)
	for(int Mh=1;Mh>=0;Mh--){
	  if(M==1 && Mh==1) continue;
	  catName = to_string(t)+"Wt"+to_string(M)+"wm"+to_string(th)+"Ht"+to_string(Mh)+"hm";
	  h_EvtTypeWH_ak84->Fill(catName,0);
	}
  //----
  TDirectory *dir_WH = oFile->mkdir("WH","WH hist");
  dir_WH->cd();
  int iHist = 0;
  for(int t=1;t>=0;t--)
    for(int M=1;M>=0;M--)
      for(int th=1;th>=0;th--)
	for(int Mh=1;Mh>=0;Mh--){
	  catName = to_string(t)+"Wt"+to_string(M)+"Wm"+to_string(th)+"Ht"+to_string(Mh)+"Hm";
	  h_MET_catWH[iHist] = new TH1D("MET_"+catName,"MET for "+catName,200,0,2000);
	  h_METvBin_catWH[iHist] = new TH1D("METvBin_"+catName,"MET for "+catName,METvbins.size()-1,&(METvbins[0]));
	  h_mT_catWH[iHist] = new TH1D("mT_"+catName,"mT for "+catName,200,0,2000);
	  h_AK8HPt_catWH[iHist] = new TH1D("AK8HPt_"+catName,"H cand AK8 Pt for "+catName,200,0,2000);
	  h_AK8HEta_catWH[iHist] = new TH1D("AK8HEta_"+catName,"H cand AK8 Eta for "+catName,120,-6,6);
	  h_AK8HM_catWH[iHist] = new TH1D("AK8HM_"+catName,"H cand AK8 Mass for "+catName,40,0,200);
	  h_dPhiMETAK8H_catWH[iHist] = new TH1D("dPhiMETAK8H_"+catName,"DPhi(H cand AK8,MET) for "+catName,80,0,4);
	  h_AK8WPt_catWH[iHist] = new TH1D("AK8WPt_"+catName,"W cand AK8 Pt for "+catName,200,0,2000);
	  h_AK8WEta_catWH[iHist] = new TH1D("AK8WEta_"+catName,"W cand AK8 Eta for "+catName,120,-6,6);
	  h_AK8WM_catWH[iHist] = new TH1D("AK8WM_"+catName,"W cand AK8 Mass for "+catName,40,0,200);
	  h_mT2J_catWH[iHist] = new TH1D("mT2J_"+catName,"mT(2ndAK8,MET) for "+catName,200,0,2000);

	  h_mTLead_catWH[iHist] = new TH1D("mTLead_"+catName,"mT(Leading AK8,MET) for "+catName,200,0,2000);
	  h_mT2Lead_catWH[iHist] = new TH1D("mT2Lead_"+catName,"mT(2nd leading AK8,MET) for "+catName,200,0,2000);

	  h_LeadNonbPairMass_catWH[iHist] = new TH1D("LeadNonbPairMass_"+catName,"M(Lead 2 Non-bjets) for "+catName,60,0,300);
	  h_LeadbPairMass_catWH[iHist] = new TH1D("LeadbPairMass_"+catName,"M(Lead 2 bjets) for "+catName,60,0,300);
	  
	  h_mEfft_catWH[iHist] = new TH1D("mEfft_"+catName,"HT + MET for "+catName,200,0,2000);
	  h_mTbMin_catWH[iHist] = new TH1D("mTbMin_"+catName,"min(mT_b1, mT_b2) for "+catName,100,0,1000);
	  h_mCT_catWH[iHist] = new TH1D("mCT_"+catName,"mCT(b1,b2) for "+catName,100,0,1000);
	  h_mTqMin_catWH[iHist] = new TH1D("mTqMin_"+catName,"min(mT_q1, mT_q2) from W for "+catName,100,0,1000);
	  h_mCTq_catWH[iHist] = new TH1D("mCTq_"+catName,"mCT(q1, q2) from W for ",100,0,1000);

	  iHist++;
	}
  for(int j=1;j>=0;j--)
    for(int W=1;W>=0;W--)
      for(int H=1;H>=0;H--){
	catName = to_string(j)+"-"+to_string(W)+"Wm"+to_string(H)+"Hm";
	h_MET_catWH[iHist] = new TH1D("MET_"+catName,"MET for "+catName,200,0,2000);
	h_METvBin_catWH[iHist] = new TH1D("METvBin_"+catName,"MET for "+catName,METvbins.size()-1,&(METvbins[0]));
	h_mT_catWH[iHist] = new TH1D("mT_"+catName,"mT for "+catName,200,0,2000);
	h_AK8HPt_catWH[iHist] = new TH1D("AK8HPt_"+catName,"H cand AK8 Pt for "+catName,200,0,2000);
	h_AK8HEta_catWH[iHist] = new TH1D("AK8HEta_"+catName,"H cand AK8 Eta for "+catName,120,-6,6);
	h_AK8HM_catWH[iHist] = new TH1D("AK8HM_"+catName,"H cand AK8 Mass for "+catName,40,0,200);
	h_dPhiMETAK8H_catWH[iHist] = new TH1D("dPhiMETAK8H_"+catName,"DPhi(H cand AK8,MET) for "+catName,80,0,4);
	h_AK8WPt_catWH[iHist] = new TH1D("AK8WPt_"+catName,"W cand AK8 Pt for "+catName,200,0,2000);
	h_AK8WEta_catWH[iHist] = new TH1D("AK8WEta_"+catName,"W cand AK8 Eta for "+catName,120,-6,6);
	h_AK8WM_catWH[iHist] = new TH1D("AK8WM_"+catName,"W cand AK8 Mass for "+catName,40,0,200);
	h_mT2J_catWH[iHist] = new TH1D("mT2J_"+catName,"mT(2ndAK8,MET) for "+catName,200,0,2000);

	h_mTLead_catWH[iHist] = new TH1D("mTLead_"+catName,"mT(Leading AK4 cand,MET) for "+catName,200,0,2000);
	h_mT2Lead_catWH[iHist] = new TH1D("mT2Lead_"+catName,"mT(2nd leading AK4 cand,MET) for "+catName,200,0,2000);
	  
	h_LeadNonbPairMass_catWH[iHist] = new TH1D("LeadNonbPairMass_"+catName,"M(Lead 2 Non-bjets) for "+catName,60,0,300);
	h_LeadbPairMass_catWH[iHist] = new TH1D("LeadbPairMass_"+catName,"M(Lead 2 bjets) for "+catName,60,0,300);
	
	h_mEfft_catWH[iHist] = new TH1D("mEfft_"+catName,"HT + MET for "+catName,200,0,2000);
	h_mTbMin_catWH[iHist] = new TH1D("mTbMin_"+catName,"min(mT_b1, mT_b2) for "+catName,100,0,1000);
	h_mCT_catWH[iHist] = new TH1D("mCT_"+catName,"mCT(b1,b2) for "+catName,100,0,1000);
	h_mTqMin_catWH[iHist] = new TH1D("mTqMin_"+catName,"min(mT_q1, mT_q2) from W for "+catName,100,0,1000);
	h_mCTq_catWH[iHist] = new TH1D("mCTq_"+catName,"mCT(q1, q2) from W for ",100,0,1000);
	iHist++;
      }
  for(int t=1;t>=0;t--)
    for(int M=1;M>=0;M--)
      for(int th=1;th>=0;th--)
	for(int Mh=1;Mh>=0;Mh--){
	  if(M==1 && Mh==1) continue;
	  catName = to_string(t)+"Wt"+to_string(M)+"wm"+to_string(th)+"Ht"+to_string(Mh)+"hm";
	  h_MET_catWH[iHist] = new TH1D("MET_"+catName,"MET for "+catName,200,0,2000);
	  h_METvBin_catWH[iHist] = new TH1D("METvBin_"+catName,"MET for "+catName,METvbins.size()-1,&(METvbins[0]));
	  h_mT_catWH[iHist] = new TH1D("mT_"+catName,"mT for "+catName,200,0,2000);
	  h_AK8HPt_catWH[iHist] = new TH1D("AK8HPt_"+catName,"H cand AK8/AK4 pair Pt for "+catName,200,0,2000);
	  h_AK8HEta_catWH[iHist] = new TH1D("AK8HEta_"+catName,"H cand AK8/AK4 pair Eta for "+catName,120,-6,6);
	  h_AK8HM_catWH[iHist] = new TH1D("AK8HM_"+catName,"H cand AK8/AK4 pair Mass for "+catName,40,0,200);
	  h_dPhiMETAK8H_catWH[iHist] = new TH1D("dPhiMETAK8H_"+catName,"DPhi(H cand AK8,MET) for "+catName,80,0,4);
	  h_AK8WPt_catWH[iHist] = new TH1D("AK8WPt_"+catName,"W cand AK8/AK4 pair Pt for "+catName,200,0,2000);
	  h_AK8WEta_catWH[iHist] = new TH1D("AK8WEta_"+catName,"W cand AK8/AK4 pair Eta for "+catName,120,-6,6);
	  h_AK8WM_catWH[iHist] = new TH1D("AK8WM_"+catName,"W cand AK8/AK4 pair Mass for "+catName,40,0,200);
	  h_mT2J_catWH[iHist] = new TH1D("mT2J_"+catName,"mT(2ndAK8,MET) for "+catName,200,0,2000);

	  h_mTLead_catWH[iHist] = new TH1D("mTLead_"+catName,"mT(Leading AK8,MET) for "+catName,200,0,2000);
	  h_mT2Lead_catWH[iHist] = new TH1D("mT2Lead_"+catName,"mT(2nd leading AK8,MET) for "+catName,200,0,2000);
	  
	  h_LeadNonbPairMass_catWH[iHist] = new TH1D("LeadNonbPairMass_"+catName,"M(Lead 2 Non-bjets) for "+catName,60,0,300);
	  h_LeadbPairMass_catWH[iHist] = new TH1D("LeadbPairMass_"+catName,"M(Lead 2 bjets) for "+catName,60,0,300);
	  
	  h_mEfft_catWH[iHist] = new TH1D("mEfft_"+catName,"HT + MET for "+catName,200,0,2000);
	  h_mTbMin_catWH[iHist] = new TH1D("mTbMin_"+catName,"min(mT_b1, mT_b2) for "+catName,100,0,1000);
	  h_mCT_catWH[iHist] = new TH1D("mCT_"+catName,"mCT(b1,b2) for "+catName,100,0,1000);
	  h_mTqMin_catWH[iHist] = new TH1D("mTqMin_"+catName,"min(mT_q1, mT_q2) from W for "+catName,100,0,1000);
	  h_mCTq_catWH[iHist] = new TH1D("mCTq_"+catName,"mCT(q1, q2) from W for ",100,0,1000);
	  
	  iHist++;
	}
  cout<<"iHist"<<iHist<<endl;
  
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

