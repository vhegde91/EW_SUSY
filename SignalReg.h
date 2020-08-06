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
#include "TF1.h"
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
  void getEventTypeFine(int);
  void getEventTypeWZW();
  TString getEventTypeWH();
  void print(Long64_t);
  void applySDmassCorrAllAK8();
  bool passHEMjetVeto(double);
  //
  TLorentzVector visa = TLorentzVector( -18.1222 , -14.4356 , 0 , 158.653);
  TLorentzVector visb = TLorentzVector( 48.2681 , 38.449 , 0 , 62.513);
  TLorentzVector met4 = TLorentzVector( -30.1459 , -24.0134 , 0 , 74.8509);
  //  ComputeMT2 mycalc;
  //
  TFile *sdCorrFile;
  TF1 *puppisd_corrGEN, *puppisd_corrRECO_cen, *puppisd_corrRECO_for;
  TH1D *CRweightsHist;
  TString catname2j;
  vector<TLorentzVector> uncorrAK8SubjSum;
  vector<double> SDmassCorrFac;
  //Variables defined
  bool isMC=true;
  double wt=0,lumiInfb=35.815165;
  double massLow = 65., massHigh = 105.; //65-90
  double massLowH = 75., massHighH = 135.; //85-135, 75-145
  double massLowZ = 75., massHighZ = 135.;
  double deepDoubleBDiscriminatorValue = 0.7; //0.3 for DoubleBDiscriminatorValue
  double dwdisValue = 0.918;
  double dzdisValue = 0.918;
  double dwdisMDvalue = 0.479;//0.245;//0.479;//0.704;
  double dzhdisMDvalue = 0.175;
  double tau21Value = 0.35;
  double tau21LPValue = 0.45;
  double tau21CRValue = 0.45;
  double deepCSVvalue = 0;
  vector<TLorentzVector> bjets;
  vector<TLorentzVector> nonbjets;

  vector<float> mainWdisc;
  float i_mainWdisc;
  TLorentzVector bestAK8J1, bestAK8J2;
  int bestAK8J1IdxInMainColl = -1, bestAK8J2IdxInMainColl = -1;
  vector<TLorentzVector> goodAk8;
  vector<double> goodAk8Mass;
  vector<int> ak8IdxInMainColl;

  vector<double> METvbinsZZMET={300,450,600,800,1000,1200,2000};
  vector<double> METvbins={200,250,300,350,400,450,550,650,800,1200};
  vector<double> ptvbins={200,300,400,500,600,800,1200,2000};
  vector<double> mTvbins={400,500,600,700,800,900,1050,1300,2000};
  vector<double> mT2Jvbins={0,100,200,300,400,500,600,800,1000,1400};
  vector<double> mTSumvbins={500,700,900,1100,1400,1700,2000,2500};
  vector<double> MT2vbins={0,50,100,150,200,250,300,350,400,450,550,650,800,1200};

  TH1D *h_filters;
  TH1D *h_RegionCat;
  TH1D *h_EvtTypeFine;
  TH1D *h_EvtTypeWH;
  TH1D *h_EvtTypeWH_ak84;
  TH1D *h_catname2j;
  TH1D *h_MET;
  TH1D *h_METvBin, *h_METvBinZZMET;
  TH1D *h_MHT;
  TH1D *h_HT;
  TH1D *h_NJets;
  TH1D *h_BTags;
  TH1D *h_LeadbPairMass;
  TH1D *h_LeadNonbPairMass;
  TH1D *h_DeepWdiscr;
  TH1D *h_DeepZdiscr;
  TH1D *h_DeepDoubleBdiscr;
  TH1D *h_DeepWdiscrMD;
  TH1D *h_DeepZHdiscrMD;
  TH1D *h_LeadAK8Mass;
  TH1D *h_AK8MassNearGenH;
  TH1D *h_AK8MassNearGenW;
  TH1D *h_AK8MassNearGenZ;
  TH1D *h_NVtx;

  TH1D *h_DeepWdisc_BG;
  TH1D *h_DeepWdiscMD_BG;
  TH1D *h_Tau21_BG;
  TH1D *h_AK8Pt_BG;

  TH1D *h_DeepWdisc_SigW;
  TH1D *h_DeepWdiscMD_SigW;
  TH1D *h_Tau21_SigW;
  TH1D *h_AK8Pt_SigW;
  TH1D *h_AK8M_SigW;
  TH1D *h_GenWmassDeepWpass;
  TH1D *h_SDmassDeepWpassGenWmatch;
  TH1D *h_DeepWdisc_SigZ;
  TH1D *h_DeepWdiscMD_SigZ;
  TH1D *h_Tau21_SigZ;
  TH1D *h_AK8Pt_SigZ;
  TH1D *h_AK8M_SigZ;
  TH1D *h_GenZmassDeepWpass;
  TH1D *h_SDmassDeepWpassGenZmatch;
  TH1D *h_dRGenWb;

  TH1D *h_MT, *h_MTvBin;
  TH1D *h_MT2J, *h_MT2JvBin;
  TH1D *h_mTRatio, *h_mTSum, *h_mTSumvBin;
  TH1D *h_mEfft;
  TH1D *h_mTbMin, *h_mCT;
  TH1D *h_dPhiMETAK8;
  TH1D *h_dPhiAK8J1J2;
  TH1D *h_RA2bBins;
  TH1D *h_dRbosons;
  TH1D *h_deepDoubleBdiscrForHcand;
  TH1D *h_TaggedHMassFailM;
  TH1D *h_TaggedWMassFailM;
  TH1D *h_dRTaggedHFailM;
  TH1D *h_dRTaggedWFailM;

  TH1D *h_Zkids;
  TH1D *h_GenPartInDeepDoubleBpassAK8;
  TH1D *h_GenPartInDeepDoubleBfailAK8;
  
  TH1D *h_dPhiLSPs;

  TH1D *h_nAK8;
  TH1D *h_AK8J1Pt, *h_AK8J1Mass, *h_AK8J1Eta, *h_AK8J1MainWdisc;
  TH1D *h_AK8J2Pt, *h_AK8J2Mass, *h_AK8J2Eta, *h_AK8J2MainWdisc;
  TH1D *h_J1Pt, *h_J2Pt, *h_J1Eta, *h_J2Eta, *h_J1Phi, *h_J2Phi;

  TH1D *h_InvMassAK8Jets,*h_dPhibJetMET,*h_dPhibJetAK8,*h_dPhibJetAK8J2;
  TH1D *h_dRbJetAK8, *h_dRbJetAK8J2, *h_AK8J1J2MassRatio;

  TH2D *h2_mTbMin_mCT;
  TH2D *h2_LeadAK4EtaPhi;
  TH2D *h2_AK8J1J2Wdisc;
  TH2D *h2_AK8J1J2Mass;
  TH2D *h2_SusyPDGMass;
  TH2D *h2_LSPPtBosonPt;

  TH1D *h_dPhi1;
  TH1D *h_dPhi2;
  TH1D *h_dPhi3;
  TH1D *h_dPhi4;

  TH1D *h_mTLep;
  TH1D *h_LepType;
  TH1D *h_LepPt;
  TH1D *h_LepEta;
  TH1D *h_LepPhi;

  TH1D *h_WTagMassJPt, *h_WAntiTagMassJPt, *h_HTagMassJPt, *h_HAntiTagMassJPt;

  TH1D *h_METvBin_2T2M_SR, *h_MET_2T2M_SR, *h_mT_2T2M_SR;
  TH1D *h_METvBin_1T2M_SR, *h_MET_1T2M_SR, *h_mT_1T2M_SR;
  TH1D *h_METvBin_1T1M_SR, *h_MET_1T1M_SR, *h_mT_1T1M_SR;
  TH1D *h_METvBin_2T2M_CR, *h_MET_2T2M_CR, *h_mT_2T2M_CR;
  TH1D *h_METvBin_0T2M_CR, *h_MET_0T2M_CR, *h_mT_0T2M_CR;
  TH1D *h_METvBin_0T1M_CR, *h_MET_0T1M_CR, *h_mT_0T1M_CR;

  TH1D *h_METvBin_FBWH_SR, *h_MET_FBWH_SR, *h_mT_FBWH_SR;
  TH1D *h_METvBin_FBW_SR, *h_MET_FBW_SR, *h_mT_FBW_SR;
  TH1D *h_METvBin_FBH_SR, *h_MET_FBH_SR, *h_mT_FBH_SR;
  TH1D *h_METvBin_FBWH_CR, *h_MET_FBWH_CR, *h_mT_FBWH_CR;
  TH1D *h_METvBin_FBW_CR, *h_MET_FBW_CR, *h_mT_FBW_CR;
  TH1D *h_METvBin_FBH_CR, *h_MET_FBH_CR, *h_mT_FBH_CR;

  TH1D *h_WmassJPt_FBWH_SR, *h_WmassJEta_FBWH_SR;
  TH1D *h_WmassJPt_FBW_SR, *h_WmassJEta_FBW_SR;
  TH1D *h_WmassJPt_FBH_SR, *h_WmassJEta_FBH_SR;
  TH1D *h_WmassJPt_FBWH_CR, *h_WmassJEta_FBWH_CR;
  TH1D *h_WmassJPt_FBW_CR, *h_WmassJEta_FBW_CR;
  TH1D *h_WmassJPt_FBH_CR, *h_WmassJEta_FBH_CR;

  TH1D *h_WmassJPt_FBWH_CRwt, *h_HmassJPt_FBWH_CRwt;
  TH1D *h_WmassJPt_FBW_CRwt, *h_HmassJPt_FBW_CRwt;
  TH1D *h_WmassJPt_FBH_CRwt, *h_HmassJPt_FBH_CRwt;
  TH1D *h_METvBin_FBWH_CRwt;
  TH1D *h_METvBin_FBW_CRwt;
  TH1D *h_METvBin_FBH_CRwt;

  TH1D *h_HmassJPt_FBWH_SR, *h_HmassJEta_FBWH_SR;
  TH1D *h_HmassJPt_FBW_SR, *h_HmassJEta_FBW_SR;
  TH1D *h_HmassJPt_FBH_SR, *h_HmassJEta_FBH_SR;
  TH1D *h_HmassJPt_FBWH_CR, *h_HmassJEta_FBWH_CR;
  TH1D *h_HmassJPt_FBW_CR, *h_HmassJEta_FBW_CR;
  TH1D *h_HmassJPt_FBH_CR, *h_HmassJEta_FBH_CR;

  TH1D *h_METvBin_FBWH;
  TH1D *h_METvBin_FBWZ;
  TH1D *h_METvBin_FBW;
  TH1D *h_METvBin_FBH;
  TH1D *h_METvBin_FBZ;
  TH1D *h_METvBin_2T2M;
  TH1D *h_METvBin_1T2M;
  TH1D *h_METvBin_1T1M;

  TH1D *h_EvtTypeWH_0AK8M;
  TH1D *h_MET_catWH[36];
  TH1D *h_HT_catWH[36];
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
  TH1D *h_nAK8_catWH[36];

  TH1D *h_LeadbPairMass_catWH[36];
  TH1D *h_LeadNonbPairMass_catWH[36];

  TH1D *h_mTLead_catWH[36];
  TH1D *h_mT2Lead_catWH[36];

  TH1D *h_deepDoubleBdiscr_catWH[36];
  TH1D *h_deepWdiscr_catWH[36];

  TH1D *h_mEfft_catWH[36];
  TH1D *h_mTbMin_catWH[36];
  TH1D *h_mCT_catWH[36];
  TH1D *h_mTqMin_catWH[36];
  TH1D *h_mCTq_catWH[36];

  //  vector<TString> wzwCats = {"Wmc_Wmc","Wmc_Wmd","Wmd_Wmd","Wmc_Tau21","Wmd_Tau21","Tau21_Tau21","Wmc_AK4"};
  TH1D *h_MET_Wmc_WWZ;
  TH1D *h_METvBin_Wmc_WWZ;
  TH1D *h_SDMass_Wmc_WWZ;
  TH1D *h_MET_Wmd_WWZ;
  TH1D *h_METvBin_Wmd_WWZ;
  TH1D *h_SDMass_Wmd_WWZ;
  TH1D *h_MET_Tau21_WWZ;
  TH1D *h_METvBin_Tau21_WWZ;
  TH1D *h_SDMass_Tau21_WWZ;
  TH1D *h_MET_Tau21AntiTag_WWZ;
  TH1D *h_METvBin_Tau21AntiTag_WWZ;
  TH1D *h_SDMass_Tau21AntiTag_WWZ;
  TH1D *h_Tau21NotDeepWTagged;

  TH1D *h_MET_WZW[17];
  TH1D *h_HT_WZW[17];
  TH1D *h_METvBin_WZW[17];
  TH1D *h_NJets_WZW[17];
  TH1D *h_nAK8_WZW[17];
  TH1D *h_mTBest_WZW[17];
  TH1D *h_LeadAK8Pt_WZW[17];
  TH1D *h_LeadAK8PtvBin_WZW[17];
  TH2D *h2_LeadAK8PtMET_WZW[17];
  TH2D *h2_nJetsMET_WZW[17];
  TH2D *h2_HTMET_WZW[17];
  TH1D *h_AK8Pt_WZW[17];
  TH1D *h_AK8Eta_WZW[17];
  TH1D *h_AK8M_WZW[17];
  TH1D *h_dPhiMETAK8_WZW[17];
  TH1D *h_deepWdiscr_WZW[17];
  TH1D *h_AK4PairPt_WZW[17];
  TH1D *h_AK4PairEta_WZW[17];
  TH1D *h_AK4PairM_WZW[17];
  TH1D *h_dPhiMETAK4Pair_WZW[17];
  TH1D *h_massAK8WAK4nonb_WZW[17];
  TH1D *h_nAK8_0T1M_shapeCorr_WZW;
  TH1D *h_METvBin_0T1M_shapeCorr_WZW;

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
  h_RegionCat = new TH1D("RegionCat","SR and CRs total yield",12,0,12);
  h_EvtTypeFine = new TH1D("EvtTypeFine","Event type fine event category",17,0,17);
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
  h_NVtx = new TH1D("NVtx","no. of vertices",100,0,100);

  h_LeadNonbPairMass = new TH1D("LeadNonbPairMass","M(Lead 2 Non-bjets)",60,0,300);
  h_LeadbPairMass = new TH1D("LeadbPairMass","M(Lead 2 bjets)",60,0,300);

  h_DeepWdiscr = new TH1D("DeepWdiscr","Deep W dicriminator for leading AK8",100,0,1);
  h_DeepZdiscr = new TH1D("DeepZdiscr","Deep Z dicriminator for leading AK8",100,0,1);
  h_DeepWdiscrMD = new TH1D("DeepWdiscrMD","Deep W dicriminator for leading AK8 mass decorrl",100,0,1);
  h_DeepZHdiscrMD = new TH1D("DeepZHdiscrMD","Deep Z/H dicriminator for leading AK8 mass decorrl",100,0,1);
  h_DeepDoubleBdiscr = new TH1D("DeepDoubleBdiscr","deepDoubleBdicsr for leading AK8",100,0,1);

  h_DeepWdisc_BG = new TH1D("DeepWdisc_BG","Deep W disc for AK8 not matched to W->qq",1000,0,1);
  h_DeepWdiscMD_BG = new TH1D("DeepWdiscMD_BG","Deep W disc MD for AK8 not matched to W->qq",1000,0,1);
  h_Tau21_BG = new TH1D("Tau21_BG","Tau21 for AK8 not matched to W->qq",1000,0,1);
  h_AK8Pt_BG = new TH1D("AK8Pt_BG","Pt for AK8 not matched to W->qq",200,0,2000);

  h_DeepWdisc_SigZ = new TH1D("DeepWdisc_SigZ","Deep W disc for AK8 matched to Z->qq, minDR < 0.1",1000,0,1);
  h_DeepWdiscMD_SigZ = new TH1D("DeepWdiscMD_SigZ","Deep W disc MD for AK8 matched to Z->qq, minDR < 0.1",1000,0,1);
  h_Tau21_SigZ = new TH1D("Tau21_SigZ","Tau21 for AK8 matched to Z->qq, minDR < 0.1",1000,0,1);
  h_AK8Pt_SigZ = new TH1D("AK8Pt_SigZ","Pt for AK8 matched to Z->qq",200,0,2000);
  h_AK8M_SigZ = new TH1D("AK8M_SigZ","Mass for AK8 matched to Z->qq",100,0,200);
  h_GenZmassDeepWpass = new TH1D("GenZmassDeepWpass","Mass of Gen Z->qq, minDR < 0.1, AK8 jet passes DeepWdisc",100,0,200);
  h_SDmassDeepWpassGenZmatch = new TH1D("SDmassDeepWpassGenZmatch","SD Mass of AK8 matched to Gen Z->qq, minDR < 0.1, AK8 jet passes DeepWdisc",100,0,200);
  h_Zkids = new TH1D("Zkids","PID of daughters of Z",25,0,25);
  h_GenPartInDeepDoubleBpassAK8 = new TH1D("GenPartInDeepDoubleBpassAK8","PDG ID particles inside (dR < 0.4) AK8 passing DeepDoubleB disc",25,0,25);
  h_GenPartInDeepDoubleBfailAK8 = new TH1D("GenPartInDeepDoubleBfailAK8","PDG ID particles inside (dR < 0.4) AK8 failing DeepDoubleB disc",25,0,25);

  h_dRGenWb = new TH1D("dRGenWb","dR(GenW,b)",40,0,2);
  h_DeepWdisc_SigW = new TH1D("DeepWdisc_SigW","Deep W disc for AK8 matched to W->qq, minDR < 0.1",1000,0,1);
  h_DeepWdiscMD_SigW = new TH1D("DeepWdiscMD_SigW","Deep W disc MD for AK8 matched to W->qq, minDR < 0.1",1000,0,1);
  h_Tau21_SigW = new TH1D("Tau21_SigW","Tau21 for AK8 matched to W->qq, minDR < 0.1",1000,0,1);
  h_AK8Pt_SigW = new TH1D("AK8Pt_SigW","Pt for AK8 matched to W->qq",200,0,2000);
  h_AK8M_SigW = new TH1D("AK8M_SigW","Mass for AK8 matched to W->qq",100,0,200);
  h_GenWmassDeepWpass = new TH1D("GenWmassDeepWpass","Mass of Gen W->qq, minDR < 0.1, AK8 jet passes DeepWdisc",100,0,200);
  h_SDmassDeepWpassGenWmatch = new TH1D("SDmassDeepWpassGenWmatch","SD Mass of AK8 matched to Gen W->qq, minDR < 0.1, AK8 jet passes DeepWdisc",100,0,200);

  h_LeadAK8Mass = new TH1D("LeadAK8Mass","Leading AK8 SD Mass",60,0,300);
  h_AK8MassNearGenH = new TH1D("AK8MassNearGenH","AK8 SD mass near GenH (dR<0.1)",60,0,300);;
  h_AK8MassNearGenW = new TH1D("AK8MassNearGenW","AK8 SD mass near GenW (dR<0.1)",60,0,300);
  h_AK8MassNearGenZ = new TH1D("AK8MassNearGenZ","AK8 SD mass near GenZ (dR<0.1)",60,0,300);

  h_dRbosons = new TH1D("dRbosons","#DeltaR b/w Gen W/Z/H",100,0,5);

  h_TaggedHMassFailM = new TH1D("TaggedHMassFailM","SD mass of AK8 tagged as H, but not within H-mass",200,0,200);
  h_dRTaggedHFailM = new TH1D("dRTaggedHFailM","dR(H tagged AK8 and failing H mass,Gen H)",60,0,3);
  h_TaggedWMassFailM = new TH1D("TaggedWMassFailM","SD mass of AK8 tagged as W, but not within W-mass",200,0,200);
  h_dRTaggedWFailM = new TH1D("dRTaggedWFailM","dR(W tagged AK8 and failing W mass,Gen W)",60,0,3);

  h_dPhiLSPs = new TH1D("dPhiLSPs","dPhi b/w LSPs",40,0,4);

  h_MT = new TH1D("mT","mT(MET,AK8J)",200,0,2000);
  h_MT2J = new TH1D("mT2J","mT(MET,AK8J2)",200,0,2000);  
  h_mTSum = new TH1D("mTSum","mT+mT2J",200,0,4000);
  h_MTvBin = new TH1D("mTvBin","mT(MET,AK8J)",mTvbins.size()-1,&(mTvbins[0]));
  h_MT2JvBin = new TH1D("mT2JvBin","mT(MET,AK8J2)",mT2Jvbins.size()-1,&(mT2Jvbins[0]));
  h_mTSumvBin = new TH1D("mTSumvBin","mT+mT2J",mTSumvbins.size()-1,&(mTSumvbins[0]));

  h_mTbMin = new TH1D("mTbMin","min(mTb1,mTb2)",200,0,2000);
  h_mCT = new TH1D("mCT","mCT(b1,b2)",200,0,2000);

  h_mTRatio = new TH1D("mTRatio","mT/mT2J: ratio of mT1 and mT2J",500,0,5);

  h_dPhiMETAK8 = new TH1D("dPhiMETAK8","dPhi(MET,AK8J)",40,0,4);
  h_dPhiAK8J1J2 = new TH1D("dPhiAK8J1J2","dPhi(AK8J1,AK8J2)",40,0,4);

  h_nAK8 = new TH1D("nAK8","no. of AK8 jets",5,0,5);
  h_AK8J1Pt = new TH1D("AK8Pt","Leading AK8 jets Pt",200,0,2000);
  h_AK8J1Eta = new TH1D("AK8Eta","AK8 Eta",120,-6,6);
  h_AK8J1Mass = new TH1D("AK8Mass","AK8 Mass",60,0,300);
  h_AK8J1MainWdisc = new TH1D("AK8MainWdisc","AK8 MainWdisc",20,0,1);

  h_mTLep = new TH1D("mTLep","mT(lepton, MET)",30,0,300);
  h_LepType = new TH1D("LepType","Lepton: 0->0 lep, 11->Ele, 13->Muon",15,0,15);
  h_LepPt = new TH1D("LepPt","Lepton Pt",100,0,1000);
  h_LepEta = new TH1D("LepEta","Lepton eta",60,-3,3);
  h_LepPhi = new TH1D("LepPhi","Lepton phi",80,-4,4);

  h_J1Pt = new TH1D("J1Pt","leading AK4 Pt",200,0,2000);
  h_J1Eta = new TH1D("J1Eta","leading AK4 Eta",120,-6,6);
  h_J1Phi = new TH1D("J1Phi","leading AK4 Phi",80,-4,4);
  h_J2Pt = new TH1D("J2Pt","2nd leading AK4 Pt",200,0,2000);
  h_J2Eta = new TH1D("J2Eta","2nd leading AK4 Eta",120,-6,6);
  h_J2Phi = new TH1D("J2Phi","J2Phi",80,-4,4);
  
  h_AK8J2Pt = new TH1D("AK8J2Pt","2nd leading AK8 jets Pt",200,0,2000);
  h_AK8J2Eta = new TH1D("AK8J2Eta","AK8J2 Eta",120,-6,6);
  h_AK8J2Mass = new TH1D("AK8J2Mass","AK8J2 Mass",60,0,300);
  h_AK8J2MainWdisc = new TH1D("AK8J2MainWdisc","AK8J2 MainWdisc",20,0,1);

  h_WTagMassJPt = new TH1D("WTagMassJPt","Pt of leading AK8 with W-mass and pass DeepW",100,0,1000);
  h_WAntiTagMassJPt = new TH1D("WAntiTagMassJPt","Pt of leading AK8 with W-mass and fail DeepW",100,0,1000);
  h_HTagMassJPt = new TH1D("HTagMassJPt","Pt of leading AK8 with H-mass and pass Deepbb",100,0,1000);
  h_HAntiTagMassJPt = new TH1D("HAntiTagMassJPt","Pt of leading AK8 with H-mass and fail Deepbb",100,0,1000);

  h_InvMassAK8Jets = new TH1D("InvMassAK8Jets","Invariant mass of leading 2 AK8 jets",200,0,2000);
  h_AK8J1J2MassRatio = new TH1D("AK8J1J2MassRatio","AK8J2 mass/AK8J1 mass",100,0,5);
  h_dRbJetAK8 = new TH1D("dRbJetAK8","#DeltaR(AK8,b-jet)",200,0,2);
  h_dRbJetAK8J2 = new TH1D("dRbJetAK8J2","#DeltaR(AK8J2,b-jet)",200,0,2);

  h_dPhibJetMET = new TH1D("dPhibJetMET","#Delta#Phi(leading b, MET)",40,0,4);
  h_dPhibJetAK8 = new TH1D("dPhibJetAK8","#Delta#Phi(leading b, AK8)",40,0,4);
  h_dPhibJetAK8J2 = new TH1D("dPhibJetAK8J2","#Delta#Phi(leading b, AK8J2)",40,0,4);
  h_deepDoubleBdiscrForHcand = new TH1D("deepDoubleBdiscrForHcand","deepDoubleBdiscriminator value for AK8 jet matched (dR < 0.3) to Gen H",100,-1,1);

  h2_mTbMin_mCT = new TH2D("mTbMin_mCT","x:min(mTb1,mTb2), mCT(b1,b2)",200,0,2000,200,0,2000);
  h2_AK8J1J2Wdisc = new TH2D("AK8J1J2Wdisc","x:AK8J1 mainWdisc vs y:AK8J2 mainWdisc",20,0,1,20,0,1);
  h2_AK8J1J2Mass = new TH2D("AK8J1J2Mass","x:AK8J1 Mass vs y:AK8J2 Mass",60,0,300,60,0,300);
  h2_SusyPDGMass = new TH2D("SusyPDGMass","x:PDG ID of SUSY particle, y:Mass of particle",50,1000000,1000050,300,0,3000);
  h2_LSPPtBosonPt = new TH2D("LSPPtBosonPt","x:LSP Pt, y:Boson Pt",200,0,2000,200,0,2000);

  h2_LeadAK4EtaPhi = new TH2D("LeadAK4EtaPhi","Leading AK4 eta (on x) vs phi (on y)",240,-6,6,160,-4,4);

  h_dPhi1 = new TH1D("DeltaPhi1","DeltaPhi1",40,0,4);
  h_dPhi2 = new TH1D("DeltaPhi2","DeltaPhi2",40,0,4);
  h_dPhi3 = new TH1D("DeltaPhi3","DeltaPhi3",40,0,4);
  h_dPhi4 = new TH1D("DeltaPhi4","DeltaPhi4",40,0,4);

  h_RA2bBins = new TH1D("RA2bBins","RA2b bins",175,0,175);  

  h_RegionCat->Fill("2T2M_SR",0);
  h_RegionCat->Fill("1T2M_SR",0);
  h_RegionCat->Fill("1T1M_SR",0);
  h_RegionCat->Fill("2T2M_CR",0);
  h_RegionCat->Fill("0T2M_CR",0);
  h_RegionCat->Fill("0T1M_CR",0);
  h_RegionCat->Fill("FBWH_SR",0);
  h_RegionCat->Fill("FBW_SR",0);
  h_RegionCat->Fill("FBH_SR",0);
  h_RegionCat->Fill("FBWH_CR",0);
  h_RegionCat->Fill("FBW_CR",0);
  h_RegionCat->Fill("FBH_CR",0);
  
  h_catname2j = (TH1D*)h_RegionCat->Clone("catname2j");

  h_METvBin_2T2M_SR = new TH1D("METvBin_2T2M_SR","MET for 0b SR with 2T2M",METvbins.size()-1,&(METvbins[0]));
  h_METvBin_1T2M_SR = new TH1D("METvBin_1T2M_SR","MET for 0b SR with 1T2M",METvbins.size()-1,&(METvbins[0]));
  h_METvBin_1T1M_SR = new TH1D("METvBin_1T1M_SR","MET for 0b SR with 1T1M",METvbins.size()-1,&(METvbins[0]));
  h_METvBin_2T2M_CR = new TH1D("METvBin_2T2M_CR","MET for 0b CR with 2T2M",METvbins.size()-1,&(METvbins[0]));
  h_METvBin_0T2M_CR = new TH1D("METvBin_0T2M_CR","MET for 0b CR with 0T2M",METvbins.size()-1,&(METvbins[0]));
  h_METvBin_0T1M_CR = new TH1D("METvBin_0T1M_CR","MET for 0b CR with 0T1M",METvbins.size()-1,&(METvbins[0]));

  h_METvBin_FBWH_SR = new TH1D("METvBin_FBWH_SR","MET for 1b SR of FBWH",METvbins.size()-1,&(METvbins[0]));
  h_METvBin_FBW_SR = new TH1D("METvBin_FBW_SR","MET for 1b SR of FBW",METvbins.size()-1,&(METvbins[0]));
  h_METvBin_FBH_SR = new TH1D("METvBin_FBH_SR","MET for 1b SR of FBH",METvbins.size()-1,&(METvbins[0]));
  h_METvBin_FBWH_CR = new TH1D("METvBin_FBWH_CR","MET for 1b CR of FBWH",METvbins.size()-1,&(METvbins[0]));
  h_METvBin_FBW_CR = new TH1D("METvBin_FBW_CR","MET for 1b CR of FBW",METvbins.size()-1,&(METvbins[0]));
  h_METvBin_FBH_CR = new TH1D("METvBin_FBH_CR","MET for 1b CR of FBH ",METvbins.size()-1,&(METvbins[0]));
  h_METvBin_FBWH_CRwt = new TH1D("METvBin_FBWH_CRwt","MET for 1b CRwt of FBWH",METvbins.size()-1,&(METvbins[0]));
  h_METvBin_FBW_CRwt = new TH1D("METvBin_FBW_CRwt","MET for 1b CRwt of FBW",METvbins.size()-1,&(METvbins[0]));
  h_METvBin_FBH_CRwt = new TH1D("METvBin_FBH_CRwt","MET for 1b CRwt of FBH ",METvbins.size()-1,&(METvbins[0]));

  h_MET_2T2M_SR = new TH1D("MET_2T2M_SR","MET for 0b SR with 2T2M",200,0,2000);
  h_MET_1T2M_SR = new TH1D("MET_1T2M_SR","MET for 0b SR with 1T2M",200,0,2000);
  h_MET_1T1M_SR = new TH1D("MET_1T1M_SR","MET for 0b SR with 1T1M",200,0,2000);
  h_MET_2T2M_CR = new TH1D("MET_2T2M_CR","MET for 0b CR with 2T2M",200,0,2000);
  h_MET_0T2M_CR = new TH1D("MET_0T2M_CR","MET for 0b CR with 0T2M",200,0,2000);
  h_MET_0T1M_CR = new TH1D("MET_0T1M_CR","MET for 0b CR with 0T1M",200,0,2000);

  h_MET_FBWH_SR = new TH1D("MET_FBWH_SR","MET for 1b SR of FBWH",200,0,2000);
  h_MET_FBW_SR = new TH1D("MET_FBW_SR","MET for 1b SR of FBW",200,0,2000);
  h_MET_FBH_SR = new TH1D("MET_FBH_SR","MET for 1b SR of FBH",200,0,2000);
  h_MET_FBWH_CR = new TH1D("MET_FBWH_CR","MET for 1b CR of FBWH",200,0,2000);
  h_MET_FBW_CR = new TH1D("MET_FBW_CR","MET for 1b CR of FBW",200,0,2000);
  h_MET_FBH_CR = new TH1D("MET_FBH_CR","MET for 1b CR of FBH ",200,0,2000);

  h_mT_2T2M_SR = new TH1D("mT_2T2M_SR","mT for 0b SR with 2T2M",200,0,2000);
  h_mT_1T2M_SR = new TH1D("mT_1T2M_SR","mT for 0b SR with 1T2M",200,0,2000);
  h_mT_1T1M_SR = new TH1D("mT_1T1M_SR","mT for 0b SR with 1T1M",200,0,2000);
  h_mT_2T2M_CR = new TH1D("mT_2T2M_CR","mT for 0b CR with 2T2M",200,0,2000);
  h_mT_0T2M_CR = new TH1D("mT_0T2M_CR","mT for 0b CR with 0T2M",200,0,2000);
  h_mT_0T1M_CR = new TH1D("mT_0T1M_CR","mT for 0b CR with 0T1M",200,0,2000);

  h_mT_FBWH_SR = new TH1D("mT_FBWH_SR","mT for 1b SR of FBWH",200,0,2000);
  h_mT_FBW_SR = new TH1D("mT_FBW_SR","mT for 1b SR of FBW",200,0,2000);
  h_mT_FBH_SR = new TH1D("mT_FBH_SR","mT for 1b SR of FBH",200,0,2000);
  h_mT_FBWH_CR = new TH1D("mT_FBWH_CR","mT for 1b CR of FBWH",200,0,2000);
  h_mT_FBW_CR = new TH1D("mT_FBW_CR","mT for 1b CR of FBW",200,0,2000);
  h_mT_FBH_CR = new TH1D("mT_FBH_CR","mT for 1b CR of FBH ",200,0,2000);

  h_WmassJPt_FBWH_SR = new TH1D("WmassJPt_FBWH_SR","Pt of leading AK8 with FBWH_SR",100,0,1000);
  h_WmassJPt_FBW_SR = new TH1D("WmassJPt_FBW_SR","Pt of leading AK8 with FBW_SR",100,0,1000);
  h_WmassJPt_FBH_SR = new TH1D("WmassJPt_FBH_SR","Pt of leading AK8 with FBH_SR",100,0,1000);
  h_WmassJPt_FBWH_CR = new TH1D("WmassJPt_FBWH_CR","Pt of leading AK8 with FBWH_CR",100,0,1000);
  h_WmassJPt_FBW_CR = new TH1D("WmassJPt_FBW_CR","Pt of leading AK8 with FBW_CR",100,0,1000);
  h_WmassJPt_FBH_CR = new TH1D("WmassJPt_FBH_CR","Pt of leading AK8 with FBH_CR",100,0,1000);
  h_WmassJPt_FBWH_CRwt = new TH1D("WmassJPt_FBWH_CRwt","Pt of leading AK8 with FBWH_CRwt",100,0,1000);
  h_WmassJPt_FBW_CRwt = new TH1D("WmassJPt_FBW_CRwt","Pt of leading AK8 with FBW_CRwt",100,0,1000);
  h_WmassJPt_FBH_CRwt = new TH1D("WmassJPt_FBH_CRwt","Pt of leading AK8 with FBH_CRwt",100,0,1000);
  h_WmassJEta_FBWH_SR = new TH1D("WmassJEta_FBWH_SR","|Eta| of leading AK8 with FBWH_SR",30,0,3);
  h_WmassJEta_FBW_SR = new TH1D("WmassJEta_FBW_SR","|Eta| of leading AK8 with FBW_SR",30,0,3);
  h_WmassJEta_FBH_SR = new TH1D("WmassJEta_FBH_SR","|Eta| of leading AK8 with FBH_SR",30,0,3);
  h_WmassJEta_FBWH_CR = new TH1D("WmassJEta_FBWH_CR","|Eta| of leading AK8 with FBWH_CR",30,0,3);
  h_WmassJEta_FBW_CR = new TH1D("WmassJEta_FBW_CR","|Eta| of leading AK8 with FBW_CR",30,0,3);
  h_WmassJEta_FBH_CR = new TH1D("WmassJEta_FBH_CR","|Eta| of leading AK8 with FBH_CR",30,0,3);
  h_HmassJPt_FBWH_SR = new TH1D("HmassJPt_FBWH_SR","Pt of leading AK8 with FBWH_SR",100,0,1000);
  h_HmassJPt_FBW_SR = new TH1D("HmassJPt_FBW_SR","Pt of leading AK8 with FBW_SR",100,0,1000);
  h_HmassJPt_FBH_SR = new TH1D("HmassJPt_FBH_SR","Pt of leading AK8 with FBH_SR",100,0,1000);
  h_HmassJPt_FBWH_CR = new TH1D("HmassJPt_FBWH_CR","Pt of leading AK8 with FBWH_CR",100,0,1000);
  h_HmassJPt_FBW_CR = new TH1D("HmassJPt_FBW_CR","Pt of leading AK8 with FBW_CR",100,0,1000);
  h_HmassJPt_FBH_CR = new TH1D("HmassJPt_FBH_CR","Pt of leading AK8 with FBH_CR",100,0,1000);
  h_HmassJPt_FBWH_CRwt = new TH1D("HmassJPt_FBWH_CRwt","Pt of leading AK8 with FBWH_CRwt",100,0,1000);
  h_HmassJPt_FBW_CRwt = new TH1D("HmassJPt_FBW_CRwt","Pt of leading AK8 with FBW_CRwt",100,0,1000);
  h_HmassJPt_FBH_CRwt = new TH1D("HmassJPt_FBH_CRwt","Pt of leading AK8 with FBH_CRwt",100,0,1000);

  h_HmassJEta_FBWH_SR = new TH1D("HmassJEta_FBWH_SR","|Eta| of leading AK8 with FBWH_SR",30,0,3);
  h_HmassJEta_FBW_SR = new TH1D("HmassJEta_FBW_SR","|Eta| of leading AK8 with FBW_SR",30,0,3);
  h_HmassJEta_FBH_SR = new TH1D("HmassJEta_FBH_SR","|Eta| of leading AK8 with FBH_SR",30,0,3);
  h_HmassJEta_FBWH_CR = new TH1D("HmassJEta_FBWH_CR","|Eta| of leading AK8 with FBWH_CR",30,0,3);
  h_HmassJEta_FBW_CR = new TH1D("HmassJEta_FBW_CR","|Eta| of leading AK8 with FBW_CR",30,0,3);
  h_HmassJEta_FBH_CR = new TH1D("HmassJEta_FBH_CR","|Eta| of leading AK8 with FBH_CR",30,0,3);

  h_METvBin_FBWH  = new TH1D("METvBin_FBWH","MET for fully boosted W & H",METvbins.size()-1,&(METvbins[0]));
  h_METvBin_FBWZ  = new TH1D("METvBin_FBWZ","MET for fully boosted W & Z",METvbins.size()-1,&(METvbins[0]));
  h_METvBin_FBW  = new TH1D("METvBin_FBW","MET for fully boosted W only",METvbins.size()-1,&(METvbins[0]));
  h_METvBin_FBH  = new TH1D("METvBin_FBH","MET for fully boosted H only",METvbins.size()-1,&(METvbins[0]));
  h_METvBin_FBZ  = new TH1D("METvBin_FBZ","MET for fully boosted Z only",METvbins.size()-1,&(METvbins[0]));
  h_METvBin_2T2M  = new TH1D("METvBin_2T2M","MET for 2T2M of WZW",METvbins.size()-1,&(METvbins[0]));
  h_METvBin_1T2M  = new TH1D("METvBin_1T2M","MET for 1T2M of WZW",METvbins.size()-1,&(METvbins[0]));
  h_METvBin_1T1M  = new TH1D("METvBin_1T1M","MET for 1T1M of WZW",METvbins.size()-1,&(METvbins[0]));
  
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
	  h_HT_catWH[iHist] = new TH1D("HT_"+catName,"HT for "+catName,300,0,3000);
	  h_METvBin_catWH[iHist] = new TH1D("METvBin_"+catName,"METvBin for "+catName,METvbins.size()-1,&(METvbins[0]));
	  h_mT_catWH[iHist] = new TH1D("mT_"+catName,"mT for "+catName,200,0,2000);
	  h_AK8HPt_catWH[iHist] = new TH1D("AK8HPt_"+catName,"H cand AK8 Pt for "+catName,200,0,2000);
	  h_AK8HEta_catWH[iHist] = new TH1D("AK8HEta_"+catName,"H cand AK8 Eta for "+catName,120,-6,6);
	  h_AK8HM_catWH[iHist] = new TH1D("AK8HM_"+catName,"H cand AK8 Mass for "+catName,40,0,200);
	  h_dPhiMETAK8H_catWH[iHist] = new TH1D("dPhiMETAK8H_"+catName,"DPhi(H cand AK8,MET) for "+catName,80,0,4);
	  h_AK8WPt_catWH[iHist] = new TH1D("AK8WPt_"+catName,"W cand AK8 Pt for "+catName,200,0,2000);
	  h_AK8WEta_catWH[iHist] = new TH1D("AK8WEta_"+catName,"W cand AK8 Eta for "+catName,120,-6,6);
	  h_AK8WM_catWH[iHist] = new TH1D("AK8WM_"+catName,"W cand AK8 Mass for "+catName,40,0,200);
	  h_mT2J_catWH[iHist] = new TH1D("mT2J_"+catName,"mT(2ndAK8,MET) for "+catName,200,0,2000);
	  h_nAK8_catWH[iHist] = new TH1D("nAK8_"+catName,"no. of AK8 jets with Pt > 200, |eta| < 2.0 for WH for "+catName,10,0,10);

	  h_mTLead_catWH[iHist] = new TH1D("mTLead_"+catName,"mT(Leading AK8,MET) for "+catName,200,0,2000);
	  h_mT2Lead_catWH[iHist] = new TH1D("mT2Lead_"+catName,"mT(2nd leading AK8,MET) for "+catName,200,0,2000);

	  h_deepWdiscr_catWH[iHist] = new TH1D("deepWdiscr_"+catName,"Best value of deep W dicriminator for "+catName,100,0,1);
	  h_deepDoubleBdiscr_catWH[iHist] = new TH1D("deepDoubleBdiscr_"+catName,"best value of deepDoubleBdicsr for "+catName,100,0,1);

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
	h_HT_catWH[iHist] = new TH1D("HT_"+catName,"HT for "+catName,300,0,3000);
	h_METvBin_catWH[iHist] = new TH1D("METvBin_"+catName,"METvBin for "+catName,METvbins.size()-1,&(METvbins[0]));
	h_mT_catWH[iHist] = new TH1D("mT_"+catName,"mT for "+catName,200,0,2000);
	h_AK8HPt_catWH[iHist] = new TH1D("AK8HPt_"+catName,"H cand AK8 Pt for "+catName,200,0,2000);
	h_AK8HEta_catWH[iHist] = new TH1D("AK8HEta_"+catName,"H cand AK8 Eta for "+catName,120,-6,6);
	h_AK8HM_catWH[iHist] = new TH1D("AK8HM_"+catName,"H cand AK8 Mass for "+catName,40,0,200);
	h_dPhiMETAK8H_catWH[iHist] = new TH1D("dPhiMETAK8H_"+catName,"DPhi(H cand AK8,MET) for "+catName,80,0,4);
	h_AK8WPt_catWH[iHist] = new TH1D("AK8WPt_"+catName,"W cand AK8 Pt for "+catName,200,0,2000);
	h_AK8WEta_catWH[iHist] = new TH1D("AK8WEta_"+catName,"W cand AK8 Eta for "+catName,120,-6,6);
	h_AK8WM_catWH[iHist] = new TH1D("AK8WM_"+catName,"W cand AK8 Mass for "+catName,40,0,200);
	h_mT2J_catWH[iHist] = new TH1D("mT2J_"+catName,"mT(2ndAK8,MET) for "+catName,200,0,2000);
	h_nAK8_catWH[iHist] = new TH1D("nAK8_"+catName,"no. of AK8 jets with Pt > 200, |eta| < 2.0 for WH for "+catName,10,0,10);

	h_mTLead_catWH[iHist] = new TH1D("mTLead_"+catName,"mT(Leading AK4 cand,MET) for "+catName,200,0,2000);
	h_mT2Lead_catWH[iHist] = new TH1D("mT2Lead_"+catName,"mT(2nd leading AK4 cand,MET) for "+catName,200,0,2000);

	  h_deepWdiscr_catWH[iHist] = new TH1D("deepWdiscr_"+catName,"Best value of deep W dicriminator for "+catName,100,0,1);
	  h_deepDoubleBdiscr_catWH[iHist] = new TH1D("deepDoubleBdiscr_"+catName,"best value of deepDoubleBdicsr for "+catName,100,0,1);
	  
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
	  h_HT_catWH[iHist] = new TH1D("HT_"+catName,"HT for "+catName,300,0,3000);
	  h_METvBin_catWH[iHist] = new TH1D("METvBin_"+catName,"METvBin for "+catName,METvbins.size()-1,&(METvbins[0]));
	  h_mT_catWH[iHist] = new TH1D("mT_"+catName,"mT for "+catName,200,0,2000);
	  h_AK8HPt_catWH[iHist] = new TH1D("AK8HPt_"+catName,"H cand AK8/AK4 pair Pt for "+catName,200,0,2000);
	  h_AK8HEta_catWH[iHist] = new TH1D("AK8HEta_"+catName,"H cand AK8/AK4 pair Eta for "+catName,120,-6,6);
	  h_AK8HM_catWH[iHist] = new TH1D("AK8HM_"+catName,"H cand AK8/AK4 pair Mass for "+catName,40,0,200);
	  h_dPhiMETAK8H_catWH[iHist] = new TH1D("dPhiMETAK8H_"+catName,"DPhi(H cand AK8,MET) for "+catName,80,0,4);
	  h_AK8WPt_catWH[iHist] = new TH1D("AK8WPt_"+catName,"W cand AK8/AK4 pair Pt for "+catName,200,0,2000);
	  h_AK8WEta_catWH[iHist] = new TH1D("AK8WEta_"+catName,"W cand AK8/AK4 pair Eta for "+catName,120,-6,6);
	  h_AK8WM_catWH[iHist] = new TH1D("AK8WM_"+catName,"W cand AK8/AK4 pair Mass for "+catName,40,0,200);
	  h_mT2J_catWH[iHist] = new TH1D("mT2J_"+catName,"mT(2ndAK8,MET) for "+catName,200,0,2000);
	  h_nAK8_catWH[iHist] = new TH1D("nAK8_"+catName,"no. of AK8 jets with Pt > 200, |eta| < 2.0 for WH for "+catName,10,0,10);

	  h_mTLead_catWH[iHist] = new TH1D("mTLead_"+catName,"mT(Leading AK8,MET) for "+catName,200,0,2000);
	  h_mT2Lead_catWH[iHist] = new TH1D("mT2Lead_"+catName,"mT(2nd leading AK8,MET) for "+catName,200,0,2000);
	  
	  h_deepWdiscr_catWH[iHist] = new TH1D("deepWdiscr_"+catName,"Best value of deep W dicriminator for "+catName,100,0,1);
	  h_deepDoubleBdiscr_catWH[iHist] = new TH1D("deepDoubleBdiscr_"+catName,"best value of deepDoubleBdicsr for "+catName,100,0,1);

	  h_LeadNonbPairMass_catWH[iHist] = new TH1D("LeadNonbPairMass_"+catName,"M(Lead 2 Non-bjets) for "+catName,60,0,300);
	  h_LeadbPairMass_catWH[iHist] = new TH1D("LeadbPairMass_"+catName,"M(Lead 2 bjets) for "+catName,60,0,300);
	  
	  h_mEfft_catWH[iHist] = new TH1D("mEfft_"+catName,"HT + MET for "+catName,200,0,2000);
	  h_mTbMin_catWH[iHist] = new TH1D("mTbMin_"+catName,"min(mT_b1, mT_b2) for "+catName,100,0,1000);
	  h_mCT_catWH[iHist] = new TH1D("mCT_"+catName,"mCT(b1,b2) for "+catName,100,0,1000);
	  h_mTqMin_catWH[iHist] = new TH1D("mTqMin_"+catName,"min(mT_q1, mT_q2) from W for "+catName,100,0,1000);
	  h_mCTq_catWH[iHist] = new TH1D("mCTq_"+catName,"mCT(q1, q2) from W for ",100,0,1000);
	  
	  iHist++;
	}
  
  TDirectory *dir_WZW = oFile->mkdir("WZW","WZW hist");
  dir_WZW->cd();
  h_SDMass_Wmc_WWZ = new TH1D("SDMass_Wmc_WWZ","SD mass of 2nd W/Z cand for 2 W/Z tagged events (Wmc_Wmc)",60,0,300);
  h_MET_Wmc_WWZ = new TH1D("MET_Wmc_WWZ","MET for 2 W-mass corrl AK8 tagged events (Wmc_Wmc)",200,0,2000);
  h_METvBin_Wmc_WWZ = new TH1D("METvBin_Wmc_WWZ","MET for 2 W-mass corrl AK8 tagged events (Wmc_Wmc)",METvbins.size()-1,&(METvbins[0]));
  h_SDMass_Wmd_WWZ = new TH1D("SDMass_Wmd_WWZ","SD mass of 2nd W/Z MD cand for 2 W/Z tagged events (Wmc_Wmd)",60,0,300);
  h_MET_Wmd_WWZ = new TH1D("MET_Wmd_WWZ","MET for 1 W-mass corrl and 1W MD AK8 tagged events (Wmc_Wmd)",200,0,2000);
  h_METvBin_Wmd_WWZ = new TH1D("METvBin_Wmd_WWZ","MET for 1 W-mass corrl and 1W MD AK8 tagged events (Wmc_Wmd)",METvbins.size()-1,&(METvbins[0]));
  h_SDMass_Tau21_WWZ = new TH1D("SDMass_Tau21_WWZ","SD mass of 2nd W/Z cand for 2 W/Z tagged events (Wmc_Tau21)",60,0,300);
  h_MET_Tau21_WWZ = new TH1D("MET_Tau21_WWZ","MET for 1 W-mass corrl and 1 Tau21 AK8 tagged events (Wmc_Tau21)",200,0,2000);
  h_METvBin_Tau21_WWZ = new TH1D("METvBin_Tau21_WWZ","MET for 1 W-mass corrl and 1 Tau21 AK8 tagged events (Wmc_Tau21)",METvbins.size()-1,&(METvbins[0]));
  h_SDMass_Tau21AntiTag_WWZ = new TH1D("SDMass_Tau21AntiTag_WWZ","SD mass of 2nd W/Z cand for 2 W/Z tagged events (Wmc_Tau21AntiTag)",60,0,300);
  h_MET_Tau21AntiTag_WWZ = new TH1D("MET_Tau21AntiTag_WWZ","MET for 1 W-mass corrl and 1 Tau21AntiTag AK8 tagged events (Wmc_Tau21AntiTag)",200,0,2000);
  h_METvBin_Tau21AntiTag_WWZ = new TH1D("METvBin_Tau21AntiTag_WWZ","MET for 1 W-mass corrl and 1 Tau21AntiTag AK8 tagged events (Wmc_Tau21AntiTag)",METvbins.size()-1,&(METvbins[0]));
  h_Tau21NotDeepWTagged = new TH1D("Tau21NotDeepWTagged","Tau21 of the AK8 not tagged by DeepW for WZW",100,0,1);

  for(int i=1;i<=h_EvtTypeFine->GetNbinsX();i++){
    catName = h_EvtTypeFine->GetXaxis()->GetBinLabel(i);
    h_MET_WZW[i-1] = new TH1D("MET_"+catName,"MET for WZW for "+catName,200,0,2000);
    h_HT_WZW[i-1] = new TH1D("HT_"+catName,"HT for WZW for "+catName,300,0,3000);
    h_METvBin_WZW[i-1] = new TH1D("METvBin_"+catName,"MET for WZW for "+catName,METvbins.size()-1,&(METvbins[0]));
    h_NJets_WZW[i-1] = new TH1D("NJets_"+catName,"Njets for WZW for "+catName,10,0,10);
    h_nAK8_WZW[i-1] = new TH1D("nAK8_"+catName,"no. of AK8 jets with Pt > 200, |eta| < 2.0 for WZW for "+catName,10,0,10);
    h_mTBest_WZW[i-1] = new TH1D("mTBest_"+catName,"mT(Best AK8/AK4pair, MET) for "+catName,200,0,2000);
    h_LeadAK8Pt_WZW[i-1] = new TH1D("LeadAK8Pt_"+catName,"Leading AK8 Pt for "+catName,200,0,2000);
    h_LeadAK8PtvBin_WZW[i-1] = new TH1D("LeadAK8PtvBin_"+catName,"Leading AK8 Pt, vbin, for "+catName,ptvbins.size()-1,&(ptvbins[0]));
    h2_LeadAK8PtMET_WZW[i-1] = new TH2D("LeadAK8PtMET_"+catName,"x:Leading AK8 Pt vs MET for "+catName,200,0,2000,200,0,2000);
    h2_nJetsMET_WZW[i-1] = new TH2D("nJetsMET_"+catName,"x:n AK4 jets vs MET for "+catName,10,0,10,METvbins.size()-1,&(METvbins[0]));
    h2_HTMET_WZW[i-1] = new TH2D("HTMET_"+catName,"x:HT jets vs MET for "+catName,300,0,3000,METvbins.size()-1,&(METvbins[0]));

    h_AK8Pt_WZW[i-1] = new TH1D("AK8Pt_"+catName,"Best AK8 Pt for "+catName,200,0,2000);
    h_AK8Eta_WZW[i-1] = new TH1D("AK8Eta_"+catName,"Best AK8 eta for "+catName,120,-6,6);
    h_AK8M_WZW[i-1] = new TH1D("AK8M_"+catName,"AK8 pair mass for "+catName,60,0,300);
    h_dPhiMETAK8_WZW[i-1] = new TH1D("dPhiMETAK8_"+catName,"DPhi(MET, Best AK8) for "+catName,40,0,4);
    h_deepWdiscr_WZW[i-1] = new TH1D("deepWdiscr_"+catName,"deepWdiscr of best AK8 for "+catName,100,0,1);

    h_AK4PairPt_WZW[i-1] = new TH1D("AK4PairPt_"+catName,"Best AK4pair Pt for "+catName,200,0,2000);
    h_AK4PairEta_WZW[i-1] = new TH1D("AK4PairEta_"+catName,"Best AK4pair eta for "+catName,120,-6,6);
    h_AK4PairM_WZW[i-1] = new TH1D("AK4PairM_"+catName,"AK4Pair pair mass for "+catName,60,0,300);
    h_dPhiMETAK4Pair_WZW[i-1] = new TH1D("dPhiMETAK4Pair_"+catName,"DPhi(MET, Best AK4Pair"+catName,40,0,4);
    h_massAK8WAK4nonb_WZW[i-1] = new TH1D("massAK8WAK4nonb_"+catName,"mass of AK8 W and nearest (0.8 < dR < 1.2) AK4 "+catName,60,0,300);
  }
  h_nAK8_0T1M_shapeCorr_WZW = new TH1D("nAK8_0T1M_shapeCorr_WZW","no. of AK8 jets with Pt > 200, |eta| < 2.0 for 0T1M with shape corr from 1T1M",10,0,10);
  h_METvBin_0T1M_shapeCorr_WZW = new TH1D("METvBin_0T1M_shapeCorr_WZW","MET for WZW for 0T1M with shape corr from 1T1M",METvbins.size()-1,&(METvbins[0]));
}

SignalReg::SignalReg(const TString &inputFileList, const char *outFileName, const char* dataset) {
  string nameData=dataset;
  TString nameData2 = nameData;
  TChain *tree = new TChain("tree");
  //  tree = new TChain("TreeMaker2/PreSelection");
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

