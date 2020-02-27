#define SignalReg_cxx
#include "SignalReg.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <iostream>
#include <vector>
#include <cstring>
#include <string>
#include <fstream>

using namespace std;

int main(int argc, char* argv[])
{

  if (argc < 2) {
    cerr << "Please give 3 arguments " << "runList " << " " << "outputFileName" << " " << "dataset" << endl;
    return -1;
  }
  const char *inputFileList = argv[1];
  const char *outFileName   = argv[2];
  const char *data          = argv[3];

  SignalReg ana(inputFileList, outFileName, data);
  cout << "dataset " << data << " " << endl;

  ana.EventLoop(data,inputFileList);

  return 0;
}

void SignalReg::EventLoop(const char *data,const char *inputFileList) {
  if (fChain == 0) return;

  Long64_t nentries = fChain->GetEntriesFast();
  cout << "nentries " << nentries << endl;
  cout << "Analyzing dataset " << data << " " << endl;

  TString s_data=data;
  Long64_t nbytes = 0, nb = 0;
  int decade = 0;

  bool isFastSim = false;
  // float xsec = 0.0, numEvents = 0.0;
  if(s_data.Contains("TChi")){
    isFastSim = true;
    //   if(s_data.Contains("TChiWZ_1000")){ xsec = 1.34352e-3; numEvents = 28771;}
    //   else if(s_data.Contains("TChiWZ_800")){ xsec = 4.75843e-3; numEvents = 34036;}
    //   cout<<"Assigning xsec as: "<<xsec<<endl;
    cout<<"Using as FastSim"<<endl;
  }
  Long64_t nEvtSurv = 0;
  int ak8J1Idx = -1;
  float deepCSVvalue = 0.;

  h_cutflow->Fill("0",0);
  h_cutflow->Fill("Weighted",0);    
  h_cutflow->Fill("MET>200",0);
  h_cutflow->Fill("dPhiCuts",0);
  h_cutflow->Fill("photonVeto",0);
  h_cutflow->Fill("LVeto",0);
  h_cutflow->Fill("Filters",0);
  h_cutflow->Fill("HTRatioDphi",0);
  h_cutflow->Fill("PassedTrigger",0);
  h_cutflow->Fill("HEMaffected",0);
  h_cutflow->Fill("NEvtsNoWtLeft",0);

  h_filters->Fill("TotEvnts",0);
  h_filters->Fill("globalSuperTightHalo2016Filter",0);
  h_filters->Fill("HBHENoiseFilter",0);
  h_filters->Fill("HBHEIsoNoiseFilter",0);
  h_filters->Fill("eeBadScFilter",0);
  h_filters->Fill("EcalDeadCellTriggerPrimitiveFilter",0);
  h_filters->Fill("BadChargedCandidateFilter",0);
  h_filters->Fill("BadPFMuonFilter",0);
  h_filters->Fill("NVtx>0",0);
  h_filters->Fill("JetID",0);
  h_filters->Fill("(MET/CaloMET<5.)",0);

  int dataRun = 0;
  if(s_data.Contains("MC_2016")){ dataRun = -2016; lumiInfb = 35.815165; deepCSVvalue = 0.6324;}
  else if(s_data.Contains("MC_2017")){ dataRun = -2017; lumiInfb = 41.486136; deepCSVvalue = 0.4184;}
  else if(s_data.Contains("MC_2018")){ dataRun = -2018; lumiInfb = 59.546381; deepCSVvalue = 0.4184;}
  
  else if(s_data.Contains("2016")){ dataRun = 2016; isMC = false; deepCSVvalue = 0.6324;}
  else if(s_data.Contains("2017")){ dataRun = 2017; isMC = false; deepCSVvalue = 0.4184;}
  else if(s_data.Contains("2018")){ dataRun = 2018; isMC = false; deepCSVvalue = 0.4184;}
  
  lumiInfb = 137.0;
  cout<<"!!!! changing intLumi to 137/fb, although you should have used 2018 intLumi...."<<endl;

  if(dataRun>0) cout<<"Processing it as "<<dataRun<<" data"<<endl;
  else if(dataRun<0) cout<<"Processing it as "<<abs(dataRun)<<" MC"<<endl;
  else cout<<"No specific data/MC year"<<endl;

  for (Long64_t jentry=0; jentry<nentries;jentry++) {
    
    // ==============print number of events done == == == == == == == =
    double progress = 10.0 * jentry / (1.0 * nentries);
    int k = int (progress);
    if (k > decade)
      cout << 10 * k << " %" <<endl;
    decade = k;
    // cout<<"j:"<<jentry<<" fcurrent:"<<fCurrent<<endl;
    // ===============read this entry == == == == == == == == == == == 
    Long64_t ientry = LoadTree(jentry);
    if (ientry < 0) break;
    nb = fChain->GetEntry(jentry);   nbytes += nb;
    //    print(jentry);  
    // if(dataRun==-2016 || dataRun==-2017) wt=Weight*1000.0*lumiInfb*NonPrefiringProb;
    // else if(dataRun <=0) wt=Weight*1000.0*lumiInfb;
    //    else wt = 1.0;
    // if(isFastSim && NumEvents==1 && CrossSection==1){
    //   CrossSection = xsec;
    //   NumEvents = numEvents;
    //   Weight = xsec/numEvents;
    // }
    //    if(jentry>100) break;
    //#################### RA2b cuts
    if(isFastSim) JetID = true;
    if(NJets < 2 || HT < 300 || MHT < 200 || NMuons!=0 || NElectrons!=0 || (MHT/HT > 1.0) || !JetID || !(DeltaPhi1 > 0.5 && DeltaPhi2 > 0.5 && DeltaPhi3 > 0.3 && DeltaPhi4 > 0.3) || isoMuonTracks!=0 || isoElectronTracks!=0 || isoPionTracks!=0) continue;  
    //####################

    wt=Weight*1000.0*lumiInfb;

    h_cutflow->Fill("0",1);
    h_cutflow->Fill("Weighted",wt);
    h_filters->Fill("TotEvnts",1);
    h_filters->Fill("globalSuperTightHalo2016Filter",globalSuperTightHalo2016Filter);
    h_filters->Fill("HBHENoiseFilter",HBHENoiseFilter);
    h_filters->Fill("HBHEIsoNoiseFilter",HBHEIsoNoiseFilter);
    h_filters->Fill("eeBadScFilter",eeBadScFilter);
    h_filters->Fill("EcalDeadCellTriggerPrimitiveFilter",EcalDeadCellTriggerPrimitiveFilter);
    h_filters->Fill("BadChargedCandidateFilter",BadChargedCandidateFilter);
    h_filters->Fill("BadPFMuonFilter",BadPFMuonFilter);
    h_filters->Fill("NVtx>0",NVtx > 0);
    h_filters->Fill("JetID",JetID);
    h_filters->Fill("(MET/CaloMET<5.)",(MET/CaloMET < 5.));

    if(!isFastSim){
      if(!(globalSuperTightHalo2016Filter==1 && HBHENoiseFilter==1 && HBHEIsoNoiseFilter==1 && eeBadScFilter==1 && EcalDeadCellTriggerPrimitiveFilter==1 && BadChargedCandidateFilter && BadPFMuonFilter && NVtx > 0 && JetID && (MET/CaloMET < 5.))) continue;
      h_cutflow->Fill("Filters",wt);
      h_cutflow->Fill("HTRatioDphi",wt);
    }
    bool HEMaffected = false;
    if(dataRun==2018 && RunNum >=319077){
      for(int i=0;i<Jets->size();i++){
	if((*Jets)[i].Pt() < 30) continue;
	if( (*Jets)[i].Eta() >= -3.20 && (*Jets)[i].Eta() <= -1.2 && 
	    (*Jets)[i].Phi() >= -1.77 && (*Jets)[i].Phi() <= -0.67 &&
	    (abs(DeltaPhi(METPhi,(*Jets)[i].Phi())) < 0.5) ){HEMaffected = true; break;}
      }
    }
    if(HEMaffected){
      h_cutflow->Fill("HEMaffected",wt);
      continue;
    }
    //--------------------------triggers
    if(!isMC){
      bool trgPass = false;
      TString trgName;
      for(int i=0;i<TriggerPass->size();i++){
	trgName = (*TriggerNames)[i];
	if(!(trgName.Contains("MET"))) continue;
	if((*TriggerPass)[i]==1 && (trgName.Contains("HLT_PFMET100_PFMHT100_IDTight_v") || trgName.Contains("HLT_PFMET110_PFMHT110_IDTight_v") ||
				    trgName.Contains("HLT_PFMET120_PFMHT120_IDTight_v") || trgName.Contains("HLT_PFMET130_PFMHT130_IDTight_v") ||
				    trgName.Contains("HLT_PFMET140_PFMHT140_IDTight_v") || trgName.Contains("HLT_PFMET90_PFMHT90_IDTight_v") || 
				    trgName.Contains("HLT_PFMETNoMu100_PFMHTNoMu100_IDTight_v") || trgName.Contains("HLT_PFMETNoMu110_PFMHTNoMu110_IDTight_v") ||
				    trgName.Contains("HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_v") || trgName.Contains("HLT_PFMETNoMu130_PFMHTNoMu130_IDTight_v") || 
				    trgName.Contains("HLT_PFMETNoMu140_PFMHTNoMu140_IDTight_v") ||  trgName.Contains("HLT_PFMETNoMu90_PFMHTNoMu90_IDTight_v"))) trgPass = true;
      }
      if(trgPass) h_cutflow->Fill("PassedTrigger",wt);
      else continue;
    }
    //--------------------------end of triggers   
    vector<TLorentzVector> nlsp,lsp,boson;
    vector<int> lspParentIdx, bosonParentIdx;
    for(int i=0;i<GenParticles->size();i++){
      if(abs((*GenParticles_ParentId)[i]) > 1e6 && abs((*GenParticles_PdgId)[i])>=23 && abs((*GenParticles_PdgId)[i])<=25){ boson.push_back((*GenParticles)[i]); bosonParentIdx.push_back((*GenParticles_ParentIdx)[i]);}
      if(abs((*GenParticles_PdgId)[i])>=1e6){
	h2_SusyPDGMass->Fill(abs((*GenParticles_PdgId)[i]),(*GenParticles)[i].M(),wt);
	if(abs((*GenParticles_PdgId)[i])==1000023 || abs((*GenParticles_PdgId)[i])==1000024) nlsp.push_back((*GenParticles)[i]);
	if(abs((*GenParticles_PdgId)[i])==1000022){ lsp.push_back((*GenParticles)[i]); lspParentIdx.push_back((*GenParticles_ParentIdx)[i]);}
      }
    }
    //--------------
    //----MET
    //    if(MET < 200) continue;
    //h_cutflow->Fill("MET>200",wt);
    float dphi1=4, dphi2=4, dphi3=4, dphi4=4;
    if(Jets->size() > 0 && (*Jets)[0].Pt() > 30 && abs((*Jets)[0].Eta()) < 6.0)
      dphi1 = (abs(DeltaPhi(METPhi,(*Jets)[0].Phi())));

    if(Jets->size() > 1 && (*Jets)[1].Pt() > 30 && abs((*Jets)[1].Eta()) < 6.0)
      dphi2 = (abs(DeltaPhi(METPhi,(*Jets)[1].Phi())));

    if(Jets->size() > 2 && (*Jets)[2].Pt() > 30 && abs((*Jets)[2].Eta()) < 6.0)
      dphi3 = (abs(DeltaPhi(METPhi,(*Jets)[2].Phi())));

    if(Jets->size() > 3 && (*Jets)[3].Pt() > 30 && abs((*Jets)[3].Eta()) < 6.0)
      dphi4 = (abs(DeltaPhi(METPhi,(*Jets)[3].Phi())));
    //    print(jentry);
    // cout<<"METPhi:"<<METPhi<<" dphi1:"<<dphi1<<" dphi2:"<<dphi2<<" dphi3:"<<dphi3<<" dphi4:"<<dphi4<<endl;
    
    //    if(!(dphi1 > 0.5 && dphi2 > 0.5 && dphi3 > 0.5 && dphi4 > 0.5)) continue;
    //    h_cutflow->Fill("dPhiCuts",wt);
    //----Photon veto
    int nPhotons=0;
    for(int i=0;i<Photons->size();i++){
      if((*Photons)[i].Pt() > 100 && (*Photons_fullID)[i] && (!(*Photons_hasPixelSeed)[i]) ){ nPhotons++; break;}
    }
    if(nPhotons>0) continue;
    //if(Photons->size()!=0) continue;
    h_cutflow->Fill("photonVeto",wt);

    if(MET < 250) continue;
    if(NJets > 6) continue;
    if(DeltaPhi1 < 1.5) continue;
    //if(BTagsDeepCSV != 0) continue;

    int evtType = getEventType();
    if(evtType==0){
      h_EvtType->Fill("0 Good AK8",wt);
      continue;
    }
    else if(evtType==200) h_EvtType->Fill("2 Boosted",wt);
    else if(evtType==100 || evtType==150) h_EvtType->Fill("1 Boosted",wt);
    else if(evtType==1) h_EvtType->Fill("1 Good AK8",wt);
    else if(evtType==50) h_EvtType->Fill("2 Good AK8",wt);
    //    if(!(evtType==100 || evtType==150)) continue;
    if(evtType!=200) continue;
    //----MT and MT2
    double mt = 0, mt2j = 0, mt2 = 0.;
    mt = sqrt(2*bestAK8J1.Pt()*MET*(1-cos(DeltaPhi(METPhi,bestAK8J1.Phi()))));
    if(bestAK8J2IdxInMainColl != -1) {
      mt2j = sqrt(2*bestAK8J2.Pt()*MET*(1-cos(DeltaPhi(METPhi,bestAK8J2.Phi()))));
      visa = bestAK8J1;
      visb = bestAK8J2;
      met4.SetPtEtaPhiE(MET,0.,METPhi,0.);
      ComputeMT2 mycalc = ComputeMT2(visa,visb,met4,0.,80.);
      mt2 = mycalc.Compute();
      //      cout<<"MET:"<<MET<<" MT2:"<<mt2<<endl;
    }
    else mt2 = 0.;
    if(mt < 500) continue;
    //    cout<<goodAk8.size()<<" :"<<ak8JBoosted.size()<<" nBoosted:"<<goodAk8.size()<<endl;
    //------------ZZMET-----------
    // if(MET < 300 || HT < 400) continue;//ZZMET
    // if(goodAk8.size()<2) continue;//ZZMET
    // if(goodAk8[0].Pt() < 200) continue;//ZZMET
    // if(goodAk8[1].Pt() < 200) continue;//ZZMET
    // if(((*JetsAK8_softDropMass)[0] < 70 || (*JetsAK8_softDropMass)[0] > 100)) continue;//ZZMET
    // if(((*JetsAK8_softDropMass)[1] < 70 || (*JetsAK8_softDropMass)[1] > 100)) continue;//ZZMET
    
    vector<TLorentzVector> bjets;
    for(int i=0;i<Jets_bJetTagDeepCSVprobb->size();i++){
      if((*Jets)[i].Pt() < 30 || abs((*Jets)[i].Eta()) > 2.4) continue;
      if((*Jets_bJetTagDeepCSVprobb)[i]+(*Jets_bJetTagDeepCSVprobbb)[i] > deepCSVvalue)
	bjets.push_back((*Jets)[i]);
    }
    sortTLorVec(&bjets);
    
    vector<TLorentzVector> ak4jNotAK8;
    for(int i=0;i<Jets->size();i++){
      if((*Jets)[i].Pt() > 30. && abs((*Jets)[i].Eta()) < 2.4){
	if((*Jets)[i].DeltaR(bestAK8J1) < 0.8) continue;
	if(bestAK8J2IdxInMainColl != -1 && (*Jets)[i].DeltaR(bestAK8J2) < 0.8) continue;
	ak4jNotAK8.push_back((*Jets)[i]);
      }
    }
    
    TLorentzVector ak4PairBosonLike;
    float minMassDiff = 1000;
    if(bestAK8J2IdxInMainColl==-1){
      if(ak4jNotAK8.size()==1) ak4PairBosonLike = ak4jNotAK8[0];//actually redundant if ak4jNotAK8.size() > 1 cut is applied.
      else if(ak4jNotAK8.size()==2) ak4PairBosonLike = ak4jNotAK8[0] + ak4jNotAK8[1];
      else if(ak4jNotAK8.size()>=3){
	for(int k=0;k<ak4jNotAK8.size();k++){
	  for(int j=k+1;j<ak4jNotAK8.size();j++){
	    if(minMassDiff > (85 - (ak4jNotAK8[k]+ak4jNotAK8[j]).M())){
	      minMassDiff = 85 - (ak4jNotAK8[k]+ak4jNotAK8[j]).M();
	      ak4PairBosonLike = ak4jNotAK8[k]+ak4jNotAK8[j];
	    }
	  }
	}
      }
      if(ak4jNotAK8.size() < 2) continue;
      if(ak4PairBosonLike.M() < 60 || ak4PairBosonLike.M() > 120) continue;
    }
    //    if(goodAk8.size()>=1 && tau21J2 > 0.45) continue;
    // if(bjets.size()>0 && goodAk8.size()>1){
    //   if(goodAk8[1].DeltaR(bjets[0]) < 0.8) continue;//ZZMET
    // }

    h_MET->Fill(MET,wt);
    h_METvBin->Fill(MET,wt);
    h_METvBinZZMET->Fill(MET,wt);
    h_MHT->Fill(MHT,wt);
    h_HT->Fill(HT,wt);
    h_NJets->Fill(NJets,wt);
    h_BTags->Fill(BTagsDeepCSV,wt);
    h_MT2->Fill(mt2,wt);
    h_MT2vBin->Fill(mt2,wt);

    int nAK8 = 0;
    for(int i=0;i<JetsAK8->size();i++){
      if((*JetsAK8)[i].Pt() > 200 && abs((*JetsAK8)[i].Eta()) < 2.4) nAK8++;
    }
    h_nAK8->Fill(nAK8,wt);
    h_dPhi1->Fill(DeltaPhi1,wt);
    h_dPhi2->Fill(DeltaPhi2,wt);
    h_dPhi3->Fill(DeltaPhi3,wt);
    h_dPhi4->Fill(DeltaPhi4,wt);
    h_RA2bBins->Fill(RA2bin,wt);
    if(bestAK8J1IdxInMainColl<0) cout<<"Bug!!!!!";

    if(boson.size()==2) h_dRbosons->Fill(boson[0].DeltaR(boson[1]),wt);
    else h_dRbosons->Fill(4.90,wt);
    h_nAk4jNotAK8->Fill(ak4jNotAK8.size(),wt);
    if(bestAK8J2IdxInMainColl==-1){
      h_dRbestAK8AK4Cand->Fill(bestAK8J1.DeltaR(ak4PairBosonLike),wt);
      h_dPhiBestAK8AK4Cand->Fill(bestAK8J1.DeltaPhi(ak4PairBosonLike),wt);
      h_AK4BosonCandMass->Fill(ak4PairBosonLike.M(),wt);
    }

    if(boson.size()==1) h_dRaK4BosonGenBoson->Fill(ak4PairBosonLike.DeltaR(boson[0]),wt);
    else if(boson.size()==2) h_dRaK4BosonGenBoson->Fill(min((ak4PairBosonLike.DeltaR(boson[0])),(ak4PairBosonLike.DeltaR(boson[1]))),wt);

    if(bestAK8J1IdxInMainColl != -1){
      h_AK8J1Pt->Fill((bestAK8J1.Pt()),wt);
      h_AK8J1Eta->Fill((bestAK8J1.Eta()),wt);
      h_AK8J1Mass->Fill((*JetsAK8_softDropMass)[bestAK8J1IdxInMainColl],wt);
      h_AK8J1Tau21->Fill(tau21[bestAK8J1IdxInMainColl],wt);
      h_MT->Fill(mt,wt);
      h_MTvBin->Fill(mt,wt);
      h_dPhiMETAK8->Fill(abs(DeltaPhi(METPhi,bestAK8J1.Phi())),wt);
    }
    if(bestAK8J2IdxInMainColl != -1){
      h_AK8J2Pt->Fill((bestAK8J2.Pt()),wt);
      h_AK8J2Eta->Fill((bestAK8J2.Eta()),wt);
      h_AK8J2Mass->Fill((*JetsAK8_softDropMass)[bestAK8J2IdxInMainColl],wt);
      h_AK8J2Tau21->Fill(tau21[bestAK8J2IdxInMainColl],wt);
      h_dPhiAK8J1J2->Fill(abs(DeltaPhi(bestAK8J1.Phi(),bestAK8J2.Phi())),wt);
      h_InvMassAK8Jets->Fill((bestAK8J1+bestAK8J2).M(),wt);
      h_AK8J1J2MassRatio->Fill((*JetsAK8_softDropMass)[bestAK8J2IdxInMainColl]/(*JetsAK8_softDropMass)[bestAK8J1IdxInMainColl],wt);
      h_MT2J->Fill(mt2j,wt);
      h_MT2JvBin->Fill(mt2j,wt);
      h_mTSum->Fill(mt+mt2j,wt);
      h_mTSumvBin->Fill(mt+mt2j,wt);
      h_mTRatio->Fill(mt/mt2j,wt);

      h2_AK8J1J2Tau21->Fill(tau21[bestAK8J1IdxInMainColl],tau21[bestAK8J2IdxInMainColl],wt);
      h2_AK8J1J2Mass->Fill((*JetsAK8_softDropMass)[bestAK8J1IdxInMainColl],(*JetsAK8_softDropMass)[bestAK8J2IdxInMainColl],wt);
      for(int i=0;i<lspParentIdx.size();i++){
	for(int j=0;j<bosonParentIdx.size();j++){
	  if(lspParentIdx[i]==bosonParentIdx[j])
	    h2_LSPPtBosonPt->Fill(lsp[i].Pt(),boson[j].Pt(),wt);
	}
      }
    }

    float bPosAK8J1 = 100, bPosAK8orAK4 = 100;
    
    if(bjets.size()==1){
      h_dRbJetAK8->Fill(bjets[0].DeltaR(bestAK8J1),wt);      
      h_dPhibJetMET->Fill(abs(DeltaPhi(METPhi,bjets[0].Phi())),wt);
      h_dPhibJetAK8->Fill(abs(DeltaPhi(METPhi,bestAK8J1.Phi())),wt);
      
      bPosAK8J1 = bjets[0].DeltaR(bestAK8J1);
      if(bestAK8J2IdxInMainColl == -1)
      bPosAK8orAK4 = bjets[0].DeltaR(ak4PairBosonLike);
      else bPosAK8orAK4 = bjets[0].DeltaR(bestAK8J2);
      
      if(bestAK8J2IdxInMainColl != -1){
	h_dPhibJetAK8J2->Fill(abs(DeltaPhi(METPhi,bestAK8J2.Phi())),wt);
	h_dRbJetAK8J2->Fill(bjets[0].DeltaR(bestAK8J2),wt);
      }
      else{
	h_dPhibJetAK4BosonCand->Fill(abs(bjets[0].DeltaPhi(ak4PairBosonLike)),wt);
	h_dRbJetAK4BosonCand->Fill(bjets[0].DeltaR(ak4PairBosonLike),wt);
      }
      h2_bJetPos->Fill(bPosAK8J1,bPosAK8orAK4,wt);
    }
    nEvtSurv++;
    h_cutflow->Fill("NEvtsNoWtLeft",1);
  } // loop over entries
  cout<<"No. of entries survived: "<<nEvtSurv<<endl;
}

int SignalReg::getEventType(){
  int evtType = -1;
  tau21.resize(0);
  i_tau21 = 1.;
  bestAK8J1.SetXYZT(0.,0.,0.,0.);
  bestAK8J2.SetXYZT(0.,0.,0.,0.);
  bestAK8J1IdxInMainColl = -1;
  bestAK8J2IdxInMainColl = -1;
  goodAk8.resize(0);
  goodAk8Mass.resize(0);
  ak8IdxInMainColl.resize(0);
  //-------------
  for(int i=0;i<JetsAK8->size();i++){
    if((*JetsAK8)[i].Pt() < 200 || abs((*JetsAK8)[i].Eta()) > 2.4) continue;
    if(((*JetsAK8_softDropMass)[i] < massLow) || ((*JetsAK8_softDropMass)[i] > massHigh)) continue;
    goodAk8.push_back((*JetsAK8)[i]);
    ak8IdxInMainColl.push_back(i);
    goodAk8Mass.push_back((*JetsAK8_softDropMass)[i]);
    i_tau21 = ((*JetsAK8_NsubjettinessTau2)[i])/((*JetsAK8_NsubjettinessTau1)[i]);
    tau21.push_back(i_tau21);
    if(abs((*JetsAK8)[i].Eta()) < 2.0){
      if(i_tau21 < 0.35){
	if(bestAK8J1IdxInMainColl == -1){
	  bestAK8J1IdxInMainColl = i;
	  bestAK8J1 = (*JetsAK8)[i];
	}
	else if(bestAK8J2IdxInMainColl == -1){
	  bestAK8J2IdxInMainColl = i;
	  bestAK8J2 = (*JetsAK8)[i];
	}
      }
    }//good ak8s
  }
  if(bestAK8J1IdxInMainColl >= 0 && bestAK8J2IdxInMainColl >= 0) return 200;
  //-----------2 boosted jets? ends here
  //-----------0 boosted jets? start
  if(tau21.size()==0) return 0;
  //-----------1 boosted jet? start
  if(bestAK8J1IdxInMainColl >=0){
    //check if there is second AK8
    for(int i=0;i<goodAk8.size();i++){
      if(bestAK8J1IdxInMainColl == ak8IdxInMainColl[i]) continue;
      bestAK8J2 = (*JetsAK8)[ak8IdxInMainColl[i]];
      bestAK8J2IdxInMainColl = ak8IdxInMainColl[i];
    }
    if(bestAK8J2IdxInMainColl >= 0) return 150;//has 1 boosted jet + 1 additional good ak8
    else return 100;//has 1 boosted jet and no additional good ak8
  }
  //--------none of the AK8 jets pass tau21, but there are ak8 jets with Pt,eta and mass. Start here
  else if(tau21.size()!=0){
    for(int i=0;i<goodAk8.size();i++){
      if(bestAK8J1IdxInMainColl == -1){
	bestAK8J1IdxInMainColl = ak8IdxInMainColl[i];
	bestAK8J1 = (*JetsAK8)[i];
      }
      else if(bestAK8J2IdxInMainColl == -1){
	bestAK8J2IdxInMainColl = i;
	bestAK8J2 = (*JetsAK8)[i];
      }
    }
    if(bestAK8J2IdxInMainColl >= 1) return 50;//there are 2 good AK8 jets (both failing Tau21).
    else return 1;//there is only one good ak8
  }
  
  return -1;  
}


void SignalReg::print(Long64_t jentry){
  //cout<<endl;
  TLorentzVector v1,photo;
  cout<<"------------------------------------------------------------"<<endl;
  cout<<"MomMass:"<<SusyMotherMass<<" Kid Mass:"<<SusyLSPMass<<endl;
  for(int i=0;i<GenParticles->size();i++){  
    //    cout<<EvtNum<<" "<<jentry<<" "<<GenParticles->size()<<" "<<i<<" PdgId:"<<(*GenParticles_PdgId)[i]<<" parentId:"<<(*GenParticles_ParentId)[i]<<" parentIndx:"<<(*GenParticles_ParentIdx)[i]<<" Status:"<<(*GenParticles_Status)[i]<<"\tPx :"<<(*GenParticles)[i].Px()<<" Py :"<<(*GenParticles)[i].Py()<<" Pz :"<<(*GenParticles)[i].Pz()<<" E: "<<(*GenParticles)[i].Energy()<<" M:"<<(*GenParticles)[i].M()<<endl;
    cout<<EvtNum<<" "<<jentry<<" "<<GenParticles->size()<<" "<<i<<" PdgId:"<<(*GenParticles_PdgId)[i]<<" parentId:"<<(*GenParticles_ParentId)[i]<<" parentIndx:"<<(*GenParticles_ParentIdx)[i]<<" Status:"<<(*GenParticles_Status)[i]<</*"\tPx:"<<(*GenParticles)[i].Px()<<" Py:"<<(*GenParticles)[i].Py()<<" Pz:"<<(*GenParticles)[i].Pz()<<*/"\tPt:"<<(*GenParticles)[i].Pt()<<" Eta:"<<(*GenParticles)[i].Eta()<<" Phi:"<<(*GenParticles)[i].Phi()<<" M:"<<(*GenParticles)[i].M()<<endl;

  }
  for(int i=0;i<Jets->size();i++){
    cout<<"JetPt:"<<(*Jets)[i].Pt()<<" JetEta:"<<(*Jets)[i].Eta()<<" JetPhi:"<<(*Jets)[i].Phi()<<endl;
  }
  cout<<"MHTPhi:"<<MHTPhi<<" DPhi1:"<<DeltaPhi1<<" DeltaPhi2:"<<DeltaPhi2<<" DeltaPhi3:"<<DeltaPhi3<<" DeltaPhi4:"<<DeltaPhi4<<endl;
  for(int i=0;i<JetsAK8->size();i++){
    cout<<"AK8 pT:"<<(*JetsAK8)[i].Pt()<<" eta:"<<(*JetsAK8)[i].Eta()<<" phi:"<<(*JetsAK8)[i].Phi()<<" softDropM:"<<(*JetsAK8_softDropMass)[i]<<" tau21:"<<(*JetsAK8_NsubjettinessTau2)[i]/(*JetsAK8_NsubjettinessTau1)[i]<<endl;
  }
}

    /*
    if(JetsAK8->size()==0) continue;
    if((*JetsAK8)[0].Pt() < 200 || abs((*JetsAK8)[0].Eta()) > 2.0) continue; //at least 1 AK8

    //    int nBoostedJets = 0;
    for(int i=0;i<JetsAK8->size();i++){
      i_tau21 = ((*JetsAK8_NsubjettinessTau2)[i])/((*JetsAK8_NsubjettinessTau1)[i]);
      tau21.push_back(i_tau21);
      if(catNum==2){//good jets with tau21 > 0.35 only
	if(((*JetsAK8)[i].Pt() > 200) && 
	   (abs((*JetsAK8)[i].Eta()) < 2.0) && 
	   (i_tau21 < 0.35) && 
	   ((*JetsAK8_softDropMass)[i] > 65) && 
	   ((*JetsAK8_softDropMass)[i] < 100)){
	  if(bestAK8J1IdxInMainColl == -1){
	    bestAK8J1IdxInMainColl = i;
	    bestAK8J1 = (*JetsAK8)[i];
	  }
	  else if(bestAK8J2IdxInMainColl == -1){
	    bestAK8J2IdxInMainColl = i;
	    bestAK8J2 = (*JetsAK8)[i];
	  }
	  nGoodAK8TightTau21++;
	  // goodAk8.push_back((*JetsAK8)[i]);
	  // goodAk8Mass.push_back((*JetsAK8_softDropMass)[i]);
	  // ak8IdxInMainColl.push_back(i);
	}
      }
      else if(catNum==1){//good jets with any tau21 value
	if(((*JetsAK8)[i].Pt() > 200) && 
	   (abs((*JetsAK8)[i].Eta()) < 2.0) && 
	   ((*JetsAK8_softDropMass)[i] > 65) && 
	   ((*JetsAK8_softDropMass)[i] < 100)){
	  if(bestAK8J1IdxInMainColl == -1 && i_tau21 < 0.35){
	    bestAK8J1IdxInMainColl = i;
	    bestAK8J1 = (*JetsAK8)[i];
	  }
	  else if(bestAK8J2IdxInMainColl == -1){
	    bestAK8J2IdxInMainColl = i;
	    bestAK8J2 = (*JetsAK8)[i];
	  }	  
	  // goodAk8.push_back((*JetsAK8)[i]);
	  // goodAk8Mass.push_back((*JetsAK8_softDropMass)[i]);
	  // ak8IdxInMainColl.push_back(i);
	  if(i_tau21 < 0.35) nGoodAK8TightTau21++;
	}
      }
    }
    //    print(jentry);

    if(catNum==2){
      if(nGoodAK8TightTau21<2) continue;
    }
    else if(catNum==1){
      if((bestAK8J1IdxInMainColl == -1) || (nGoodAK8TightTau21 > 1)) continue;
    }
    else continue;
    */
