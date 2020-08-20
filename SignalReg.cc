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

  bool isFastSim = false, applySDmassCorr = 1;
  bool isDataSR = true, isSLepCR = 0;

  TFile *fFR = TFile::Open("TFfile_VJets_MC2018.root");
  TH1D *h_WFR = (TH1D*)fFR->Get("WFR");
  TH1D *h_HFR = (TH1D*)fFR->Get("HFR");
  double CRrewtFac = 1;
  bool scaleCR = 1;
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

  h_cutflow->Fill("0",0);
  h_cutflow->Fill("Weighted",0);    
  h_cutflow->Fill("Filters",wt);
  h_cutflow->Fill("HTRatioDphi",wt);
  h_cutflow->Fill("LVeto",0);
  h_cutflow->Fill("PhotonVeto",0);
  h_cutflow->Fill("MET>200",0);
  h_cutflow->Fill("HT>300",0);
  h_cutflow->Fill("NJets2-6",0);
  h_cutflow->Fill("dPhiCuts",0);
  h_cutflow->Fill("nAK8Jets>0",0);
  h_cutflow->Fill("!HEMaffected",0);
  h_cutflow->Fill("PassedTrigger",0);
  h_cutflow->Fill("nAK8>0",0);
  h_cutflow->Fill("mT>500",0);
  h_cutflow->Fill("mTbMin>200",0);
  h_cutflow->Fill("0b_GoodAK8",0);
  h_cutflow->Fill("1b_GoodAK8",0);
  h_cutflow->Fill("NEvtsNoWtLeft",0);
  h_cutflow->Fill("NEvtsLeft",0);

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
  if(s_data.Contains("MC_2016")){ dataRun = -2016; lumiInfb = 35.815165; deepCSVvalue = 0.6321;}
  else if(s_data.Contains("MC_2017")){ dataRun = -2017; lumiInfb = 41.486136; deepCSVvalue = 0.4941;}
  else if(s_data.Contains("MC_2018")){ dataRun = -2018; lumiInfb = 59.546381; deepCSVvalue = 0.4184;}
  
  else if(s_data.Contains("2016")){ dataRun = 2016; isMC = false; deepCSVvalue = 0.6321;}
  else if(s_data.Contains("2017")){ dataRun = 2017; isMC = false; deepCSVvalue = 0.4941;}
  else if(s_data.Contains("2018")){ dataRun = 2018; isMC = false; deepCSVvalue = 0.4184;}

  if(isMC) isDataSR = false;
  if(s_data.Contains("Ele") || s_data.Contains("Muon")){ isDataSR = false; isSLepCR = true;}

  // lumiInfb = 137.0;
  // cout<<"!!!! changing intLumi to 137/fb, although you should have used 2018 intLumi...."<<endl;

  if(dataRun>0) cout<<"Processing it as "<<dataRun<<" data"<<endl;
  else if(dataRun<0) cout<<"Processing it as "<<abs(dataRun)<<" MC"<<endl;
  else cout<<"No specific data/MC year"<<endl;

  if(applySDmassCorr){
    if(abs(dataRun) == 2016 || abs(dataRun) == 2017 || abs(dataRun) == 2018){
      sdCorrFile = TFile::Open("puppiCorr.root");//https://github.com/cms-jet/PuppiSoftdropMassCorrections/tree/80X
      if(sdCorrFile->IsOpen()){
	puppisd_corrGEN      = (TF1*)sdCorrFile->Get("puppiJECcorr_gen");
	puppisd_corrRECO_cen = (TF1*)sdCorrFile->Get("puppiJECcorr_reco_0eta1v3");
	puppisd_corrRECO_for = (TF1*)sdCorrFile->Get("puppiJECcorr_reco_1v3eta2v5");
      }
      else{
	cout<<"!!!!! can't get SD mass corr file. Will not apply any corrn."<<endl;
	applySDmassCorr = false;
      }
    }
  }

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
    if(dataRun==-2016 || dataRun==-2017) wt=Weight*1000.0*lumiInfb*NonPrefiringProb;
    else if(dataRun <=0) wt=Weight*1000.0*lumiInfb;
    else wt = 1.0;
    // if(isFastSim && NumEvents==1 && CrossSection==1){
    //   CrossSection = xsec;
    //   NumEvents = numEvents;
    //   Weight = xsec/numEvents;
    // }
    //    if(jentry>200) break;
    //#################### RA2b cuts
    h_cutflow->Fill("0",1);
    h_cutflow->Fill("Weighted",wt);    

    if(isFastSim) JetID = true;
    if(!isSLepCR){
      if(NJets < 2 || HT < 300 || MHT < 200 || NMuons!=0 || NElectrons!=0 || (MHT/HT > 1.0) || !JetID || !(DeltaPhi1 > 0.5 && DeltaPhi2 > 0.5 && DeltaPhi3 > 0.3 && DeltaPhi4 > 0.3) || isoMuonTracks!=0 || isoElectronTracks!=0 || isoPionTracks!=0) continue;
    }
    else if(s_data.Contains("Ele")){
      if(NJets < 2 || HT < 300 || MHT < 200 || NMuons!=0 || (MHT/HT > 1.0) || !JetID || !(DeltaPhi1 > 0.5 && DeltaPhi2 > 0.5 && DeltaPhi3 > 0.3 && DeltaPhi4 > 0.3) || isoMuonTracks!=0 || isoPionTracks!=0) continue;
    }
    else if(s_data.Contains("Muon")){
      if(NJets < 2 || HT < 300 || MHT < 200 || NElectrons!=0 || (MHT/HT > 1.0) || !JetID || !(DeltaPhi1 > 0.5 && DeltaPhi2 > 0.5 && DeltaPhi3 > 0.3 && DeltaPhi4 > 0.3) || isoElectronTracks!=0 || isoPionTracks!=0) continue;
    }
    //####################

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
    if(!isSLepCR && (NMuons!=0 || NElectrons!=0)) continue; h_cutflow->Fill("LVeto",wt);
    //----Photon veto
    int nPhotons=0;
    for(int i=0;i<Photons->size();i++){
      if((*Photons)[i].Pt() > 100 && (*Photons_fullID)[i] && (!(*Photons_hasPixelSeed)[i]) ){ nPhotons++; break;}
    }
    if(nPhotons>0) continue;
    h_cutflow->Fill("PhotonVeto",wt);

    if(MET < 200) continue; h_cutflow->Fill("MET>200",wt);
    if(HT  < 300) continue; h_cutflow->Fill("HT>300",wt);
    if(NJets > 6) continue; h_cutflow->Fill("NJets2-6",wt);
  
    if(DeltaPhi1 < 1.5) continue; h_cutflow->Fill("dPhiCuts",wt);
    float dphi1=4, dphi2=4, dphi3=4, dphi4=4;
    if(Jets->size() > 0 && (*Jets)[0].Pt() > 30 && abs((*Jets)[0].Eta()) < 6.0)
      dphi1 = (abs(DeltaPhi(METPhi,(*Jets)[0].Phi())));
    if(Jets->size() > 1 && (*Jets)[1].Pt() > 30 && abs((*Jets)[1].Eta()) < 6.0)
      dphi2 = (abs(DeltaPhi(METPhi,(*Jets)[1].Phi())));
    if(Jets->size() > 2 && (*Jets)[2].Pt() > 30 && abs((*Jets)[2].Eta()) < 6.0)
      dphi3 = (abs(DeltaPhi(METPhi,(*Jets)[2].Phi())));
    if(Jets->size() > 3 && (*Jets)[3].Pt() > 30 && abs((*Jets)[3].Eta()) < 6.0)
      dphi4 = (abs(DeltaPhi(METPhi,(*Jets)[3].Phi())));

    if(!(JetsAK8->size()>0 && (*JetsAK8)[0].Pt() > 200 && abs((*JetsAK8)[0].Eta()) < 2.4)) continue;
    h_cutflow->Fill("nAK8Jets>0",wt);

    bool HEMaffected = false;
    if(dataRun==2018 && RunNum >=319077){
      for(int i=0;i<Jets->size();i++){
	if((*Jets)[i].Pt() < 30) continue;
	if( (*Jets)[i].Eta() >= -3.20 && (*Jets)[i].Eta() <= -1.2 && 
	    (*Jets)[i].Phi() >= -1.77 && (*Jets)[i].Phi() <= -0.67 &&
	    (abs(DeltaPhi(METPhi,(*Jets)[i].Phi())) < 0.5) ){HEMaffected = true; break;}
      }
    }
    if(dataRun==-2018){//for MC 2018
      if( (EvtNum % 1000 > 1000*21.0/59.6) && !passHEMjetVeto(30.)) HEMaffected = true;
    }
    if(HEMaffected) continue;
    h_cutflow->Fill("!HEMaffected",wt);
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
    double mTlep = 0;
    int lepType = 0;//11 == ele, 13 == mu
    TLorentzVector lep1;
    if(s_data.Contains("Ele")){
      for(int i=0;i<Electrons->size();i++){
    	if((*Electrons_passIso)[i] == true){lep1 = (*Electrons)[i]; lepType = 11; break;}
      }
      mTlep = sqrt(2*lep1.Pt()*MET*(1-cos(DeltaPhi(METPhi,lep1.Phi()))));
      if(mTlep > 100) continue;
    }
    else if(s_data.Contains("Muon")){

      for(int i=0;i<Muons->size();i++){
      	if((*Muons_passIso)[i] == true){lep1 = (*Muons)[i]; lepType = 13; break;}
      }
      mTlep = sqrt(2*lep1.Pt()*MET*(1-cos(DeltaPhi(METPhi,lep1.Phi()))));
      if(mTlep > 100) continue;
    }

    if(applySDmassCorr) applySDmassCorrAllAK8();//SD mass corrections
    //---------------------Gen infos
    if(isMC){
      vector<TLorentzVector> nlsp,lsp,boson,higgs,genW,genZ;
      vector<int> lspParentIdx, bosonParentIdx;
      for(int i=0;i<GenParticles->size();i++){
	if(abs((*GenParticles_PdgId)[i])>=23 && abs((*GenParticles_PdgId)[i])<=25){ boson.push_back((*GenParticles)[i]); bosonParentIdx.push_back((*GenParticles_ParentIdx)[i]);
	  if(abs((*GenParticles_PdgId)[i])==23) genZ.push_back((*GenParticles)[i]);
	  if(abs((*GenParticles_PdgId)[i])==24) genW.push_back((*GenParticles)[i]);
	  if(abs((*GenParticles_PdgId)[i])==25) higgs.push_back((*GenParticles)[i]);
	}
	if(abs((*GenParticles_PdgId)[i])>=1e6){
	  h2_SusyPDGMass->Fill(abs((*GenParticles_PdgId)[i]),(*GenParticles)[i].M(),wt);
	  if(abs((*GenParticles_PdgId)[i])==1000023 || abs((*GenParticles_PdgId)[i])==1000024) nlsp.push_back((*GenParticles)[i]);
	  if(abs((*GenParticles_PdgId)[i])==1000022){ lsp.push_back((*GenParticles)[i]); lspParentIdx.push_back((*GenParticles_ParentIdx)[i]);}
	}
	if(abs((*GenParticles_PdgId)[i]) < 26 && (*GenParticles_ParentId)[i]==23) h_Zkids->Fill(abs((*GenParticles_PdgId)[i]),wt);
      }
      if(lsp.size()==2) h_dPhiLSPs->Fill(abs(lsp[0].DeltaPhi(lsp[1])),wt);
      else h_dPhiLSPs->Fill(3.5,wt);

      double mindr=1000;
      int mindrIdx = -1;
      for(int i=0;i<JetsAK8->size();i++){
	for(int j=0;j<higgs.size();j++){
	  if(mindr > higgs[j].DeltaR((*JetsAK8)[i])){ mindr = higgs[j].DeltaR((*JetsAK8)[i]); mindrIdx = i;}
	}
      }
      if(mindr<0.3) h_deepDoubleBdiscrForHcand->Fill((*JetsAK8_deepDoubleBDiscriminatorH)[mindrIdx],wt);

      if(genW.size()==0 && genZ.size()==0){
	for(int i=0;i<JetsAK8->size();i++){
	  if((*JetsAK8_softDropMass)[i] < massLow || (*JetsAK8_softDropMass)[i] > massHigh || (*JetsAK8)[i].Pt() < 200 || abs((*JetsAK8)[i].Eta()) > 2.4) continue;
	  h_DeepWdisc_BG->Fill((*JetsAK8_wDiscriminatorDeep)[i],wt);
	  h_DeepWdiscMD_BG->Fill((*JetsAK8_wDiscriminatorDeepDecorrel)[i],wt);
	  h_Tau21_BG->Fill(((*JetsAK8_NsubjettinessTau2)[i])/((*JetsAK8_NsubjettinessTau1)[i]),wt);
	  h_AK8Pt_BG->Fill((*JetsAK8)[i].Pt(),wt);
	}
      }
      else{
	for(int i=0;i<JetsAK8->size();i++){
	  if((*JetsAK8_softDropMass)[i] < massLow || (*JetsAK8_softDropMass)[i] > massHigh || (*JetsAK8)[i].Pt() < 200 || abs((*JetsAK8)[i].Eta()) > 2.4) continue;
	  float minDR_boson = 9999.;
	  int pdgBoson = 0, bosonIdx = -1;
	  for(int j=0;j<GenParticles->size();j++){
	    if((abs((*GenParticles_PdgId)[j]) < 6) && (*GenParticles_ParentIdx)[j] >= 0 &&
	       (minDR_boson > (*JetsAK8)[i].DeltaR((*GenParticles)[(*GenParticles_ParentIdx)[j]])) ){
	      minDR_boson = (*JetsAK8)[i].DeltaR((*GenParticles)[(*GenParticles_ParentIdx)[j]]);
	      pdgBoson = abs((*GenParticles_PdgId)[(*GenParticles_ParentIdx)[j]]);
	      bosonIdx = (*GenParticles_ParentIdx)[j];
	    }
	  }
	  if(minDR_boson < 0.1){
	    if(pdgBoson == 23){
	      h_DeepWdisc_SigZ->Fill((*JetsAK8_wDiscriminatorDeep)[i],wt);
	      h_DeepWdiscMD_SigZ->Fill((*JetsAK8_wDiscriminatorDeepDecorrel)[i],wt);
	      h_Tau21_SigZ->Fill(((*JetsAK8_NsubjettinessTau2)[i])/((*JetsAK8_NsubjettinessTau1)[i]),wt);
	      h_AK8Pt_SigZ->Fill((*JetsAK8)[i].Pt(),wt);
	      h_AK8M_SigZ->Fill((*JetsAK8_softDropMass)[i],wt);
	      if((*JetsAK8_wDiscriminatorDeep)[i] > dwdisValue){
		h_GenZmassDeepWpass->Fill((*GenParticles)[bosonIdx].M(),wt);
		h_SDmassDeepWpassGenZmatch->Fill((*JetsAK8_softDropMass)[i],wt);
	      }
	      if((*JetsAK8_deepDoubleBDiscriminatorH)[i] > deepDoubleBDiscriminatorValue){
		for(int j=0;j<GenParticles->size();j++)
		  if((abs((*GenParticles_PdgId)[j]) < 26) && (*JetsAK8)[i].DeltaR((*GenParticles)[j]) < 0.4) h_GenPartInDeepDoubleBpassAK8->Fill(abs((*GenParticles_PdgId)[j]),wt);
	      }
	      else{
		for(int j=0;j<GenParticles->size();j++)
		  if((abs((*GenParticles_PdgId)[j]) < 26) && (*JetsAK8)[i].DeltaR((*GenParticles)[j]) < 0.4) h_GenPartInDeepDoubleBfailAK8->Fill(abs((*GenParticles_PdgId)[j]),wt);
	      }
	    }
	    else if(pdgBoson == 24){
	      h_DeepWdisc_SigW->Fill((*JetsAK8_wDiscriminatorDeep)[i],wt);
	      h_DeepWdiscMD_SigW->Fill((*JetsAK8_wDiscriminatorDeepDecorrel)[i],wt);
	      h_Tau21_SigW->Fill(((*JetsAK8_NsubjettinessTau2)[i])/((*JetsAK8_NsubjettinessTau1)[i]),wt);
	      h_AK8Pt_SigW->Fill((*JetsAK8)[i].Pt(),wt);
	      h_AK8M_SigW->Fill((*JetsAK8_softDropMass)[i],wt);
	      if((*JetsAK8_wDiscriminatorDeep)[i] > dwdisValue){
		h_GenWmassDeepWpass->Fill((*GenParticles)[bosonIdx].M(),wt);
		h_SDmassDeepWpassGenWmatch->Fill((*JetsAK8_softDropMass)[i],wt);
	      }
	    }
	  }
	  else if(minDR_boson > 1.0){
	    h_DeepWdisc_BG->Fill((*JetsAK8_wDiscriminatorDeep)[i],wt);
	    h_DeepWdiscMD_BG->Fill((*JetsAK8_wDiscriminatorDeepDecorrel)[i],wt);
	    h_Tau21_BG->Fill(((*JetsAK8_NsubjettinessTau2)[i])/((*JetsAK8_NsubjettinessTau1)[i]),wt);
	    h_AK8Pt_BG->Fill((*JetsAK8)[i].Pt(),wt);
	  }
	}
      }
      for(int i=0;i<genW.size();i++){
	for(int j=0;j<GenParticles->size();j++){
	  if(abs((*GenParticles_PdgId)[j]) == 5)
	    h_dRGenWb->Fill(genW[i].DeltaR((*GenParticles)[j]),wt);
	}
      }
      for(int i=0;i<JetsAK8->size();i++){
	for(int j=0;j<genZ.size();j++)
	  if(genZ[j].DeltaR((*JetsAK8)[i]) < 0.1) h_AK8MassNearGenZ->Fill((*JetsAK8_softDropMass)[i],wt);
	for(int j=0;j<genW.size();j++)
	  if(genW[j].DeltaR((*JetsAK8)[i]) < 0.1) h_AK8MassNearGenW->Fill((*JetsAK8_softDropMass)[i],wt);
	for(int j=0;j<higgs.size();j++)
	  if(higgs[j].DeltaR((*JetsAK8)[i]) < 0.1) h_AK8MassNearGenH->Fill((*JetsAK8_softDropMass)[i],wt);      
      }
      if(boson.size()==2) h_dRbosons->Fill(boson[0].DeltaR(boson[1]),wt);
      else h_dRbosons->Fill(4.90,wt);
      for(int i=0;i<lspParentIdx.size();i++){
	for(int j=0;j<bosonParentIdx.size();j++){
	  if(lspParentIdx[i]==bosonParentIdx[j])
	    h2_LSPPtBosonPt->Fill(lsp[i].Pt(),boson[j].Pt(),wt);
	}
      }
    }//isMC
    //------------end of gen info
    int nAK8 = 0;
    for(int i=0;i<JetsAK8->size();i++){
      if((*JetsAK8)[i].Pt() > 200 && abs((*JetsAK8)[i].Eta()) < 2.0) nAK8++;
    }
    if(nAK8==0) continue;
    h_cutflow->Fill("nAK8>0",wt);
    double mt = 0, mt2j = 0;
    if(nAK8 > 0)
      mt = sqrt(2*(*JetsAK8)[0].Pt()*MET*(1-cos(DeltaPhi(METPhi,(*JetsAK8)[0].Phi()))));
    if(nAK8 > 1)
      mt2j = sqrt(2*(*JetsAK8)[1].Pt()*MET*(1-cos(DeltaPhi(METPhi,(*JetsAK8)[1].Phi()))));
    if(mt < 500) continue;
    h_cutflow->Fill("mT>500",wt);
    //---- make a collection of b-jets and non-bjets----------
    nonbjets.resize(0);
    bjets.resize(0);
    for(int i=0;i<Jets_bJetTagDeepCSVprobb->size();i++){
      if((*Jets)[i].Pt() < 30 || abs((*Jets)[i].Eta()) > 2.4) continue;
      if((*Jets_bJetTagDeepCSVprobb)[i]+(*Jets_bJetTagDeepCSVprobbb)[i] > deepCSVvalue)
	bjets.push_back((*Jets)[i]);
      else nonbjets.push_back((*Jets)[i]);
    }
    sortTLorVec(&bjets);
    sortTLorVec(&nonbjets);
    //--------------check if AK8 is near AK4 bjet
    vector<bool> jetsAK8hasb;
    for(int i=0;i<JetsAK8->size();i++){
      bool foundbInAK8 = 0;
      for(int b=0;b<bjets.size();b++){
	if(bjets[b].DeltaR((*JetsAK8)[i]) < 0.8){
	  foundbInAK8 = 1;
	  break;
	}
      }
      jetsAK8hasb.push_back(foundbInAK8);    
    }// AK8 loop for b content

    if(!isMC && !isSLepCR && MET > 200) isDataSR = true;

    TLorentzVector b1,b2,q1,q2;
    double mtbmin = 0, mct = 0, mtqmin = 0, mctq = 0.;
    if(!isDataSR){
      if(bjets.size()==1){
	b1 = bjets[0];
	mtbmin = sqrt(2*b1.Pt()*MET*(1-cos(DeltaPhi(METPhi,b1.Phi()))));
      }
      else if(bjets.size()>=2){ 
	b1 = bjets[0]; 
	b2 = bjets[1];
	mtbmin = min( sqrt(2*b1.Pt()*MET*(1-cos(DeltaPhi(METPhi,b1.Phi())))),
		      sqrt(2*b2.Pt()*MET*(1-cos(DeltaPhi(METPhi,b2.Phi())))));
	mct = sqrt(2*b1.Pt()*b2.Pt()*(1 + cos(b1.DeltaPhi(b2))));
	h2_mTbMin_mCT->Fill(mtbmin,mct,wt);
      }
      //      if(mtbmin < 200) continue;
      if(bjets.size()>0){
	if(mtbmin < 200 && !isSLepCR) continue;
	h_mTbMin->Fill(mtbmin,wt);
	h_mCT->Fill(mct,wt);
	h_cutflow->Fill("mTbMin>200",wt);
      }

      if(bjets.size() >=2) h_LeadbPairMass->Fill((bjets[0]+bjets[1]).M(),wt);
      else h_LeadbPairMass->Fill(0.,wt);
      if(nonbjets.size() >=2) h_LeadNonbPairMass->Fill((nonbjets[0]+nonbjets[1]).M(),wt);
      else h_LeadNonbPairMass->Fill(0.,wt);
    
      if(JetsAK8->size() > 0){
	h_LeadAK8Mass->Fill((*JetsAK8_softDropMass)[0],wt);
	h_DeepWdiscr->Fill((*JetsAK8_wDiscriminatorDeep)[0],wt);
	h_DeepWdiscrMD->Fill((*JetsAK8_wDiscriminatorDeepDecorrel)[0],wt);
	h_DeepZdiscr->Fill((*JetsAK8_zDiscriminatorDeep)[0],wt);
	h_DeepZHdiscrMD->Fill((*JetsAK8_zhDiscriminatorDeepDecorrel)[0],wt);
	h_DeepDoubleBdiscr->Fill((*JetsAK8_deepDoubleBDiscriminatorH)[0],wt);
      }
           
      h_MET->Fill(MET,wt);
      h_METvBin->Fill(MET,wt);
      h_METvBinZZMET->Fill(MET,wt);
      h_MHT->Fill(MHT,wt);
      h_HT->Fill(HT,wt);
      h_NJets->Fill(NJets,wt);
      h_BTags->Fill(BTagsDeepCSV,wt);
      h_NVtx->Fill(NVtx,wt);

      h_mTLep->Fill(mTlep,wt);
      h_LepType->Fill(lepType,wt);
      h_LepPt->Fill(lep1.Pt(),wt);
      h_LepEta->Fill(lep1.Eta(),wt);
      h_LepPhi->Fill(lep1.Phi(),wt);

      h_J1Pt->Fill((*Jets)[0].Pt(),wt);
      h_J1Eta->Fill((*Jets)[0].Eta(),wt);
      h_J1Phi->Fill((*Jets)[0].Phi(),wt);
      h_J2Pt->Fill((*Jets)[1].Pt(),wt);
      h_J2Eta->Fill((*Jets)[1].Eta(),wt);
      h_J2Phi->Fill((*Jets)[1].Phi(),wt);

      h2_LeadAK4EtaPhi->Fill((*Jets)[0].Eta(),(*Jets)[0].Phi(),wt);   
    
      h_nAK8->Fill(nAK8,wt);
      h_dPhi1->Fill(DeltaPhi1,wt);
      h_dPhi2->Fill(DeltaPhi2,wt);
      h_dPhi3->Fill(DeltaPhi3,wt);
      h_dPhi4->Fill(DeltaPhi4,wt);
      h_RA2bBins->Fill(RA2bin,wt);

      h_AK8J1Pt->Fill(((*JetsAK8)[0].Pt()),wt);
      h_AK8J1Eta->Fill(((*JetsAK8)[0].Eta()),wt);
      h_AK8J1Mass->Fill((*JetsAK8_softDropMass)[0],wt);
      h_AK8J1MainWdisc->Fill((*JetsAK8_wDiscriminatorDeep)[0],wt);
      h_MT->Fill(mt,wt);
      h_MTvBin->Fill(mt,wt);
      h_dPhiMETAK8->Fill(abs(DeltaPhi(METPhi,(*JetsAK8)[0].Phi())),wt);

      if(nAK8 > 1){    
	h_AK8J2Pt->Fill(((*JetsAK8)[1].Pt()),wt);
	h_AK8J2Eta->Fill(((*JetsAK8)[1].Eta()),wt);
	h_AK8J2Mass->Fill((*JetsAK8_softDropMass)[1],wt);
	h_AK8J2MainWdisc->Fill((*JetsAK8_wDiscriminatorDeep)[1],wt);
	h_dPhiAK8J1J2->Fill(abs(DeltaPhi((*JetsAK8)[0].Phi(),(*JetsAK8)[1].Phi())),wt);
	h_InvMassAK8Jets->Fill(((*JetsAK8)[0]+(*JetsAK8)[1]).M(),wt);
	h_AK8J1J2MassRatio->Fill((*JetsAK8_softDropMass)[1]/(*JetsAK8_softDropMass)[0],wt);
	h_MT2J->Fill(mt2j,wt);
	h_MT2JvBin->Fill(mt2j,wt);
	h_mTSum->Fill(mt+mt2j,wt);
	h_mTSumvBin->Fill(mt+mt2j,wt);
	h_mTRatio->Fill(mt/mt2j,wt);

	h2_AK8J1J2Wdisc->Fill((*JetsAK8_wDiscriminatorDeep)[0],(*JetsAK8_wDiscriminatorDeep)[1],wt);
	h2_AK8J1J2Mass->Fill((*JetsAK8_softDropMass)[0],(*JetsAK8_softDropMass)[1],wt);
      }
   
      if(bjets.size()==1){
	h_dRbJetAK8->Fill(bjets[0].DeltaR((*JetsAK8)[0]),wt);      
	h_dPhibJetMET->Fill(abs(DeltaPhi(METPhi,bjets[0].Phi())),wt);
	h_dPhibJetAK8->Fill(abs(DeltaPhi(bjets[0].Phi(),(*JetsAK8)[0].Phi())),wt);
      
	if(nAK8 > 1){
	  h_dPhibJetAK8J2->Fill(abs(DeltaPhi(bjets[0].Phi(),(*JetsAK8)[1].Phi())),wt);
	  h_dRbJetAK8J2->Fill(bjets[0].DeltaR(bestAK8J2),wt);
	}
      }
    }//?Data SR

    vector<TLorentzVector> ak8jets;
    vector<double> sdMass;
    vector<int> idxInMainAK8;
    vector<bool> passWmass, passHmass;
    int nDeepW = 0, nDeepbb = 0, nTTau21SR = 0, nMTau21SR = 0, nLTau21CR = 0, nDeepWTau21 = 0, nDeepWTau21LP = 0, nDeepWTau21CR = 0;
    TLorentzVector bestW, bestH;
    int bestWidx = -1, nAK8PtEta = 0, bestHidx = -1, nHmass = 0, nWmass = 0;
    TString catname;

    if(BTagsDeepCSV==0){
      for(int i=0;i<JetsAK8->size();i++){
	if((*JetsAK8)[i].Pt() < 200 || abs((*JetsAK8)[i].Eta()) > 2.0) continue;
	nAK8PtEta++;
	double temp1 = (*JetsAK8_softDropMass)[i];
	if(temp1 < massLow || temp1 > massHigh) continue;
	ak8jets.push_back((*JetsAK8)[i]);
	sdMass.push_back(temp1);
	idxInMainAK8.push_back(i);
      }
      if(nAK8PtEta < 2) continue;
      if(ak8jets.size()==0) continue;
      h_cutflow->Fill("0b_GoodAK8",wt);

      for(int j=0;j<ak8jets.size();j++){
	int i = idxInMainAK8[j];
	float dWv = (*JetsAK8_wDiscriminatorDeep)[i], tau21V = ((*JetsAK8_NsubjettinessTau2)[i])/((*JetsAK8_NsubjettinessTau1)[i]);
	if(dWv > dwdisValue){
	  nDeepW++;
	  if(bestWidx < 0){ bestWidx = i; bestW = ak8jets[j];}
	}

	if(dWv > dwdisValue || tau21V < tau21Value) nDeepWTau21++;
	if(dWv > dwdisValue || tau21V > tau21Value && tau21V < tau21LPValue) nDeepWTau21LP++;
	if(dWv > dwdisValue || tau21V > tau21LPValue) nDeepWTau21CR++;
	
	if(tau21V < tau21Value) nTTau21SR++;
	else if(tau21V < tau21LPValue) nMTau21SR++;
	else nLTau21CR++;
      }
      int nTag = 0, nTagLP = 0;
      if(nDeepW>=1){
	if(nDeepWTau21>=2 && nTTau21SR>=1){ nTag = 2; catname = "2T2M_SR";}
	else if(nDeepWTau21LP>=2 && nMTau21SR>=1){ nTag = 1; catname = "1T2M_SR";}
	else if(nDeepWTau21CR>=2 && nLTau21CR>=1){ nTag = 1; catname = "2T2M_CR";}
	else {nTag = 1; catname = "1T1M_SR";}
      }
      else{
	if(ak8jets.size()>=2) catname = "0T2M_CR";
	else if(ak8jets.size()==1) catname = "0T1M_CR";
	else cout<<"ERROR .... This event with 0 good AK8 must have been rejected before."<<endl;
      }
      if(catname.Contains("_SR") && !isMC && !isSLepCR) isDataSR = true;
      if(isDataSR) continue;

      h_RegionCat->Fill(catname,wt);
      if(catname=="2T2M_SR"){
	h_MET_2T2M_SR->Fill(MET,wt);
	h_METvBin_2T2M_SR->Fill(MET,wt);
	h_mT_2T2M_SR->Fill(MET,wt);
      }
      else if(catname=="1T2M_SR"){
	h_MET_1T2M_SR->Fill(MET,wt);
	h_METvBin_1T2M_SR->Fill(MET,wt);
	h_mT_1T2M_SR->Fill(MET,wt);
      }
      else if(catname=="1T1M_SR"){
	h_MET_1T1M_SR->Fill(MET,wt);
	h_METvBin_1T1M_SR->Fill(MET,wt);
	h_mT_1T1M_SR->Fill(MET,wt);
      }
      else if(catname=="2T2M_CR"){
	h_MET_2T2M_CR->Fill(MET,wt);
	h_METvBin_2T2M_CR->Fill(MET,wt);
	h_mT_2T2M_CR->Fill(MET,wt);
      }
      else if(catname=="0T2M_CR"){
	h_MET_0T2M_CR->Fill(MET,wt);
	h_METvBin_0T2M_CR->Fill(MET,wt);
	h_mT_0T2M_CR->Fill(MET,wt);
      }
      else if(catname=="0T1M_CR"){
	h_MET_0T1M_CR->Fill(MET,wt);
	h_METvBin_0T1M_CR->Fill(MET,wt);
	h_mT_0T1M_CR->Fill(MET,wt);
      }
      else cout<<"!!!Cannot find category 0b"<<endl;      
    }//if 0b
    else{
      int leadPassHmassIdx = -1, leadPassWmassIdx = -1;
      for(int i=0;i<JetsAK8->size();i++){
	if((*JetsAK8)[i].Pt() < 200 || abs((*JetsAK8)[i].Eta()) > 2.0) continue;
	nAK8PtEta++;
	bool foundbInAK8 = jetsAK8hasb[i];
	double temp1 = (*JetsAK8_softDropMass)[i];
	if(!foundbInAK8 && temp1 > massLow && temp1 < massHigh){ nWmass++; passWmass.push_back(1);}
	else passWmass.push_back(0);
	if(foundbInAK8 && temp1 > massLowH && temp1 < massHighH){ nHmass++; passHmass.push_back(1);}
	else passHmass.push_back(0);
	ak8jets.push_back((*JetsAK8)[i]);
	sdMass.push_back(temp1);
	idxInMainAK8.push_back(i);
	if(passWmass[i]==1 && nWmass==1 && (*JetsAK8_wDiscriminatorDeep)[i] > dwdisValue){
	  h_WTagMassJPt->Fill((*JetsAK8)[i].Pt(),wt);
	  leadPassWmassIdx = i;
	}
	else if((passWmass[i]==1 && nWmass==1 && (*JetsAK8_wDiscriminatorDeep)[i] < dwdisValue)){
	  h_WAntiTagMassJPt->Fill((*JetsAK8)[i].Pt(),wt);
	  leadPassWmassIdx = i;
	}
	if(passHmass[i]==1 && nHmass==1 && (*JetsAK8_deepDoubleBDiscriminatorH)[i] > deepDoubleBDiscriminatorValue){
	  h_HTagMassJPt->Fill((*JetsAK8)[i].Pt(),wt);
	  leadPassHmassIdx = i;
	}
	else if(passHmass[i]==1 && nHmass==1 && (*JetsAK8_deepDoubleBDiscriminatorH)[i] < deepDoubleBDiscriminatorValue){
	  h_HAntiTagMassJPt->Fill((*JetsAK8)[i].Pt(),wt);
	  leadPassHmassIdx = i;
	}
      }
      
      if(nAK8PtEta==1){
	int i = idxInMainAK8[0];
	if(jetsAK8hasb[i] && nHmass==1){
	  if((*JetsAK8_deepDoubleBDiscriminatorH)[i] > deepDoubleBDiscriminatorValue){
	    h_AK8Pt_HtHm1AK8->Fill((*JetsAK8)[i].Pt(),wt);
	    h_METvBin_HtHm1AK8->Fill(MET,wt);
	  }
	  else{
	    h_AK8Pt_HatHm1AK8->Fill((*JetsAK8)[i].Pt(),wt);
            h_METvBin_HatHm1AK8->Fill(MET,wt);
	  }
	}
	if(!jetsAK8hasb[i] && nWmass==1){
	  if((*JetsAK8_wDiscriminatorDeep)[i] > dwdisValue){
	    h_AK8Pt_WtWm1AK8->Fill((*JetsAK8)[i].Pt(),wt);
	    h_METvBin_WtWm1AK8->Fill(MET,wt);
	  }
	  else{
	    h_AK8Pt_WatWm1AK8->Fill((*JetsAK8)[i].Pt(),wt);
            h_METvBin_WatWm1AK8->Fill(MET,wt);
	  }
	}	
      }

      if(nAK8PtEta < 2) continue;
      if(ak8jets.size()==0) continue;
      h_cutflow->Fill("1b_GoodAK8",wt);
      int nFBW = 0, nFBH = 0;
      for(int j=0;j<ak8jets.size();j++){
	int i = idxInMainAK8[j];
	float dWv = (*JetsAK8_wDiscriminatorDeep)[i], dbbV = (*JetsAK8_deepDoubleBDiscriminatorH)[i];
	if(dWv > dwdisValue && !(jetsAK8hasb[i])){
	  nDeepW++;
	  if(passWmass[j]){ 
	    nFBW++;
	    if(bestWidx < 0){ bestWidx = i; bestW = ak8jets[j];}
	  }
	}
	if(dbbV > deepDoubleBDiscriminatorValue && jetsAK8hasb[i]){
	  nDeepbb++;
	  if(passHmass[j]){
	    nFBH++;
	    if(bestHidx < 0){ bestHidx = i; bestH = ak8jets[j];}
	  }
	}
      }

      if(nFBW>=1 && nFBH>=1) catname = "FBWH_SR";
      else if(nFBW>=1) catname = "FBW_SR";
      else if(nFBH>=1) catname = "FBH_SR";
      else if(nFBW==0 && nFBH==0 && nHmass>=1 && nWmass>=1) catname = "FBWH_CR";
      else if(nFBH==0 && nHmass>=1) catname = "FBH_CR";
      else if(nFBW==0 && nWmass>=1) catname = "FBW_CR";
      // else if(nDeepbb==0 && nDeepW==0 && nHmass>=1 && nWmass>=1) catname = "FBWH_CR";
      // else if(nDeepbb==0 && nHmass>=1) catname = "FBH_CR";
      // else if(nDeepW==0 && nWmass>=1) catname = "FBW_CR";
      else continue;

      if(catname.Contains("_SR") && !isMC && !isSLepCR) isDataSR = true;
      if(isDataSR) continue;
      h_RegionCat->Fill(catname,wt);
      if(catname.Contains("_CR") && scaleCR){
	CRrewtFac = 1.0;
	if(leadPassWmassIdx >=0){
	  if((*JetsAK8)[leadPassWmassIdx].Pt() < 500) CRrewtFac = h_WFR->GetBinContent(h_WFR->FindBin((*JetsAK8)[leadPassWmassIdx].Pt()));
	  else CRrewtFac = 0.32;
	}
	if(leadPassHmassIdx >=0){
	  if((*JetsAK8)[leadPassHmassIdx].Pt() < 500) CRrewtFac = CRrewtFac * h_HFR->GetBinContent(h_HFR->FindBin((*JetsAK8)[leadPassHmassIdx].Pt()));
	  else CRrewtFac = CRrewtFac*0.64;
	}
      }
      if(catname=="FBWH_SR"){ 
	h_MET_FBWH_SR->Fill(MET,wt);
	h_METvBin_FBWH_SR->Fill(MET,wt);
	h_mT_FBWH_SR->Fill(MET,wt);
	if(bestWidx >=0){
	  h_WmassJPt_FBWH_SR->Fill(bestW.Pt(),wt);
	  h_WmassJEta_FBWH_SR->Fill(abs(bestW.Eta()),wt);	}
	else{h_WmassJPt_FBWH_SR->Fill(0.,wt); h_WmassJEta_FBWH_SR->Fill(2.9,wt);}
	if(bestHidx >=0){
	  h_HmassJPt_FBWH_SR->Fill(bestH.Pt(),wt);
	  h_HmassJEta_FBWH_SR->Fill(abs(bestH.Eta()),wt);	}
	else{h_HmassJPt_FBWH_SR->Fill(0.,wt); h_HmassJEta_FBWH_SR->Fill(2.9,wt);}
      }
      else if(catname=="FBW_SR"){  
	h_MET_FBW_SR->Fill(MET,wt);
	h_METvBin_FBW_SR->Fill(MET,wt);
	h_mT_FBW_SR->Fill(MET,wt);
	if(bestWidx >=0){
	  h_WmassJPt_FBW_SR->Fill(bestW.Pt(),wt);
	  h_WmassJEta_FBW_SR->Fill(abs(bestW.Eta()),wt);	}
	else{h_WmassJPt_FBW_SR->Fill(0.,wt); h_WmassJEta_FBW_SR->Fill(2.9,wt);}
	if(bestHidx >=0){
	  h_HmassJPt_FBW_SR->Fill(bestH.Pt(),wt);
	  h_HmassJEta_FBW_SR->Fill(abs(bestH.Eta()),wt);	}
	else{h_HmassJPt_FBW_SR->Fill(0.,wt); h_HmassJEta_FBW_SR->Fill(2.9,wt);}
      }
      else if(catname=="FBH_SR"){  
	h_MET_FBH_SR->Fill(MET,wt);
	h_METvBin_FBH_SR->Fill(MET,wt);
	h_mT_FBH_SR->Fill(MET,wt);
	if(bestWidx >=0){
	  h_WmassJPt_FBH_SR->Fill(bestW.Pt(),wt);
	  h_WmassJEta_FBH_SR->Fill(abs(bestW.Eta()),wt);	}
	else{h_WmassJPt_FBH_SR->Fill(0.,wt); h_WmassJEta_FBH_SR->Fill(2.9,wt);}
	if(bestHidx >=0){
	  h_HmassJPt_FBH_SR->Fill(bestH.Pt(),wt);
	  h_HmassJEta_FBH_SR->Fill(abs(bestH.Eta()),wt);	}
	else{h_HmassJPt_FBH_SR->Fill(0.,wt); h_HmassJEta_FBH_SR->Fill(2.9,wt);}
      }
      else if(catname=="FBWH_CR"){    
	h_MET_FBWH_CR->Fill(MET,wt);
	h_METvBin_FBWH_CR->Fill(MET,wt);
	h_mT_FBWH_CR->Fill(MET,wt);
	h_METvBin_FBWH_CRwt->Fill(MET,wt*CRrewtFac);
	if(leadPassWmassIdx >=0){
	  h_WmassJPt_FBWH_CR->Fill((*JetsAK8)[leadPassWmassIdx].Pt(),wt);
	  h_WmassJEta_FBWH_CR->Fill(abs((*JetsAK8)[leadPassWmassIdx].Eta()),wt);
	  h_WmassJPt_FBWH_CRwt->Fill((*JetsAK8)[leadPassWmassIdx].Pt(),wt*CRrewtFac); }
	else{h_WmassJPt_FBWH_CR->Fill(0.,wt); h_WmassJEta_FBWH_CR->Fill(2.9,wt);}
	if(leadPassHmassIdx >=0){
	  h_HmassJPt_FBWH_CR->Fill((*JetsAK8)[leadPassHmassIdx].Pt(),wt);
	  h_HmassJEta_FBWH_CR->Fill(abs((*JetsAK8)[leadPassHmassIdx].Eta()),wt);
	  h_HmassJPt_FBWH_CRwt->Fill((*JetsAK8)[leadPassHmassIdx].Pt(),wt*CRrewtFac); }
	else{h_HmassJPt_FBWH_CR->Fill(0.,wt); h_HmassJEta_FBWH_CR->Fill(2.9,wt);}
      }
      else if(catname=="FBW_CR"){     
	h_MET_FBW_CR->Fill(MET,wt);
	h_METvBin_FBW_CR->Fill(MET,wt);
	h_mT_FBW_CR->Fill(MET,wt);
	h_METvBin_FBW_CRwt->Fill(MET,wt*CRrewtFac);
	if(leadPassWmassIdx >=0){
	  h_WmassJPt_FBW_CR->Fill((*JetsAK8)[leadPassWmassIdx].Pt(),wt);
	  h_WmassJEta_FBW_CR->Fill(abs((*JetsAK8)[leadPassWmassIdx].Eta()),wt);
	  h_WmassJPt_FBW_CRwt->Fill(abs((*JetsAK8)[leadPassWmassIdx].Pt()),wt*CRrewtFac); }
	else{h_WmassJPt_FBW_CR->Fill(0.,wt); h_WmassJEta_FBW_CR->Fill(2.9,wt);}
	if(leadPassHmassIdx >=0){
	  h_HmassJPt_FBW_CR->Fill((*JetsAK8)[leadPassHmassIdx].Pt(),wt);
	  h_HmassJEta_FBW_CR->Fill(abs((*JetsAK8)[leadPassHmassIdx].Eta()),wt);
	  h_HmassJPt_FBW_CRwt->Fill(abs((*JetsAK8)[leadPassHmassIdx].Pt()),wt*CRrewtFac); }	
	else{h_HmassJPt_FBW_CR->Fill(0.,wt); h_HmassJEta_FBW_CR->Fill(2.9,wt);}
      }
      else if(catname=="FBH_CR"){
	h_MET_FBH_CR->Fill(MET,wt);
	h_METvBin_FBH_CR->Fill(MET,wt);
	h_mT_FBH_CR->Fill(MET,wt);
	h_METvBin_FBH_CRwt->Fill(MET,wt*CRrewtFac);
	if(leadPassWmassIdx >=0){
	  h_WmassJPt_FBH_CR->Fill((*JetsAK8)[leadPassWmassIdx].Pt(),wt);
	  h_WmassJEta_FBH_CR->Fill(abs((*JetsAK8)[leadPassWmassIdx].Eta()),wt);
	  h_WmassJPt_FBH_CRwt->Fill(abs((*JetsAK8)[leadPassWmassIdx].Pt()),wt*CRrewtFac);	}
	else{h_WmassJPt_FBH_CR->Fill(0.,wt); h_WmassJEta_FBH_CR->Fill(2.9,wt);}
	if(leadPassHmassIdx >=0){
	  h_HmassJPt_FBH_CR->Fill((*JetsAK8)[leadPassHmassIdx].Pt(),wt);
	  h_HmassJEta_FBH_CR->Fill(abs((*JetsAK8)[leadPassHmassIdx].Eta()),wt);
	  h_HmassJPt_FBH_CRwt->Fill(abs((*JetsAK8)[leadPassHmassIdx].Pt()),wt*CRrewtFac);	}
	else{h_HmassJPt_FBH_CR->Fill(0.,wt); h_HmassJEta_FBH_CR->Fill(2.9,wt);}
      }
    }//>=1b

    TString WH_catName;
    int evtType = 0;
    if(BTagsDeepCSV==0){
      evtType = getEventType();
      getEventTypeFine(evtType);//to be called only after calling getEventType();
      getEventTypeWZW();
    }
    else{
      WH_catName = getEventTypeWH();
      if(catname!=catname2j && nAK8PtEta==2 && nHmass==1 && nWmass==1){cout<<catname<<" catname2j:"<<catname2j<<endl;}
    }

    nEvtSurv++;
    h_cutflow->Fill("NEvtsNoWtLeft",1);
    h_cutflow->Fill("NEvtsLeft",wt);
  } // loop over entries
  cout<<"No. of entries survived: "<<nEvtSurv<<endl;
}

int SignalReg::getEventType(){
  if(BTagsDeepCSV > 0) return 0;
  int evtType = -1;
  mainWdisc.resize(0);
  i_mainWdisc = 1.;
  bestAK8J1.SetXYZT(0.,0.,0.,0.);
  bestAK8J2.SetXYZT(0.,0.,0.,0.);
  bestAK8J1IdxInMainColl = -1;
  bestAK8J2IdxInMainColl = -1;
  goodAk8.resize(0);
  goodAk8Mass.resize(0);
  ak8IdxInMainColl.resize(0);
  //-------------
  for(int i=0;i<JetsAK8->size();i++){
    if((*JetsAK8)[i].Pt() < 200 || abs((*JetsAK8)[i].Eta()) > 2.0) continue;
    if(((*JetsAK8_softDropMass)[i] < massLow) || ((*JetsAK8_softDropMass)[i] > massHigh)) continue;
    goodAk8.push_back((*JetsAK8)[i]);
    ak8IdxInMainColl.push_back(i);
    goodAk8Mass.push_back((*JetsAK8_softDropMass)[i]);
    //    i_mainWdisc = ((*JetsAK8_NsubjettinessTau2)[i])/((*JetsAK8_NsubjettinessTau1)[i]);
    i_mainWdisc = (*JetsAK8_wDiscriminatorDeep)[i];
    mainWdisc.push_back(i_mainWdisc);
    if(abs((*JetsAK8)[i].Eta()) < 2.0){
      if(i_mainWdisc > dwdisValue){
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
  //-------------tagging variables
  vector<int> deepWIdx, tau21Idx;
  double sdMass = 0, maxDeepW = 0., minTau21 = 100.;
  int bestTau21ak8Idx = -2, bestDeepWIdx = -1;
  for(int i=0;i<JetsAK8->size();i++){
    if((*JetsAK8)[i].Pt() < 200 || abs((*JetsAK8)[i].Eta()) > 2.0) continue;
    sdMass = (*JetsAK8_softDropMass)[i];
    if(sdMass > massLow && sdMass < massHigh){
      if((*JetsAK8_wDiscriminatorDeep)[i] > dwdisValue){
	deepWIdx.push_back(i);
	if(maxDeepW < (*JetsAK8_wDiscriminatorDeep)[i]){
	  maxDeepW = (*JetsAK8_wDiscriminatorDeep)[i];
	  bestDeepWIdx = i;
	}
      }
      if(((*JetsAK8_NsubjettinessTau2)[i])/((*JetsAK8_NsubjettinessTau1)[i]) < tau21Value){
	tau21Idx.push_back(i);
	if(minTau21 > ((*JetsAK8_NsubjettinessTau2)[i])/((*JetsAK8_NsubjettinessTau1)[i])){
	  minTau21 = ((*JetsAK8_NsubjettinessTau2)[i])/((*JetsAK8_NsubjettinessTau1)[i]);
	  bestTau21ak8Idx = i;
	}	  
      }
    }
  }
  if(bestDeepWIdx == bestTau21ak8Idx){
    if(deepWIdx.size()==1 && tau21Idx.size()==1) bestTau21ak8Idx = -2;//forget about tau21, potential 1TxM
    if(deepWIdx.size()>=1 && tau21Idx.size()>=2){ bestTau21ak8Idx = (bestTau21ak8Idx==tau21Idx[0]) ? tau21Idx[1] : tau21Idx[0];}
    if(deepWIdx.size()>=2 && tau21Idx.size()==1){ bestDeepWIdx = (bestDeepWIdx==deepWIdx[0]) ? deepWIdx[1] : deepWIdx[0];}
  }

  bestAK8J1IdxInMainColl = bestDeepWIdx;
  if(bestDeepWIdx >=0) bestAK8J1 = (*JetsAK8)[bestDeepWIdx];
  //-----------------------------------------------
  if(bestAK8J1IdxInMainColl >= 0 && bestAK8J2IdxInMainColl >= 0) evtType =  200;
  //-----------2 boosted jets? ends here
  //-----------0 boosted jets? start
  if(mainWdisc.size()==0) evtType =  0;
  //-----------1 boosted jet? start
  if(bestAK8J1IdxInMainColl >=0 && bestAK8J2IdxInMainColl < 0){
    //check if there is second AK8
    for(int i=0;i<goodAk8.size();i++){
      if(bestAK8J1IdxInMainColl == ak8IdxInMainColl[i]) continue;
      bestAK8J2 = (*JetsAK8)[ak8IdxInMainColl[i]];
      bestAK8J2IdxInMainColl = ak8IdxInMainColl[i];
      break;
    }
    if(bestAK8J2IdxInMainColl >= 0) evtType =  150;//has 1 boosted jet + 1 additional good ak8. 1T2M
    else evtType =  100;//has 1 boosted jet and no additional good ak8. 1T1M
  }
  //--------none of the AK8 jets pass mainWdisc, but there are ak8 jets with Pt,eta and mass. Start here
  // else if(mainWdisc.size()!=0){
  //   for(int i=0;i<goodAk8.size();i++){
  //     if(bestAK8J1IdxInMainColl == -1){
  // 	bestAK8J1IdxInMainColl = ak8IdxInMainColl[i];
  // 	bestAK8J1 = (*JetsAK8)[i];
  //     }
  //     else if(bestAK8J2IdxInMainColl == -1){
  // 	bestAK8J2IdxInMainColl = i;
  // 	bestAK8J2 = (*JetsAK8)[i];
  //     }
  //   }
  //   if(bestAK8J2IdxInMainColl >= 1) evtType =  50;//there are 2 good AK8 jets (both failing MainWdisc). 0T2M
  //   else evtType =  1;//there is only one good ak8. 0T1M
  // }
  return evtType;  
}

void SignalReg::getEventTypeWZW(){
  if(BTagsDeepCSV > 0) return;
  vector<TString> name;
  int bestWidx = -1, nWmd = 0, nZmd = 0, nWmc = 0, nZmc = 0, nTau21 = 0, nTau21AntiTag = 0;
  double mt = 0., sdMass_Wmd = 0., sdMass_Zmd = 0., sdMass_Wmc = 0., sdMass_Zmc = 0., sdMass_tau21 = 0., sdMass_tau21AntiTag = 0., sdMass = 0.;

  if(JetsAK8->size() > 0 && (*JetsAK8)[0].Pt() > 200 && abs((*JetsAK8)[0].Eta()) < 2.0) 
    mt = sqrt(2*(*JetsAK8)[0].Pt()*MET*(1-cos(DeltaPhi(METPhi,(*JetsAK8)[0].Phi()))));
  if(mt < 500) return;
  if(bestAK8J1IdxInMainColl < 0) return;
  for(int i=0; i < JetsAK8->size(); i++){
    if((*JetsAK8)[i].Pt() < 200 || abs((*JetsAK8)[i].Eta()) > 2.0 || i==bestAK8J1IdxInMainColl) continue;
    sdMass = (*JetsAK8_softDropMass)[i];
    if(sdMass < massLow || sdMass > massHigh) continue;
    if((*JetsAK8_wDiscriminatorDeep)[i] > dwdisValue){
      nWmc++;
      if(nWmc==1) sdMass_Wmc = sdMass;
    }
    if((*JetsAK8_zDiscriminatorDeep)[i] > dzdisValue){ 
      nZmc++;
      if(nZmc==1) sdMass_Zmc = sdMass;
    }
    if((*JetsAK8_wDiscriminatorDeepDecorrel)[i] > dwdisMDvalue){
      nWmd++;
      if(nWmd==1) sdMass_Wmd = sdMass;
    }
    if((*JetsAK8_zhDiscriminatorDeepDecorrel)[i] > dzhdisMDvalue){ 
      nZmd++;
      if(nZmd==1) sdMass_Zmd = sdMass;
    }
    if(((*JetsAK8_NsubjettinessTau2)[i])/((*JetsAK8_NsubjettinessTau1)[i]) < tau21Value){
      nTau21++;
      if(nTau21==1) sdMass_tau21 = sdMass;
    }
    if(((*JetsAK8_NsubjettinessTau2)[i])/((*JetsAK8_NsubjettinessTau1)[i]) > tau21CRValue){
      nTau21AntiTag++;
      if(nTau21AntiTag==1) sdMass_tau21AntiTag = sdMass;
    }
    h_Tau21NotDeepWTagged->Fill(((*JetsAK8_NsubjettinessTau2)[i])/((*JetsAK8_NsubjettinessTau1)[i]),wt);
  }
  // if(EvtNum==42835){
  //   print(0);
  //   cout<<"bestAK8J1IdxInMainColl:"<<bestAK8J1IdxInMainColl<<" nWmc:"<<nWmc<<" nWmd:"<<nWmd<<" tau21:"<<nTau21<<" sdMass_Wmc:"<<sdMass_Wmc<<" sdMass_Wmd:"<<sdMass_Wmd<<" sdMass_tau21:"<<sdMass_tau21<<endl;
  // }
  if(nWmc>=1){
    h_SDMass_Wmc_WWZ->Fill(sdMass_Wmc,wt);
    h_MET_Wmc_WWZ->Fill(MET,wt);
    h_METvBin_Wmc_WWZ->Fill(MET,wt);
  }
  if(nWmd>=1){
    h_SDMass_Wmd_WWZ->Fill(sdMass_Wmd,wt);
    h_MET_Wmd_WWZ->Fill(MET,wt);
    h_METvBin_Wmd_WWZ->Fill(MET,wt);
  }
  if(nTau21>=1){
    h_SDMass_Tau21_WWZ->Fill(sdMass_tau21,wt);
    h_MET_Tau21_WWZ->Fill(MET,wt);
    h_METvBin_Tau21_WWZ->Fill(MET,wt);
  }
  else if(nTau21AntiTag>=1){
    h_SDMass_Tau21AntiTag_WWZ->Fill(sdMass_tau21AntiTag,wt);
    h_MET_Tau21AntiTag_WWZ->Fill(MET,wt);
    h_METvBin_Tau21AntiTag_WWZ->Fill(MET,wt);
  }
}

void SignalReg::getEventTypeFine(int evtType){
  int nTag = 0, nMass = 0, nmass = 0, nTagTau21 = 0, nTagTau21LP = 0, nDeepWTag = 0, nAK8 = 0;
  bool has2AK8 = false; bool printEvt = 0;
  if(JetsAK8->size()==0){
    nTag = 0;
    nMass = 0;
    nTagTau21 = 0;
  }
  else if(JetsAK8->size()>0){
    for(int i=0;i<JetsAK8->size();i++){
      if((*JetsAK8)[i].Pt() > 200 && abs((*JetsAK8)[i].Eta()) < 2.0){
	nAK8++;
	if(((*JetsAK8_softDropMass)[i] > massLow) && ((*JetsAK8_softDropMass)[i] < massHigh)){
	  nMass++;
	  //	if( (((*JetsAK8_NsubjettinessTau2)[i])/((*JetsAK8_NsubjettinessTau1)[i])) < mainWdiscValue)
	  if( ((*JetsAK8_wDiscriminatorDeep)[i] > dwdisValue))
	    nTag++;
	  if(bestAK8J1IdxInMainColl >=0 && bestAK8J1IdxInMainColl != i && (((*JetsAK8_NsubjettinessTau2)[i])/((*JetsAK8_NsubjettinessTau1)[i]) < tau21Value))
	    nTagTau21++;
	  if(bestAK8J1IdxInMainColl >=0 && bestAK8J1IdxInMainColl != i && (((*JetsAK8_NsubjettinessTau2)[i])/((*JetsAK8_NsubjettinessTau1)[i]) > tau21Value) && (((*JetsAK8_NsubjettinessTau2)[i])/((*JetsAK8_NsubjettinessTau1)[i]) < tau21LPValue))
	    nTagTau21LP++;
	  if( ((*JetsAK8_wDiscriminatorDeep)[i] > dwdisValue)) nDeepWTag++;//Should this be inside mass cut? 0T1M vvvv
	}
	//	if( ((*JetsAK8_wDiscriminatorDeep)[i] > dwdisValue)) nDeepWTag++;//Should this be inside mass cut? 0T1M vvv
	//	cout<<((*JetsAK8_wDiscriminatorDeep)[i])<<"----";
	if(i>0) has2AK8 = true;
      }
    }
    if(nTag>=1)
      nTag = (nTagTau21 == 0) ? 1 : 2;
    if(nMass>=3) nMass = 2;
  }

  //--------- AK4 ---------
  vector<TLorentzVector> ak4jNotAK8;
  TLorentzVector topCand; topCand.SetPtEtaPhiE(0.,0.,0.,0.);
  for(int i=0;i<Jets->size();i++){
    if((*Jets)[i].Pt() > 30. && abs((*Jets)[i].Eta()) < 2.4){
      if(bestAK8J1IdxInMainColl>=0){
	if(topCand.Pt() < 10 && bestAK8J1.DeltaR((*Jets)[i]) > 0.8 && bestAK8J1.DeltaR((*Jets)[i]) < 1.2)
	  topCand = bestAK8J1 + (*Jets)[i];
      }
      else if(JetsAK8->size() > 0 && (*JetsAK8)[0].Pt() > 200 && abs((*JetsAK8)[0].Eta()) < 2.0){
	if(topCand.Pt() < 10 && (*JetsAK8)[0].DeltaR((*Jets)[i]) > 0.8 && (*JetsAK8)[0].DeltaR((*Jets)[i]) < 1.2)
	  topCand = (*JetsAK8)[0] + (*Jets)[i];
      }
      if((*Jets)[i].DeltaR(bestAK8J1) < 0.8) continue;
      if(bestAK8J2IdxInMainColl != -1 && (*Jets)[i].DeltaR(bestAK8J2) < 0.8) continue;
      ak4jNotAK8.push_back((*Jets)[i]);
    }
  }
    
  TLorentzVector ak4PairBosonLike;
  float minMassDiff = 1000;
  if(bestAK8J2IdxInMainColl==-1){
    //    if(ak4jNotAK8.size()==1) ak4PairBosonLike = ak4jNotAK8[0];//actually redundant if ak4jNotAK8.size() > 1 cut is applied.
    if(ak4jNotAK8.size()>=2) ak4PairBosonLike = ak4jNotAK8[0] + ak4jNotAK8[1];
    // else if(ak4jNotAK8.size()>=3){
    //   for(int k=0;k<ak4jNotAK8.size();k++){
    // 	for(int j=k+1;j<ak4jNotAK8.size();j++){
    // 	  if(minMassDiff > (80 - (ak4jNotAK8[k]+ak4jNotAK8[j]).M())){
    // 	    minMassDiff = 80 - (ak4jNotAK8[k]+ak4jNotAK8[j]).M();
    // 	    ak4PairBosonLike = ak4jNotAK8[k]+ak4jNotAK8[j];
    // 	  }
    // 	}
    //   }
    // }
    if(ak4jNotAK8.size() >= 2 && ak4PairBosonLike.M() > 65 && ak4PairBosonLike.M() < 105) nmass = 1;
    else nmass = 0;
  }
  double mt = 0, mtak4 = 0.;
  TString catName;
  int iHist = 0;
  bool accEvt = 1;

  if(has2AK8){
    if(nTag==1 && nMass==2 && nTagTau21LP==0) accEvt = 0;
    if(nTag==0 && nDeepWTag!=0) accEvt = 0;//for CRs, antitagged
    catName = to_string(nTag)+"T"+to_string(nMass)+"M";
    if(accEvt) h_EvtTypeFine->Fill(catName,wt);
    if(JetsAK8->size() > 0 && (*JetsAK8)[0].Pt() > 200. && abs((*JetsAK8)[0].Eta()) < 2.0) 
      mt = sqrt(2*(*JetsAK8)[0].Pt()*MET*(1-cos(DeltaPhi(METPhi,(*JetsAK8)[0].Phi()))));
    //      mt = sqrt(2*bestAK8J1.Pt()*MET*(1-cos(DeltaPhi(METPhi,bestAK8J1.Phi()))));
    //      if(nTag==1 && nMass==2){print(0);cout<<"nTagTau21LP:"<<nTagTau21LP<<" accEvt:"<<accEvt<<endl;}
    //    if(!accEvt) {cout<<"!!! oooo nMass:"<<nMass<<" nTag:"<<nTag<<" has2ndAK8:"<<has2AK8<<" accEvt:"<<accEvt<<" catName:"<<catName<<" nDeepWTag:"<<nDeepWTag<<" bestAK8J1IdxInMainColl:"<<bestAK8J1IdxInMainColl<<endl; print(0);}
    if(mt > 500 && accEvt){
      iHist = h_EvtTypeFine->GetXaxis()->FindBin(catName) - 1;
      if(iHist >= 17) cout<<catName<<":Overflow in EvtTypeFine WZW."<<endl;

      h_MET_WZW[iHist]->Fill(MET,wt);
      h_HT_WZW[iHist]->Fill(HT,wt);
      h_METvBin_WZW[iHist]->Fill(MET,wt);
      h_NJets_WZW[iHist]->Fill(NJets,wt);
      h_mTBest_WZW[iHist]->Fill(mt,wt);
      h_LeadAK8Pt_WZW[iHist]->Fill((*JetsAK8)[0].Pt(),wt);
      h_LeadAK8PtvBin_WZW[iHist]->Fill((*JetsAK8)[0].Pt(),wt);
      h2_LeadAK8PtMET_WZW[iHist]->Fill((*JetsAK8)[0].Pt(),MET,wt);
      h2_nJetsMET_WZW[iHist]->Fill(NJets,MET,wt);
      h2_HTMET_WZW[iHist]->Fill(HT,MET,wt);
      h_nAK8_WZW[iHist]->Fill(nAK8,wt);
      // if(nTag==0 && nMass==1){
      // 	double nak8wti = 1.;// nak8wts[3] = {0.31894, 0.295489, 0.26246315};
      // 	// if(nAK8>=4) nak8wti = nak8wts[2];
      // 	// else nak8wti = nak8wts[nAK8-2];
      // 	nak8wti = CRweightsHist->GetBinContent(CRweightsHist->FindBin((*JetsAK8)[0].Pt()));
      // 	h_nAK8_0T1M_shapeCorr_WZW->Fill(nAK8,wt*nak8wti);
      // 	h_METvBin_0T1M_shapeCorr_WZW->Fill(MET,wt*nak8wti);
      // }

      if(bestAK8J1IdxInMainColl >= 0){
	h_AK8Pt_WZW[iHist]->Fill(bestAK8J1.Pt(),wt);
	h_AK8Eta_WZW[iHist]->Fill(bestAK8J1.Eta(),wt);
	h_AK8M_WZW[iHist]->Fill((*JetsAK8_softDropMass)[bestAK8J1IdxInMainColl],wt);
	h_dPhiMETAK8_WZW[iHist]->Fill(abs(DeltaPhi(METPhi,bestAK8J1.Phi())),wt);
	h_deepWdiscr_WZW[iHist]->Fill((*JetsAK8_wDiscriminatorDeep)[bestAK8J1IdxInMainColl],wt);
      }
      else{
	h_AK8Pt_WZW[iHist]->Fill((*JetsAK8)[0].Pt(),wt);
	h_AK8Eta_WZW[iHist]->Fill((*JetsAK8)[0].Eta(),wt);
	h_AK8M_WZW[iHist]->Fill((*JetsAK8_softDropMass)[0],wt);
	h_dPhiMETAK8_WZW[iHist]->Fill(abs(DeltaPhi(METPhi,(*JetsAK8)[0].Phi())),wt);
	h_deepWdiscr_WZW[iHist]->Fill((*JetsAK8_wDiscriminatorDeep)[0],wt);
      }
	
      h_AK4PairPt_WZW[iHist]->Fill(ak4PairBosonLike.Pt(),wt);
      h_AK4PairEta_WZW[iHist]->Fill(ak4PairBosonLike.Eta(),wt);
      h_AK4PairM_WZW[iHist]->Fill(ak4PairBosonLike.M(),wt);
      h_dPhiMETAK4Pair_WZW[iHist]->Fill(abs(DeltaPhi(METPhi,ak4PairBosonLike.Phi())),wt);
      h_massAK8WAK4nonb_WZW[iHist]->Fill(topCand.M(),wt);
      if(catName=="2T2M") h_METvBin_2T2M->Fill(MET,wt);
      if(catName=="1T2M") h_METvBin_1T2M->Fill(MET,wt);
      if(catName=="1T1M") h_METvBin_1T1M->Fill(MET,wt);
    }
  }
  else{
    if(JetsAK8->size() > 0){
      if(JetsAK8->size() > 0 && (*JetsAK8)[0].Pt() > 200. && abs((*JetsAK8)[0].Eta()) < 2.0) 
      	mt = sqrt(2*(*JetsAK8)[0].Pt()*MET*(1-cos(DeltaPhi(METPhi,(*JetsAK8)[0].Phi()))));
      catName = to_string(nTag)+"T"+to_string(nMass)+"M"+to_string(nmass)+"m";
      h_EvtTypeFine->Fill(catName,wt);
      if(mt > 500){
	iHist = h_EvtTypeFine->GetXaxis()->FindBin(catName) - 1;
	if(iHist >= 17) cout<<catName<<":Overflow in EvtTypeFine WZW 2."<<endl;
      
	h_MET_WZW[iHist]->Fill(MET,wt);
	h_HT_WZW[iHist]->Fill(HT,wt);
	h_METvBin_WZW[iHist]->Fill(MET,wt);
	h_NJets_WZW[iHist]->Fill(NJets,wt);
	h_mTBest_WZW[iHist]->Fill(mt,wt);
	h_LeadAK8Pt_WZW[iHist]->Fill((*JetsAK8)[0].Pt(),wt);
	h_LeadAK8PtvBin_WZW[iHist]->Fill((*JetsAK8)[0].Pt(),wt);
	h_nAK8_WZW[iHist]->Fill(nAK8,wt);
	h2_LeadAK8PtMET_WZW[iHist]->Fill((*JetsAK8)[0].Pt(),MET,wt);
	h2_nJetsMET_WZW[iHist]->Fill(NJets,MET,wt);
	h2_HTMET_WZW[iHist]->Fill(HT,MET,wt);

	if(bestAK8J1IdxInMainColl >= 0){
	  h_AK8Pt_WZW[iHist]->Fill(bestAK8J1.Pt(),wt);
	  h_AK8Eta_WZW[iHist]->Fill(bestAK8J1.Eta(),wt);
	  h_AK8M_WZW[iHist]->Fill((*JetsAK8_softDropMass)[bestAK8J1IdxInMainColl],wt);
	  h_dPhiMETAK8_WZW[iHist]->Fill(abs(DeltaPhi(METPhi,bestAK8J1.Phi())),wt);
	  h_deepWdiscr_WZW[iHist]->Fill((*JetsAK8_wDiscriminatorDeep)[bestAK8J1IdxInMainColl],wt);
	}
	else{
	  h_AK8Pt_WZW[iHist]->Fill((*JetsAK8)[0].Pt(),wt);
	  h_AK8Eta_WZW[iHist]->Fill((*JetsAK8)[0].Eta(),wt);
	  h_AK8M_WZW[iHist]->Fill((*JetsAK8_softDropMass)[0],wt);
	  h_dPhiMETAK8_WZW[iHist]->Fill(abs(DeltaPhi(METPhi,(*JetsAK8)[0].Phi())),wt);
	  h_deepWdiscr_WZW[iHist]->Fill((*JetsAK8_wDiscriminatorDeep)[0],wt);
	}

	h_AK4PairPt_WZW[iHist]->Fill(ak4PairBosonLike.Pt(),wt);
	h_AK4PairEta_WZW[iHist]->Fill(ak4PairBosonLike.Eta(),wt);
	h_AK4PairM_WZW[iHist]->Fill(ak4PairBosonLike.M(),wt);
	h_dPhiMETAK4Pair_WZW[iHist]->Fill(abs(DeltaPhi(METPhi,ak4PairBosonLike.Phi())),wt);
	h_massAK8WAK4nonb_WZW[iHist]->Fill(topCand.M(),wt);
      }
    }
    else{
      mt = sqrt(2*ak4PairBosonLike.Pt()*MET*(1-cos(DeltaPhi(METPhi,ak4PairBosonLike.Phi()))));
      catName = to_string(nTag)+"T"+to_string(nMass)+"M"+to_string(nmass)+"m";
      h_EvtTypeFine->Fill(catName,wt);
      if(mt > 500){
	iHist = h_EvtTypeFine->GetXaxis()->FindBin(catName) - 1;
	if(iHist >= 17) cout<<catName<<":Overflow in EvtTypeFine WZW 3."<<nTagTau21<<endl;

	h_MET_WZW[iHist]->Fill(MET,wt);
	h_HT_WZW[iHist]->Fill(HT,wt);
	h_METvBin_WZW[iHist]->Fill(MET,wt);
	h_NJets_WZW[iHist]->Fill(NJets,wt);
	h_nAK8_WZW[iHist]->Fill(nAK8,wt);
	h2_nJetsMET_WZW[iHist]->Fill(NJets,MET,wt);
	h2_HTMET_WZW[iHist]->Fill(HT,MET,wt);
	h_mTBest_WZW[iHist]->Fill(mt,wt);
	h_AK4PairPt_WZW[iHist]->Fill(ak4PairBosonLike.Pt(),wt);
	h_AK4PairEta_WZW[iHist]->Fill(ak4PairBosonLike.Eta(),wt);
	h_AK4PairM_WZW[iHist]->Fill(ak4PairBosonLike.M(),wt);
	h_dPhiMETAK4Pair_WZW[iHist]->Fill(abs(DeltaPhi(METPhi,ak4PairBosonLike.Phi())),wt);
	h_massAK8WAK4nonb_WZW[iHist]->Fill(topCand.M(),wt);
      }
    }
  }
}

TString SignalReg::getEventTypeWH(){
  int nTag = 0, nMass = 0, nmass = 0, nAK8 = 0;
  int nTagH = 0, nMassH = 0, nmassH = 0;
  vector<bool> jetsAK8hasb;
  TString catName, catName0,catName1;

  for(int i=0;i<JetsAK8->size();i++){
    bool foundbInAK8 = 0;
    for(int b=0;b<bjets.size();b++){
      if(bjets[b].DeltaR((*JetsAK8)[i]) < 0.8){
	foundbInAK8 = 1;
	break;
      }
    }
    jetsAK8hasb.push_back(foundbInAK8);    
  }// AK8 loop for b content

  catname2j = "";
  vector<int> ak8Idx1;
  for(int i=0;i<JetsAK8->size();i++){
    if((*JetsAK8)[i].Pt() > 200 && abs((*JetsAK8)[i].Eta()) < 2.0){
      ak8Idx1.push_back(i);
    }
  }
  if(ak8Idx1.size()==2){
    int j1Idx = ak8Idx1[0], j2Idx = ak8Idx1[1];
    bool j1Wtag = 0, j1Wmass = 0, j1Htag = 0, j1Hmass = 0, j1nearb = jetsAK8hasb[j1Idx];
    bool j2Wtag = 0, j2Wmass = 0, j2Htag = 0, j2Hmass = 0, j2nearb = jetsAK8hasb[j2Idx];
    int j1FBH = -1, j1FBW = -1, j2FBH = -1, j2FBW = -1, aFBW = 0, aFBH = 0;
    int nwmass = 0, nhmass = 0;
    double j1mass = (*JetsAK8_softDropMass)[j1Idx], j2mass = (*JetsAK8_softDropMass)[j2Idx];
    double j1DeepWval = (*JetsAK8_wDiscriminatorDeep)[j1Idx], j2DeepWval = (*JetsAK8_wDiscriminatorDeep)[j2Idx];
    double j1Deepbbval = (*JetsAK8_deepDoubleBDiscriminatorH)[j1Idx], j2Deepbbval = (*JetsAK8_deepDoubleBDiscriminatorH)[j2Idx];
    if(!j1nearb){
      if(j1mass > massLow && j1mass < massHigh) j1Wmass = 1;
      else j1Wmass = 0;
      if(j1DeepWval > dwdisValue) j1Wtag = 1;
      else j1Wtag = 0;
    }
    else{
      if(j1mass > massLowH && j1mass < massHighH) j1Hmass = 1;
      else j1Hmass = 0;
      if(j1Deepbbval > deepDoubleBDiscriminatorValue)
	j1Htag = 1;
      else j1Htag = 0;
    }
    if(!j2nearb){
      if(j2mass > massLow && j2mass < massHigh) j2Wmass = 1;
      else j2Wmass = 0;
      if(j2DeepWval > dwdisValue) j2Wtag = 1;
      else j2Wtag = 0;
    }
    else{
      if(j2mass > massLowH && j2mass < massHighH) j2Hmass = 1;
      else j2Hmass = 0;
      if(j2Deepbbval > deepDoubleBDiscriminatorValue)
	j2Htag = 1;
      else j2Htag = 0;
    }
    if(!j1nearb){
      if(j1Wtag && j1Wmass){ j1FBW = 1; aFBW++;}
      else j1FBW = 0;
      j1FBH = -1;
      if(j1Wmass) nwmass++;
    }
    else{
      j1FBW = -1;
      if(j1Htag && j1Hmass){ j1FBH = 1; aFBH++;}
      else j1FBH = 0;
      if(j1Hmass) nhmass++;
    }
    if(!j2nearb){
      if(j2Wtag && j2Wmass){ j2FBW = 1; aFBW++;}
      else j2FBW = 0;
      j2FBH = -1;
      if(j2Wmass) nwmass++;
    }
    else{
      j2FBW = -1;
      if(j2Htag && j2Hmass){ j2FBH = 1; aFBH++;}
      else j2FBH = 0;
      if(j2Hmass) nhmass++;
    }
    
    if(j1FBH==1 && j1FBW==1) cout<<"!!!!j1FBH==1&& j1FBW==1"<<endl;
    //    if(aFBH>1 || aFBW>1) cout<<"!!!!aFBH>1 || aFBW>1"<<aFBH<<" "<<aFBW<<endl;
    if(aFBW>1 || aFBH>1) return "";
    if(nwmass>1 || nhmass >1) return "";

    if(aFBW==1 && aFBH==1) catname2j = catname2j+"FBWH_SR";
    else{
      if(aFBW==1) catname2j = catname2j+"FBW_SR";
      if(aFBH==1) catname2j = catname2j+"FBH_SR";
      if(aFBW==0 && aFBH==0){
	if(nwmass==1 && nhmass==1) catname2j = catname2j+"FBWH_CR";
	if(nwmass==1 && nhmass==0) catname2j = catname2j+"FBW_CR";
	if(nwmass==0 && nhmass==1) catname2j = catname2j+"FBH_CR";
      }
    }
    h_catname2j->Fill(catname2j,wt);
  }

  bool has2AK8 = false, fbH = 0, fbW = 0, fbZ = 0;
  double ak8Hmass = 0, ak8Wmass = 0, nMassZ =0, ddbdis = 0, dwdis = 0;
  TLorentzVector ak8H, ak8W;
  for(int i=0;i<JetsAK8->size();i++){
    if((*JetsAK8)[i].Pt() > 200 && abs((*JetsAK8)[i].Eta()) < 2.0){
      nAK8++;
      if(jetsAK8hasb[i]){
	if(ddbdis < (*JetsAK8_deepDoubleBDiscriminatorH)[i]) ddbdis = (*JetsAK8_deepDoubleBDiscriminatorH)[i];
	if(ak8H.Pt() < (*JetsAK8)[i].Pt()){ ak8H = (*JetsAK8)[i]; ak8Hmass = (*JetsAK8_softDropMass)[i];}
	if(((*JetsAK8_softDropMass)[i] > massLowH) && ((*JetsAK8_softDropMass)[i] < massHighH)){
	  nMassH++;
	  if((*JetsAK8_deepDoubleBDiscriminatorH)[i] > deepDoubleBDiscriminatorValue) fbH = true;
	}      
	if((*JetsAK8_deepDoubleBDiscriminatorH)[i] > deepDoubleBDiscriminatorValue){
	  nTagH++;
	  if(!(((*JetsAK8_softDropMass)[i] > massLowH) && ((*JetsAK8_softDropMass)[i] < massHighH))){
	    h_TaggedHMassFailM->Fill((*JetsAK8_softDropMass)[i],wt);
	    for(int j=0;isMC && j<GenParticles->size();j++){
	      if(abs((*GenParticles_PdgId)[j])==25) h_dRTaggedHFailM->Fill((*GenParticles)[j].DeltaR((*JetsAK8)[i]),wt);
	    }
	  }
	}
	if(((*JetsAK8_softDropMass)[i] > massLowZ) && ((*JetsAK8_softDropMass)[i] < massHighZ)){
	  nMassZ++;
	  if((*JetsAK8_deepDoubleBDiscriminatorH)[i] > deepDoubleBDiscriminatorValue) fbZ = true;
	}
      }
      else{
	if(dwdis < (*JetsAK8_wDiscriminatorDeep)[i]) dwdis = (*JetsAK8_wDiscriminatorDeep)[i];
	if(ak8W.Pt() < (*JetsAK8)[i].Pt()){ ak8W = (*JetsAK8)[i]; ak8Wmass = (*JetsAK8_softDropMass)[i];}
	if(((*JetsAK8_softDropMass)[i] > massLow) && ((*JetsAK8_softDropMass)[i] < massHigh)){
	  nMass++;
	  if((*JetsAK8_wDiscriminatorDeep)[i] > dwdisValue) fbW = true;
	}
	//	if( (((*JetsAK8_NsubjettinessTau2)[i])/((*JetsAK8_NsubjettinessTau1)[i])) < tau21Value){
	if( (*JetsAK8_wDiscriminatorDeep)[i] > dwdisValue){
	  nTag++;
	  if(!(((*JetsAK8_softDropMass)[i] > massLow) && ((*JetsAK8_softDropMass)[i] < massHigh))){
	    h_TaggedWMassFailM->Fill((*JetsAK8_softDropMass)[i],wt);
	    for(int j=0;isMC && j<GenParticles->size();j++){
	      if(abs((*GenParticles_PdgId)[j])==24) h_dRTaggedWFailM->Fill((*GenParticles)[j].DeltaR((*JetsAK8)[i]),wt);
	    }
	  }
	}
      }
    }
  }//for AK8
  if(nMass>1)  nMass = 1;
  if(nMassH>1) nMassH = 1;
  if(nTag>1)   nTag = 1;
  if(nTagH>1)  nTagH = 1;

  TLorentzVector ak4Hcand;
  if(bjets.size() == 2) ak4Hcand = (bjets[0] + bjets[1]);
  else if(bjets.size()>2){
    ak4Hcand = bjets[0] + bjets[1];
    for(int i=0;i<bjets.size();i++){
      for(int j=i+1;j<bjets.size();j++){
	if(abs(ak4Hcand.M()-125) > abs((bjets[i]+bjets[j]).M()-125)) ak4Hcand = bjets[i]+bjets[j];
      }
    }
  }
  if((ak4Hcand.M() > massLowH) && (ak4Hcand.M() < massHighH)) nmassH++;
  
  TLorentzVector ak4Wcand;
  for(int i=0;i<Jets_bJetTagDeepCSVprobb->size();i++){
    if((*Jets)[i].Pt() < 30 || abs((*Jets)[i].Eta()) > 2.4) continue;
    if((*Jets_bJetTagDeepCSVprobb)[i]+(*Jets_bJetTagDeepCSVprobbb)[i] > deepCSVvalue)	continue;

    for(int j=i+1;j<Jets->size();j++){
      if((*Jets)[j].Pt() < 30 || abs((*Jets)[j].Eta()) > 2.4) continue;
      if((*Jets_bJetTagDeepCSVprobb)[j]+(*Jets_bJetTagDeepCSVprobbb)[j] > deepCSVvalue) continue;
      if(abs(ak4Wcand.M()-80) > abs(((*Jets)[i]+(*Jets)[j]).M()-80)) ak4Wcand = (*Jets)[i]+(*Jets)[j];
    }
  }
  if((ak4Wcand.M() > massLow) && (ak4Wcand.M() < massHigh)) nmass++;
  
  TLorentzVector b1,b2,q1,q2;
  double mtbmin = 0, mct = 0, mtqmin = 0, mctq = 0.;
  if(bjets.size()==1){
    b1 = bjets[0];
    mtbmin = sqrt(2*b1.Pt()*MET*(1-cos(DeltaPhi(METPhi,b1.Phi()))));
  }
  else if(bjets.size()>=2){ 
    b1 = bjets[0]; 
    b2 = bjets[1];
    mtbmin = min( sqrt(2*b1.Pt()*MET*(1-cos(DeltaPhi(METPhi,b1.Phi())))),
		  sqrt(2*b2.Pt()*MET*(1-cos(DeltaPhi(METPhi,b2.Phi())))));
    mct = sqrt(2*b1.Pt()*b2.Pt()*(1 + cos(b1.DeltaPhi(b2))));
  }

  for(int j=0;j<Jets->size();j++){
    if((*Jets)[j].Pt() < 30 || abs((*Jets)[j].Eta()) > 2.4) continue;
    if((*Jets_bJetTagDeepCSVprobb)[j]+(*Jets_bJetTagDeepCSVprobbb)[j] > deepCSVvalue) continue;
    if(q1.Pt() > 1. && q2.Pt() > 1.0) break;
    else if(q1.Pt() <= 1.0) q1 = (*Jets)[j];
    else if(q2.Pt() <= 1.0) q2 = (*Jets)[j];    
  }
  
  if(q1.Pt() > 1.0 && q2.Pt() < 1.0) mtqmin = sqrt(2*q1.Pt()*MET*(1-cos(DeltaPhi(METPhi,q1.Phi()))));
  else if(q1.Pt() > 1.0 && q2.Pt() > 1.0){
    mtqmin = min( sqrt(2*q1.Pt()*MET*(1-cos(DeltaPhi(METPhi,q1.Phi())))),
                  sqrt(2*q2.Pt()*MET*(1-cos(DeltaPhi(METPhi,q2.Phi())))));
    mctq = sqrt(2*q1.Pt()*q2.Pt()*(1 + cos(q1.DeltaPhi(q2))));
  }
  

  // if(JetsAK8->size()>=1 && (*JetsAK8)[0].Pt() > 200 && abs((*JetsAK8)[0].Eta()) < 2.0)
  //   catName0 = "1-"+to_string(nmass)+"Wm"+to_string(nmassH)+"Hm";
  // else  catName0 = "0-"+to_string(nmass)+"Wm"+to_string(nmassH)+"Hm";

  catName0 = "0-"+to_string(nmass)+"Wm"+to_string(nmassH)+"Hm";

  int iHist = 15 + h_EvtTypeWH_0AK8M->GetXaxis()->FindBin(catName0);
  if(iHist < 24){
    double mt=0,mt2j=0, mtLead = 0, mt2Lead = 0.;
    mt = sqrt(2*ak4Hcand.Pt()*MET*(1-cos(DeltaPhi(METPhi,ak4Hcand.Phi()))));
    mt2j = sqrt(2*ak4Wcand.Pt()*MET*(1-cos(DeltaPhi(METPhi,ak4Wcand.Phi()))));

    //    if(JetsAK8->size() > 0 && (*JetsAK8)[0].Pt() > 200 && abs((*JetsAK8)[0].Eta()) < 2.0){ mtLead = sqrt(2*(*JetsAK8)[0].Pt()*MET*(1-cos(DeltaPhi(METPhi,(*JetsAK8)[0].Phi())))); }
    if(ak4Hcand.Pt() > ak4Wcand.Pt()){ mtLead = mt; mt2Lead = mt2j;}
    else if(ak4Hcand.Pt() < ak4Wcand.Pt()){ mtLead = mt2j; mt2Lead = mt;}
    
    if(nmass==0 && nmassH==0 && NJets >=2) mtLead = sqrt(2*(((*Jets)[0]+(*Jets)[1]).Pt())*MET*(1-cos(DeltaPhi(METPhi,((*Jets)[0]+(*Jets)[1]).Phi()))));
    
    if(mtLead > 500. && mtbmin > 200 && (nMass+nMassH+nTag+nTagH < 4)){
      h_EvtTypeWH_0AK8M->Fill(catName0,wt);

      h_MET_catWH[iHist]->Fill(MET,wt);
      h_HT_catWH[iHist]->Fill(HT,wt);
      h_METvBin_catWH[iHist]->Fill(MET,wt);
      h_mT_catWH[iHist]->Fill(mt,wt);
      h_mTLead_catWH[iHist]->Fill(mtLead,wt);
      h_mT2Lead_catWH[iHist]->Fill(mt2Lead,wt);
      h_AK8HPt_catWH[iHist]->Fill(ak4Hcand.Pt(),wt);
      h_AK8HEta_catWH[iHist]->Fill(ak4Hcand.Eta(),wt);
      h_AK8HM_catWH[iHist]->Fill(ak4Hcand.M(),wt);
      h_dPhiMETAK8H_catWH[iHist]->Fill(abs(DeltaPhi(METPhi,ak4Hcand.Phi())),wt);
      h_AK8WPt_catWH[iHist]->Fill(ak4Wcand.Pt(),wt);
      h_AK8WEta_catWH[iHist]->Fill(ak4Wcand.Eta(),wt);
      h_AK8WM_catWH[iHist]->Fill(ak4Wcand.M(),wt);
      h_mT2J_catWH[iHist]->Fill(mt2j,wt);
      h_nAK8_catWH[iHist]->Fill(nAK8,wt);

      if(bjets.size() >=2) h_LeadbPairMass_catWH[iHist]->Fill((bjets[0]+bjets[1]).M(),wt);
      else h_LeadbPairMass_catWH[iHist]->Fill(0.,wt);
      if(nonbjets.size() >=2) h_LeadNonbPairMass_catWH[iHist]->Fill((nonbjets[0]+nonbjets[1]).M(),wt);
      else h_LeadNonbPairMass_catWH[iHist]->Fill(0.,wt);
      h_deepDoubleBdiscr_catWH[iHist]->Fill(ddbdis,wt);
      h_deepWdiscr_catWH[iHist]->Fill(dwdis,wt);

      h_mEfft_catWH[iHist]->Fill(MET+HT,wt);;
      h_mTbMin_catWH[iHist]->Fill(mtbmin,wt);
      h_mCT_catWH[iHist]->Fill(mct,wt);
      h_mTqMin_catWH[iHist]->Fill(mtqmin,wt);
      h_mCTq_catWH[iHist]->Fill(mctq,wt);
    }
  }
  else cout<<"Oveflow in cat for WH (0AK8)"<<endl;
  
  
  catName = to_string(nTag)+"Wt"+to_string(nMass)+"Wm"+to_string(nTagH)+"Ht"+to_string(nMassH)+"Hm";
  //  if(nTagH!=0 || nMassH==0 || (nTag>=1 && nMass>=1)){ cout<<"nHm:"<<nMassH<<" nWm:"<<nMass<<" nHT:"<<nTagH<<" nWT:"<<nTag<<endl; print(0);}
  double mt=0,mt2j=0, mtLead = 0, mt2Lead = 0.;
  mt = sqrt(2*ak8H.Pt()*MET*(1-cos(DeltaPhi(METPhi,ak8H.Phi()))));
  mt2j = sqrt(2*ak8W.Pt()*MET*(1-cos(DeltaPhi(METPhi,ak8W.Phi()))));

  if(JetsAK8->size() > 0 && (*JetsAK8)[0].Pt() > 200 && abs((*JetsAK8)[0].Eta()) < 2.0){ mtLead = sqrt(2*(*JetsAK8)[0].Pt()*MET*(1-cos(DeltaPhi(METPhi,(*JetsAK8)[0].Phi())))); }
  if(JetsAK8->size() > 1 && (*JetsAK8)[1].Pt() > 200 && abs((*JetsAK8)[1].Eta()) < 2.0){ mt2Lead = sqrt(2*(*JetsAK8)[1].Pt()*MET*(1-cos(DeltaPhi(METPhi,(*JetsAK8)[1].Phi())))); }
    
  iHist = h_EvtTypeWH->GetXaxis()->FindBin(catName)-1;
  if(iHist < 16){
    if((catName=="0Wt0Wm0Ht0Hm" || mtLead > 500.) && mtbmin > 200){
      h_EvtTypeWH->Fill(catName,wt);

      h_MET_catWH[iHist]->Fill(MET,wt);
      h_HT_catWH[iHist]->Fill(HT,wt);
      h_METvBin_catWH[iHist]->Fill(MET,wt);
      h_mT_catWH[iHist]->Fill(mt,wt);
      h_AK8HPt_catWH[iHist]->Fill(ak8H.Pt(),wt);
      h_AK8HEta_catWH[iHist]->Fill(ak8H.Eta(),wt);
      h_AK8HM_catWH[iHist]->Fill(ak8Hmass,wt);
      h_dPhiMETAK8H_catWH[iHist]->Fill(abs(DeltaPhi(METPhi,ak8H.Phi())),wt);
      h_AK8WPt_catWH[iHist]->Fill(ak8W.Pt(),wt);
      h_AK8WEta_catWH[iHist]->Fill(ak8W.Eta(),wt);
      h_AK8WM_catWH[iHist]->Fill(ak8Wmass,wt);
      h_mT2J_catWH[iHist]->Fill(mt2j,wt);
      h_nAK8_catWH[iHist]->Fill(nAK8,wt);

      if(bjets.size() >=2) h_LeadbPairMass_catWH[iHist]->Fill((bjets[0]+bjets[1]).M(),wt);
      else h_LeadbPairMass_catWH[iHist]->Fill(0.,wt);
      if(nonbjets.size() >=2) h_LeadNonbPairMass_catWH[iHist]->Fill((nonbjets[0]+nonbjets[1]).M(),wt);
      else h_LeadNonbPairMass_catWH[iHist]->Fill(0.,wt);

      h_deepDoubleBdiscr_catWH[iHist]->Fill(ddbdis,wt);
      h_deepWdiscr_catWH[iHist]->Fill(dwdis,wt);
      
      h_mTLead_catWH[iHist]->Fill(mtLead,wt);
      h_mT2Lead_catWH[iHist]->Fill(mt2Lead,wt);
      h_mEfft_catWH[iHist]->Fill(MET+HT,wt);
      h_mTbMin_catWH[iHist]->Fill(mtbmin,wt);
      h_mCT_catWH[iHist]->Fill(mct,wt);
      h_mTqMin_catWH[iHist]->Fill(mtqmin,wt);
      h_mCTq_catWH[iHist]->Fill(mctq,wt);
      if(fbW && fbH) h_METvBin_FBWH->Fill(MET,wt);
      if(fbW && fbZ) h_METvBin_FBWZ->Fill(MET,wt);
      if(fbW && !fbH && !fbZ) h_METvBin_FBW->Fill(MET,wt);
      if(!fbW && fbH) h_METvBin_FBH->Fill(MET,wt);
      if(!fbW && fbZ) h_METvBin_FBZ->Fill(MET,wt);
      //      if(fbH==1 && fbZ==1) cout<<"!!!! doubleB mass windows overlap"<<endl;
    }
    if(nMass==0 || nMassH==0){
      catName1 = to_string(nTag)+"Wt"+to_string(nMass)+"wm"+to_string(nTagH)+"Ht"+to_string(nMassH)+"hm";
      double mt=0,mt2j=0, mtLead = 0, mt2Lead = 0.;
      if(nMass == 0) ak8Wmass = ak4Wcand.M();
      else ak4Wcand = ak8W;
      if(nMassH== 0) ak8Hmass = ak4Hcand.M();
      else ak4Hcand = ak8H;

      mt = sqrt(2*ak4Hcand.Pt()*MET*(1-cos(DeltaPhi(METPhi,ak4Hcand.Phi()))));
      mt2j = sqrt(2*ak4Wcand.Pt()*MET*(1-cos(DeltaPhi(METPhi,ak4Wcand.Phi()))));
      
      if(JetsAK8->size() > 0 && (*JetsAK8)[0].Pt() > 200 && abs((*JetsAK8)[0].Eta()) < 2.0){ mtLead = sqrt(2*(*JetsAK8)[0].Pt()*MET*(1-cos(DeltaPhi(METPhi,(*JetsAK8)[0].Phi())))); }
      else if(ak4Hcand.Pt() > ak4Wcand.Pt()){ mtLead = mt; mt2Lead = mt2j;}
      else if(ak4Hcand.Pt() < ak4Wcand.Pt()){ mtLead = mt2j; mt2Lead = mt;}
      
      if(mtLead > 500. && mtbmin > 200 && (ak8Wmass > massLow) && (ak8Wmass < massHigh) && (ak8Hmass > massLowH) && (ak8Hmass < massHighH) ){
	iHist = 23 + h_EvtTypeWH_ak84->GetXaxis()->FindBin(catName1);
	if(iHist < 36){
	  h_EvtTypeWH_ak84->Fill(catName1,wt);
	  
	  h_MET_catWH[iHist]->Fill(MET,wt);
	  h_HT_catWH[iHist]->Fill(HT,wt);
	  h_METvBin_catWH[iHist]->Fill(MET,wt);
	  h_mT_catWH[iHist]->Fill(mt,wt);
	  h_AK8HPt_catWH[iHist]->Fill(ak4Hcand.Pt(),wt);
	  h_AK8HEta_catWH[iHist]->Fill(ak4Hcand.Eta(),wt);
	  h_AK8HM_catWH[iHist]->Fill(ak8Hmass,wt);
	  h_dPhiMETAK8H_catWH[iHist]->Fill(abs(DeltaPhi(METPhi,ak4Hcand.Phi())),wt);
	  h_AK8WPt_catWH[iHist]->Fill(ak4Wcand.Pt(),wt);
	  h_AK8WEta_catWH[iHist]->Fill(ak4Wcand.Eta(),wt);
	  h_AK8WM_catWH[iHist]->Fill(ak8Wmass,wt);
	  h_mT2J_catWH[iHist]->Fill(mt2j,wt);
	  h_nAK8_catWH[iHist]->Fill(nAK8,wt);

	  if(bjets.size() >=2) h_LeadbPairMass_catWH[iHist]->Fill((bjets[0]+bjets[1]).M(),wt);
	  else h_LeadbPairMass_catWH[iHist]->Fill(0.,wt);
	  if(nonbjets.size() >=2) h_LeadNonbPairMass_catWH[iHist]->Fill((nonbjets[0]+nonbjets[1]).M(),wt);
	  else h_LeadNonbPairMass_catWH[iHist]->Fill(0.,wt);

	  h_deepDoubleBdiscr_catWH[iHist]->Fill(ddbdis,wt);
	  h_deepWdiscr_catWH[iHist]->Fill(dwdis,wt);
	  
	  h_mTLead_catWH[iHist]->Fill(mtLead,wt);
	  h_mT2Lead_catWH[iHist]->Fill(mt2Lead,wt);
	  h_mEfft_catWH[iHist]->Fill(MET+HT,wt);
	  h_mTbMin_catWH[iHist]->Fill(mtbmin,wt);
	  h_mCT_catWH[iHist]->Fill(mct,wt);
	  h_mTqMin_catWH[iHist]->Fill(mtqmin,wt);
	  h_mCTq_catWH[iHist]->Fill(mctq,wt);	  
	}
	else cout<<"Overflow in cat for WH with 0Hm or 0Wm"<<endl;
      }
    }      
  }
  else cout<<"Oveflow in cat for WH: "<<catName<<endl;
  return catName+catName0;
}

void SignalReg::applySDmassCorrAllAK8(){
  float genCorr  = 1.;
  float recoCorr = 1.;
  float totalWeight = 1.;

  TLorentzVector subjetSum;
  uncorrAK8SubjSum.resize(0);
  SDmassCorrFac.resize(0);

  for(int i=0;i< JetsAK8->size();i++){
    for(int j=0;j < (*JetsAK8_subjets)[i].size();j++){
      if(j==0) subjetSum = (*JetsAK8_subjets)[i][j];
      else subjetSum = subjetSum + (*JetsAK8_subjets)[i][j];
      //      cout<<"subj pt,eta,phi,M:"<<(*JetsAK8_subjets)[i][j].Pt()<<", "<<(*JetsAK8_subjets)[i][j].Eta()<<", "<<(*JetsAK8_subjets)[i][j].Phi()<<" ,"<<(*JetsAK8_subjets)[i][j].M()<<endl;
    }
    uncorrAK8SubjSum.push_back(subjetSum);
    SDmassCorrFac.push_back(1);
    
    genCorr =  puppisd_corrGEN->Eval( (*JetsAK8)[i].Pt() );
    if( abs((*JetsAK8)[i].Eta())  <= 1.3 ) 
      recoCorr = puppisd_corrRECO_cen->Eval( (*JetsAK8)[i].Pt() );
    else 
      recoCorr = puppisd_corrRECO_for->Eval( (*JetsAK8)[i].Pt() );
    
    SDmassCorrFac.push_back(genCorr * recoCorr);
    (*JetsAK8_softDropMass)[i] = (subjetSum.M())*genCorr*recoCorr;
    //    cout<<"genCorr:"<<genCorr<<" recoCorr:"<<recoCorr<<" subjetSum.M():"<<subjetSum.M()<<endl;
  }
}
bool SignalReg::passHEMjetVeto(double ptThresh) {
  for (int p = 0; p < Jets->size(); p++){
    if (-3.2 <= Jets->at(p).Eta() && Jets->at(p).Eta() <= -1.2 &&
	-1.77 <= Jets->at(p).Phi() && Jets->at(p).Phi() <= -0.67 &&
	Jets->at(p).Pt() > ptThresh && abs(DeltaPhi(Jets->at(p).Phi(),METPhi)) < 0.5)
      return false;
  }
  return true;
}
void SignalReg::print(Long64_t jentry){
  //cout<<endl;
  TLorentzVector v1,photo;
  cout<<"------------------------------------------------------------"<<endl;
  cout<<"MomMass:"<<SusyMotherMass<<" Kid Mass:"<<SusyLSPMass<<endl;
  for(int i=0;isMC && i<GenParticles->size();i++){  
    cout<<EvtNum<<" "<<jentry<<" "<<GenParticles->size()<<" "<<i<<" PdgId:"<<(*GenParticles_PdgId)[i]<<" parentId:"<<(*GenParticles_ParentId)[i]<<" parentIndx:"<<(*GenParticles_ParentIdx)[i]<<" Status:"<<(*GenParticles_Status)[i]<</*"\tPx:"<<(*GenParticles)[i].Px()<<" Py:"<<(*GenParticles)[i].Py()<<" Pz:"<<(*GenParticles)[i].Pz()<<*/"\tPt:"<<(*GenParticles)[i].Pt()<<" Eta:"<<(*GenParticles)[i].Eta()<<" Phi:"<<(*GenParticles)[i].Phi()<<" M:"<<(*GenParticles)[i].M()<<endl;

  }
  for(int i=0;i<GenJets->size();i++){
    cout<<"GenJetPt:"<<(*GenJets)[i].Pt()<<" GenJetEta:"<<(*GenJets)[i].Eta()<<" GenJetPhi:"<<(*GenJets)[i].Phi()<<endl;
  }
  for(int i=0;i<Jets->size();i++){
    cout<<"JetPt:"<<(*Jets)[i].Pt()<<" JetEta:"<<(*Jets)[i].Eta()<<" JetPhi:"<<(*Jets)[i].Phi()<<" DeepCSV:"<<(*Jets_bJetTagDeepCSVprobb)[i]+(*Jets_bJetTagDeepCSVprobbb)[i]<<endl;
  }
  cout<<"MET:"<<MET<<" METPhi:"<<METPhi<<" MHTPhi:"<<MHTPhi<<" DPhi1:"<<DeltaPhi1<<" DeltaPhi2:"<<DeltaPhi2<<" DeltaPhi3:"<<DeltaPhi3<<" DeltaPhi4:"<<DeltaPhi4<<endl;
  for(int i=0;i<JetsAK8->size();i++){
    cout<<"AK8 pT:"<<(*JetsAK8)[i].Pt()<<" eta:"<<(*JetsAK8)[i].Eta()<<" phi:"<<(*JetsAK8)[i].Phi()<<" softDropM:"<<(*JetsAK8_softDropMass)[i]<<" tau21:"<<(*JetsAK8_NsubjettinessTau2)[i]/(*JetsAK8_NsubjettinessTau1)[i]<<" DeepDoubleBDiscr:"<<(*JetsAK8_deepDoubleBDiscriminatorH)[i]<<" WDiscr:"<<(*JetsAK8_wDiscriminatorDeep)[i]<<" WdiscrDecorr:"<<(*JetsAK8_wDiscriminatorDeepDecorrel)[i]<<endl;
  }
}
