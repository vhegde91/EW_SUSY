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
  if(s_data.Contains("MC_2016")){ dataRun = -2016; lumiInfb = 35.815165; deepCSVvalue = 0.6321;}
  else if(s_data.Contains("MC_2017")){ dataRun = -2017; lumiInfb = 41.486136; deepCSVvalue = 0.4941;}
  else if(s_data.Contains("MC_2018")){ dataRun = -2018; lumiInfb = 59.546381; deepCSVvalue = 0.4184;}
  
  else if(s_data.Contains("2016")){ dataRun = 2016; isMC = false; deepCSVvalue = 0.6321;}
  else if(s_data.Contains("2017")){ dataRun = 2017; isMC = false; deepCSVvalue = 0.4941;}
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
    
    nonbjets.resize(0);
    bjets.resize(0);
    for(int i=0;i<Jets_bJetTagDeepCSVprobb->size();i++){
      if((*Jets)[i].Pt() < 30 || abs((*Jets)[i].Eta()) > 2.4) continue;
      if((*Jets_bJetTagDeepCSVprobb)[i]+(*Jets_bJetTagDeepCSVprobbb)[i] > deepCSVvalue)
	bjets.push_back((*Jets)[i]);
      else nonbjets.push_back((*Jets)[i]);
    }
    sortTLorVec(&bjets);
    
    TString WH_catName;// = getEventTypeWH();
    if(BTagsDeepCSV != 0) continue;
    int evtType = getEventType();
    getEventTypeFine(evtType);//to be called only after calling getEventType();
    getEventTypeWZW();

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
    for(int i=0;i<JetsAK8->size();i++){
      for(int j=0;j<genZ.size();j++)
	if(genZ[j].DeltaR((*JetsAK8)[i]) < 0.3) h_AK8MassNearGenZ->Fill((*JetsAK8_softDropMass)[i],wt);
      for(int j=0;j<genW.size();j++)
	if(genW[j].DeltaR((*JetsAK8)[i]) < 0.3) h_AK8MassNearGenW->Fill((*JetsAK8_softDropMass)[i],wt);
      for(int j=0;j<higgs.size();j++)
	if(higgs[j].DeltaR((*JetsAK8)[i]) < 0.3) h_AK8MassNearGenH->Fill((*JetsAK8_softDropMass)[i],wt);      
    }
    
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
      h2_mTbMin_mCT->Fill(mtbmin,mct,wt);
    }
    h_mTbMin->Fill(mtbmin,wt);
    h_mCT->Fill(mct,wt);

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
    //    cout<<goodAk8.size()<<" :"<<ak8JBoosted.size()<<" nBoosted:"<<goodAk8.size()<<endl;
    //------------ZZMET-----------
    // if(MET < 300 || HT < 400) continue;//ZZMET
    // if(goodAk8.size()<2) continue;//ZZMET
    // if(goodAk8[0].Pt() < 200) continue;//ZZMET
    // if(goodAk8[1].Pt() < 200) continue;//ZZMET
    // if(((*JetsAK8_softDropMass)[0] < 70 || (*JetsAK8_softDropMass)[0] > 100)) continue;//ZZMET
    // if(((*JetsAK8_softDropMass)[1] < 70 || (*JetsAK8_softDropMass)[1] > 100)) continue;//ZZMET
        
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
	    if(minMassDiff > (80 - (ak4jNotAK8[k]+ak4jNotAK8[j]).M())){
	      minMassDiff = 80 - (ak4jNotAK8[k]+ak4jNotAK8[j]).M();
	      ak4PairBosonLike = ak4jNotAK8[k]+ak4jNotAK8[j];
	    }
	  }
	}
      }
      if(ak4jNotAK8.size() < 2) continue;
      if(ak4PairBosonLike.M() < 60 || ak4PairBosonLike.M() > 120) continue;
    }
    //    if(evtType!=200) continue;
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

    if(evtType==0){
      h_EvtType->Fill("0 Good AK8",wt);
      continue;
    }
    else if(evtType==200) h_EvtType->Fill("2 Boosted",wt);
    else if(evtType==100 || evtType==150) h_EvtType->Fill("1 Boosted",wt);
    else if(evtType==1) h_EvtType->Fill("1 Good AK8",wt);
    else if(evtType==50) h_EvtType->Fill("2 Good AK8",wt);
    if(!(evtType==100 || evtType==150)) continue;
    if(mt < 500) continue;

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
      h_AK8J1MainWdisc->Fill(mainWdisc[bestAK8J1IdxInMainColl],wt);
      h_MT->Fill(mt,wt);
      h_MTvBin->Fill(mt,wt);
      h_dPhiMETAK8->Fill(abs(DeltaPhi(METPhi,bestAK8J1.Phi())),wt);
    }
    if(bestAK8J2IdxInMainColl != -1){
      h_AK8J2Pt->Fill((bestAK8J2.Pt()),wt);
      h_AK8J2Eta->Fill((bestAK8J2.Eta()),wt);
      h_AK8J2Mass->Fill((*JetsAK8_softDropMass)[bestAK8J2IdxInMainColl],wt);
      h_AK8J2MainWdisc->Fill(mainWdisc[bestAK8J2IdxInMainColl],wt);
      h_dPhiAK8J1J2->Fill(abs(DeltaPhi(bestAK8J1.Phi(),bestAK8J2.Phi())),wt);
      h_InvMassAK8Jets->Fill((bestAK8J1+bestAK8J2).M(),wt);
      h_AK8J1J2MassRatio->Fill((*JetsAK8_softDropMass)[bestAK8J2IdxInMainColl]/(*JetsAK8_softDropMass)[bestAK8J1IdxInMainColl],wt);
      h_MT2J->Fill(mt2j,wt);
      h_MT2JvBin->Fill(mt2j,wt);
      h_mTSum->Fill(mt+mt2j,wt);
      h_mTSumvBin->Fill(mt+mt2j,wt);
      h_mTRatio->Fill(mt/mt2j,wt);

      h2_AK8J1J2MainWdisc->Fill(mainWdisc[bestAK8J1IdxInMainColl],mainWdisc[bestAK8J2IdxInMainColl],wt);
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
  else if(mainWdisc.size()!=0){
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
    if(bestAK8J2IdxInMainColl >= 1) evtType =  50;//there are 2 good AK8 jets (both failing MainWdisc). 0T2M
    else evtType =  1;//there is only one good ak8. 0T1M
  }
  return evtType;  
}

void SignalReg::getEventTypeWZW(){
  if(BTagsDeepCSV > 0) return;
  vector<TString> name;
  int bestWidx = -1, nWmd = 0, nZmd = 0, nWmc = 0, nZmc = 0, nTau21 = 0;
  double mt = 0., sdMass_Wmd = 0., sdMass_Zmd = 0., sdMass_Wmc = 0., sdMass_Zmc = 0., sdMass_tau21 = 0., sdMass = 0.;
  for(int i=0; i < JetsAK8->size(); i++){
    if(mt < 0.1) mt = sqrt(2*(*JetsAK8)[i].Pt()*MET*(1-cos(DeltaPhi(METPhi,(*JetsAK8)[i].Phi()))));
    sdMass = (*JetsAK8_softDropMass)[i];
    if(sdMass < massLow || sdMass > massHigh) continue;
    if((*JetsAK8)[i].Pt() < 200 || abs((*JetsAK8)[i].Eta()) > 2.0) continue;
    if((*JetsAK8_wDiscriminatorDeep)[i] > dwdisValue && bestWidx < 0) bestWidx = i;
  }
  if(bestWidx < 0) return;
  if(mt < 500) return;

  for(int i=0; i < JetsAK8->size(); i++){
    if((*JetsAK8)[i].Pt() < 200 || abs((*JetsAK8)[i].Eta()) > 2.0 || i==bestWidx) continue;
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
  }
  // if(nWmc){
  //   print(0);
  //   cout<<"bestWidx:"<<bestWidx<<" nWmc:"<<nWmc<<" nWmd:"<<nWmd<<" tau21:"<<nTau21<<" sdMass_Wmc:"<<sdMass_Wmc<<" sdMass_Wmd:"<<sdMass_Wmd<<" sdMass_tau21:"<<sdMass_tau21<<endl;
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
}

TString SignalReg::getEventTypeWH(){
  int nTag = 0, nMass = 0, nmass = 0;
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
  
  bool has2AK8 = false;
  double ak8Hmass = 0, ak8Wmass = 0, ddbdis = 0, dwdis = 0;
  TLorentzVector ak8H, ak8W;
  for(int i=0;i<JetsAK8->size();i++){
    if((*JetsAK8)[i].Pt() > 200 && abs((*JetsAK8)[i].Eta()) < 2.0){
      if(jetsAK8hasb[i]){
	if(ddbdis < (*JetsAK8_deepDoubleBDiscriminatorH)[i]) ddbdis = (*JetsAK8_deepDoubleBDiscriminatorH)[i];
	if(ak8H.Pt() < (*JetsAK8)[i].Pt()){ ak8H = (*JetsAK8)[i]; ak8Hmass = (*JetsAK8_softDropMass)[i];}
	if(((*JetsAK8_softDropMass)[i] > massLowH) && ((*JetsAK8_softDropMass)[i] < massHighH)){
	  nMassH++;
	}
	if((*JetsAK8_deepDoubleBDiscriminatorH)[i] > deepDoubleBDiscriminatorValue){
	  nTagH++;
	  if(!(((*JetsAK8_softDropMass)[i] > massLowH) && ((*JetsAK8_softDropMass)[i] < massHighH))){
	    h_TaggedHMassFailM->Fill((*JetsAK8_softDropMass)[i],wt);
	    for(int j=0;j<GenParticles->size();j++){
	      if(abs((*GenParticles_PdgId)[j])==25) h_dRTaggedHFailM->Fill((*GenParticles)[j].DeltaR((*JetsAK8)[i]),wt);
	    }
	  }
	}
      }
      else{
	if(dwdis < (*JetsAK8_wDiscriminatorDeep)[i]) dwdis = (*JetsAK8_wDiscriminatorDeep)[i];
	if(ak8W.Pt() < (*JetsAK8)[i].Pt()){ ak8W = (*JetsAK8)[i]; ak8Wmass = (*JetsAK8_softDropMass)[i];}
	if(((*JetsAK8_softDropMass)[i] > massLow) && ((*JetsAK8_softDropMass)[i] < massHigh)){
	  nMass++;
	}
	//	if( (((*JetsAK8_NsubjettinessTau2)[i])/((*JetsAK8_NsubjettinessTau1)[i])) < tau21Value){
	if( (*JetsAK8_wDiscriminatorDeep)[i] > dwdisValue){
	  nTag++;
	  if(!(((*JetsAK8_softDropMass)[i] > massLow) && ((*JetsAK8_softDropMass)[i] < massHigh))){
	    h_TaggedWMassFailM->Fill((*JetsAK8_softDropMass)[i],wt);
	    for(int j=0;j<GenParticles->size();j++){
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
  // cout<<"nHm:"<<nMassH<<" nWm:"<<nMass<<" nHT:"<<nTagH<<" nWT:"<<nTag<<endl;
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

void SignalReg::getEventTypeFine(int evtType){
  int nTag = 0, nMass = 0, nmass = 0, nTagTau21 = 0;
  bool has2AK8 = false;
  if(JetsAK8->size()==0){
    nTag = 0;
    nMass = 0;
    nTagTau21 = 0;
  }
  else if(JetsAK8->size()>0){
    for(int i=0;i<JetsAK8->size();i++){
      if((*JetsAK8)[i].Pt() > 200 && abs((*JetsAK8)[i].Eta()) < 2.0){
	if(((*JetsAK8_softDropMass)[i] > massLow) && ((*JetsAK8_softDropMass)[i] < massHigh))
	  nMass++;
	//	if( (((*JetsAK8_NsubjettinessTau2)[i])/((*JetsAK8_NsubjettinessTau1)[i])) < mainWdiscValue)
	if( ((*JetsAK8_wDiscriminatorDeep)[i] > dwdisValue))
	  nTag++;
	if(bestAK8J1IdxInMainColl >=0 && bestAK8J1IdxInMainColl != i && (((*JetsAK8_NsubjettinessTau2)[i])/((*JetsAK8_NsubjettinessTau1)[i]) < tau21Value))
	  nTagTau21++;
	if(i>0) has2AK8 = true;
      }
    }
    if(nTag>=1)
      nTag = (nTagTau21 == 0) ? 1 : 2;
    if(nMass>=3) nMass = 2;
  }
  //--------- AK4 ---------
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

  if(has2AK8){
    if(bestAK8J1IdxInMainColl >=0){
      catName = to_string(nTag)+"T"+to_string(nMass)+"M";
      h_EvtTypeFine->Fill(catName,wt);
      if(JetsAK8->size() > 0 && (*JetsAK8)[0].Pt() > 200. && abs((*JetsAK8)[0].Eta()) < 2.0) 
      	mt = sqrt(2*(*JetsAK8)[0].Pt()*MET*(1-cos(DeltaPhi(METPhi,(*JetsAK8)[0].Phi()))));
      //      mt = sqrt(2*bestAK8J1.Pt()*MET*(1-cos(DeltaPhi(METPhi,bestAK8J1.Phi()))));

      if(mt > 500){
	iHist = h_EvtTypeFine->GetXaxis()->FindBin(catName) - 1;
	if(iHist >= 17) cout<<catName<<":Overflow in EvtTypeFine WZW."<<endl;

	h_MET_WZW[iHist]->Fill(MET,wt);
	h_METvBin_WZW[iHist]->Fill(MET,wt);
	h_NJets_WZW[iHist]->Fill(NJets,wt);
	h_mTBest_WZW[iHist]->Fill(mt,wt);
	h_AK8Pt_WZW[iHist]->Fill(bestAK8J1.Pt(),wt);
	h_AK8Eta_WZW[iHist]->Fill(bestAK8J1.Eta(),wt);
	h_AK8M_WZW[iHist]->Fill((*JetsAK8_softDropMass)[bestAK8J1IdxInMainColl],wt);
	h_dPhiMETAK8_WZW[iHist]->Fill(abs(DeltaPhi(METPhi,bestAK8J1.Phi())),wt);
	h_deepWdiscr_WZW[iHist]->Fill((*JetsAK8_wDiscriminatorDeep)[bestAK8J1IdxInMainColl],wt);

	h_AK4PairPt_WZW[iHist]->Fill(ak4PairBosonLike.Pt(),wt);
	h_AK4PairEta_WZW[iHist]->Fill(ak4PairBosonLike.Eta(),wt);
	h_AK4PairM_WZW[iHist]->Fill(ak4PairBosonLike.M(),wt);
	h_dPhiMETAK4Pair_WZW[iHist]->Fill(abs(DeltaPhi(METPhi,ak4PairBosonLike.Phi())),wt);	
      }
    }
  }
  else{
    if(bestAK8J1IdxInMainColl >= 0){
      if(JetsAK8->size() > 0 && (*JetsAK8)[0].Pt() > 200. && abs((*JetsAK8)[0].Eta()) < 2.0) 
      	mt = sqrt(2*(*JetsAK8)[0].Pt()*MET*(1-cos(DeltaPhi(METPhi,(*JetsAK8)[0].Phi()))));
      catName = to_string(nTag)+"T"+to_string(nMass)+"M"+to_string(nmass)+"m";
      h_EvtTypeFine->Fill(catName,wt);
      if(mt > 500){
	iHist = h_EvtTypeFine->GetXaxis()->FindBin(catName) - 1;
	if(iHist >= 17) cout<<catName<<":Overflow in EvtTypeFine WZW 2."<<endl;
      
	h_MET_WZW[iHist]->Fill(MET,wt);
	h_METvBin_WZW[iHist]->Fill(MET,wt);
	h_NJets_WZW[iHist]->Fill(NJets,wt);
	h_mTBest_WZW[iHist]->Fill(mt,wt);
	h_AK8Pt_WZW[iHist]->Fill(bestAK8J1.Pt(),wt);
	h_AK8Eta_WZW[iHist]->Fill(bestAK8J1.Eta(),wt);
	h_AK8M_WZW[iHist]->Fill((*JetsAK8_softDropMass)[bestAK8J1IdxInMainColl],wt);
	h_dPhiMETAK8_WZW[iHist]->Fill(abs(DeltaPhi(METPhi,bestAK8J1.Phi())),wt);
	h_deepWdiscr_WZW[iHist]->Fill((*JetsAK8_wDiscriminatorDeep)[bestAK8J1IdxInMainColl],wt);

	h_AK4PairPt_WZW[iHist]->Fill(ak4PairBosonLike.Pt(),wt);
	h_AK4PairEta_WZW[iHist]->Fill(ak4PairBosonLike.Eta(),wt);
	h_AK4PairM_WZW[iHist]->Fill(ak4PairBosonLike.M(),wt);
	h_dPhiMETAK4Pair_WZW[iHist]->Fill(abs(DeltaPhi(METPhi,ak4PairBosonLike.Phi())),wt);
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
	h_METvBin_WZW[iHist]->Fill(MET,wt);
	h_NJets_WZW[iHist]->Fill(NJets,wt);
	h_mTBest_WZW[iHist]->Fill(mt,wt);
	h_AK4PairPt_WZW[iHist]->Fill(ak4PairBosonLike.Pt(),wt);
	h_AK4PairEta_WZW[iHist]->Fill(ak4PairBosonLike.Eta(),wt);
	h_AK4PairM_WZW[iHist]->Fill(ak4PairBosonLike.M(),wt);
	h_dPhiMETAK4Pair_WZW[iHist]->Fill(abs(DeltaPhi(METPhi,ak4PairBosonLike.Phi())),wt);
      }
    }
  }
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
  for(int i=0;i<GenJets->size();i++){
    cout<<"JetPt:"<<(*GenJets)[i].Pt()<<" JetEta:"<<(*GenJets)[i].Eta()<<" JetPhi:"<<(*GenJets)[i].Phi()<<endl;
  }
  for(int i=0;i<Jets->size();i++){
    cout<<"JetPt:"<<(*Jets)[i].Pt()<<" JetEta:"<<(*Jets)[i].Eta()<<" JetPhi:"<<(*Jets)[i].Phi()<<endl;
  }
  cout<<"MET:"<<MET<<" METPhi:"<<METPhi<<" MHTPhi:"<<MHTPhi<<" DPhi1:"<<DeltaPhi1<<" DeltaPhi2:"<<DeltaPhi2<<" DeltaPhi3:"<<DeltaPhi3<<" DeltaPhi4:"<<DeltaPhi4<<endl;
  for(int i=0;i<JetsAK8->size();i++){
    cout<<"AK8 pT:"<<(*JetsAK8)[i].Pt()<<" eta:"<<(*JetsAK8)[i].Eta()<<" phi:"<<(*JetsAK8)[i].Phi()<<" softDropM:"<<(*JetsAK8_softDropMass)[i]<<" tau21:"<<(*JetsAK8_NsubjettinessTau2)[i]/(*JetsAK8_NsubjettinessTau1)[i]<<" DeepDoubleBDiscr:"<<(*JetsAK8_deepDoubleBDiscriminatorH)[i]<<" WDiscr:"<<(*JetsAK8_wDiscriminatorDeep)[i]<<" WdiscrDecorr:"<<(*JetsAK8_wDiscriminatorDeepDecorrel)[i]<<endl;
  }
}
