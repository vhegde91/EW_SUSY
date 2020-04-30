void makeCutFlowSignal(){
  gStyle->SetOptStat(0);
  TCanvas *c = new TCanvas("cutflowSig","cutflowSig",900,500);
  TFile *f1 = TFile::Open("root://cmseos.fnal.gov//store/user/vhegde/RA2bTrees_V17/EW_SUSY/skims/tree_signal/tree_TChiWH_NLSP800_LSP100_MC2018_private.root");
  TFile *f = new TFile("signalCutflowsRA2Tree.root","recreate");
  TH1D *h1 =  (TH1D*)f1->Get("cutflow");
  h1->SetName("TChiWH_800_100");
  h1->SetLineColor(kRed);
  f->cd();
  h1->Write();h1->Draw();

  f1 = TFile::Open("root://cmseos.fnal.gov//store/user/vhegde/RA2bTrees_V17/EW_SUSY/skims/tree_signal/tree_TChiWH_NLSP800_LSP400_MC2018_private.root");
  h1 =  (TH1D*)f1->Get("cutflow");
  h1->SetName("TChiWH_800_400");
  h1->SetLineColor(kBlue);
  f->cd();
  h1->Write();h1->Draw("same");
 
 
  f1 = TFile::Open("root://cmseos.fnal.gov//store/user/vhegde/RA2bTrees_V17/EW_SUSY/skims/tree_signal/tree_TChiWH_NLSP800_LSP600_MC2018_private.root");
  h1 =  (TH1D*)f1->Get("cutflow");
  h1->SetName("TChiWH_800_600");
  h1->SetLineColor(kMagenta);
  f->cd();
  h1->Write();h1->Draw("same");
 
 
  f1 = TFile::Open("root://cmseos.fnal.gov//store/user/vhegde/RA2bTrees_V17/EW_SUSY/skims/tree_signal/tree_TChiWH_NLSP600_LSP400_MC2018_private.root");
  h1 =  (TH1D*)f1->Get("cutflow");h1->SetName("TChiWH_600_400");h1->SetLineColor(kOrange);
  f->cd();h1->Write();h1->Draw("same");
 
 
  f1 = TFile::Open("root://cmseos.fnal.gov//store/user/vhegde/RA2bTrees_V17/EW_SUSY/skims/tree_signal/tree_TChiWH_NLSP600_LSP100_MC2018_private.root");
  h1 =  (TH1D*)f1->Get("cutflow");h1->SetName("TChiWH_600_100");h1->SetLineColor(kTeal+9);
  f->cd();h1->Write();h1->Draw("same");
 
 
  f1 = TFile::Open("root://cmseos.fnal.gov//store/user/vhegde/RA2bTrees_V17/EW_SUSY/skims/tree_signal/tree_TChiWH_NLSP400_LSP1_MC2018_private.root");
  h1 =  (TH1D*)f1->Get("cutflow");h1->SetName("TChiWH_400_1");h1->SetLineColor(kCyan);
  f->cd();h1->Write();h1->Draw("same");
  
  c->SaveAs("cutflowSig.pdf");
}
