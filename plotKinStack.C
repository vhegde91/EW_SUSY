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

char name[100];
char name2[100];
TString name3;
TLatex textOnTop,intLumiE;
const int nfiles=8,nBG=5;    //Specify no. of files
TFile *f[nfiles];
//int col[11]={kPink+1,kTeal+9,kGreen,kYellow,kOrange,kBlue,kCyan,kRed,kMagenta,kBlue+2,kPink+1};  //Specify Colors b's
vector<int> col={kTeal+9,kGreen,kYellow,kOrange,kGray+1,kBlue,kCyan,kRed,kMagenta,kBlue+2};  //Specify Colors b's

void decorate(TH1D*,int,const char*);
void decorate(THStack*,int,const char*);
void drawlegend(TH1D*,int,const char*);
void printInt(TH1D*,int,const char*);
TString getXlabel(TString);
TLegend *legend1=new TLegend(0.4501, 0.65,  0.88, 0.88);

//TLegend *legend2=new TLegend(0.7, 0.9,  0.90, 0.65);
//TLegend *legend2=new TLegend(0.6, 0.90,  0.98, 0.455);
void setLastBinAsOverFlow(TH1D*);
void plotKinStack(TString,int);
void plotKinStack(TString);
void plotKinStack(TString varName){plotKinStack(varName,1);}
void plotKinStack(TString varName, int rebin){
  double sr_Integral=0,cr_Integral=0;
  TH1::SetDefaultSumw2(1);
  gStyle->SetOptStat(0);
  gStyle->SetTitle(0);
  //  TString varName = "AK8Tau21";
  TString xLabel = getXlabel(varName);

  f[1] = new TFile("ST__MC2018.root");
  f[0] = new TFile("Rare_MC2018MC2017MC2016.root");
  //  f[0] = new TFile("QCD_HT_MC2018.root");
  f[2] = new TFile("TTJets_MC2018.root");
  f[3] = new TFile("WJetsToLNu_HT_MC2018.root");
  f[4] = new TFile("ZJetsToNuNu_HT_MC2018.root");
  f[5] = new TFile("TChiWZ_800_100_MC2018.root");
  //f[6] = new TFile("TChiWZ_800_400_MC2018.root");
  f[6] = new TFile("TChiWH_800_100_MC2018.root");
  f[7] = new TFile("TChiWW_800_100_MC2018.root");
  // f[9] = new TFile("TChiWZ_1000_1_MC2018.root");

  TCanvas *c_cA=new TCanvas(varName,"plot of a kin var",1500,900);  
  gStyle->SetTextSize(2);
  THStack *hs_var=new THStack("var_Stack","MET Stacked");
  TH1D *h_total;
  // TH1D *h_MET_R[nfiles];
  // for(int i=0;i<nfiles;i++){
  //   sprintf(name,"hist_file%i",i);
  //   h_MET_R[i]=new TH1D(name,name,21,0.5,21.5);
  // }
  vector<double> Bcnt;
  double intLumi=137;
  TLatex tl1;
  for(int i=0;i<nfiles;i++){
        
    TH1D *h_MET=(TH1D*)f[i]->FindObjectAny(varName);
    h_MET->Rebin(rebin);
    //    h_MET->GetYaxis()->SetRangeUser(100.5,20000);
    //    h_MET->SetMinimum(100);
    decorate(h_MET,i,f[i]->GetName());
    
    if(i<=(nBG-1)){
      hs_var->Add(h_MET);
      if(i==0) h_total = (TH1D*)h_MET->Clone("totalHist");
      else h_total->Add(h_MET);
    }

    if(i==nBG-1) {
      c_cA->cd();
      hs_var->Draw("BAR HIST");
      hs_var->Draw("HIST");
      hs_var->SetMinimum(0.8);
      hs_var->SetMaximum(hs_var->GetMaximum()*10);
      decorate(hs_var,i,f[i]->GetName()); 
      //hs_var->GetYaxis()->SetRangeUser(100.5,20000);
      h_total->SetFillStyle(3013);
      h_total->SetFillColor(kGray+2);
      h_total->Draw("e2 same");
    }
    if(i>=nBG){ 
      c_cA->cd(); 
      h_MET->SetMarkerStyle(20);
      h_MET->SetMarkerColor(col[i]);
      h_MET->SetLineColor(col[i]);
      h_MET->SetLineWidth(3);
      h_MET->Draw("hist same");
      //      h_MET->GetYaxis()->SetRangeUser(0.5,20000);
      //      h_MET->GetYaxis()->SetRangeUser(100.5,20000);
    }
    name3 = f[i]->GetName();
    name3.ReplaceAll("_MC2018.root","");
    drawlegend(h_MET,i,name3);
    if(i==nfiles-1){ 
      hs_var->GetXaxis()->SetTitleOffset(.90);
      hs_var->GetXaxis()->SetTitle(xLabel); hs_var->GetYaxis()->SetTitle("Events");hs_var->SetTitle(0);
    }
  }
  legend1->SetNColumns(2);
  legend1->SetBorderSize(0);
  c_cA->cd(); gPad->SetLogy();legend1->Draw();
  //  gPad->RedrawAxis();
  //  hs_var->GetXaxis()->SetTitle(xLabel);
  
  textOnTop.SetTextSize(0.04);
  intLumiE.SetTextSize(0.04);
  textOnTop.DrawLatexNDC(0.12,0.91,"CMS #it{#bf{Simulation Preliminary}}");
  sprintf(name2,"#bf{%0.1f fb^{-1} (13 TeV)}",intLumi);
  intLumiE.DrawLatexNDC(0.68,0.91,name2);

  name3 = varName+".png";
  c_cA->SaveAs(name3);
  cout<<"*****************************************************"<<endl;
  cout<<"Int Lumi(inv.fb) for file1:"<<setprecision(4)<<intLumi<<endl;
    
}

void decorate(THStack *hs,int i,const char* fname){
  //  hs->SetMinimum(0.5);
  //hs->SetTitle(0);
  hs->GetXaxis()->SetLabelSize(.05);
  hs->GetYaxis()->SetLabelSize(.05);
  hs->GetXaxis()->SetTitleSize(0.05);
  hs->GetYaxis()->SetTitleSize(0.05);
  //  drawlegend(hist,i,fname);
  //  gPad->Update();
  gStyle->SetOptStat(0);
}
void decorate(TH1D* hist,int i,const char* fname){
  hist->SetLineColor(col[i]);
  if(i<nBG) {
    hist->SetFillColor(col[i]);
    hist->SetLineColor(kBlack);
    hist->SetLineWidth(1);
  }
  else hist->SetLineWidth(2);
  hist->SetTitle(0);
  hist->GetXaxis()->SetLabelSize(.06);
  hist->GetYaxis()->SetLabelSize(.06);
  //hist->SetXLabelSize(0.05);
  hist->GetXaxis()->SetTitleSize(0.06);
  // drawlegend(hist,i,fname);
  //  gPad->Update();
  setLastBinAsOverFlow(hist);
  gStyle->SetOptStat(0);
  //Hlist.Add(hist);
}

void drawlegend(TH1D *hist,int i,const char* fname){
  gStyle->SetLegendBorderSize(0);
 
  TString lName=fname;
  lName.ReplaceAll(".root","");
  
  if(lName.Contains("ST_")){lName="Single t";}
  else if(lName.Contains("DYJetsToLL")){lName="DY(l^{+}l^{-})";}
  else if(lName.Contains("WJetsToLNu")){lName="W(l#nu)+ jets";}
  else if(lName.Contains("Rare")){lName="Rare";}
  else if(lName.Contains("TTJets")){lName="t #bar{t}";}
  else if(lName.Contains("QCD")){lName="QCD";}
  else if(lName.Contains("WGJetsToLNuG")){lName="W(l#nu)+ #gamma";}
  else if(lName.Contains("ZJetsToNuNu")){lName="Z(#nu#bar{#nu})+ jets";}
  else if(lName.Contains("TTGJets")){lName="t #bar{t}+ #gamma";}
  else if(lName.Contains("GJets")){lName="#gamma +jets";}
  else if(lName.Contains("Run2016")){lName="Data";}
  else if(lName.Contains("TChiWZ_1000_100")){lName="TChiWZ_1000_100";}
  else if(lName.Contains("TChiWZ_800_100")){lName="TChiWZ_800_100";}
  else if(lName.Contains("TChiWZ_1000_1")){lName="TChiWZ_1000_1";}
  else if(lName.Contains("TChiWZ_800_1")){lName="TChiWZ_800_1";}
  else if(lName.Contains("TChiWZ_600_1")){lName="TChiWZ_600_1";}
  //  else if(lName.Contains("T5bbbbZg_1600_150")){lName="T5bbbbZG(1.6,0.15)";}
  //  else if(lName.Contains("T5bbbbZg_1600_150")){lName="#tilde{g}_{1600}#rightarrow b#bar{b}#tilde{#chi}_{1,150}^{0}";}
  else if(lName.Contains("T5bbbbZg_1600_1550")){lName="T5bbbb_ZG_1550";}

  // const char *l_name=lName.c_str();
  if(i<nBG)legend1->AddEntry(hist,lName,"f");
  else legend1->AddEntry(hist,lName,"l");
  // legend1->SetTextSize(0.04);
}


void setLastBinAsOverFlow(TH1D* h_hist){
  double lastBinCt =h_hist->GetBinContent(h_hist->GetNbinsX()),overflCt =h_hist->GetBinContent(h_hist->GetNbinsX()+1);
  double lastBinErr=h_hist->GetBinError(h_hist->GetNbinsX()),  overflErr=h_hist->GetBinError(h_hist->GetNbinsX()+1);
  
  if(lastBinCt!=0 && overflCt!=0)
    lastBinErr = (lastBinCt+overflCt)* (sqrt( ((lastBinErr/lastBinCt)*(lastBinErr/lastBinCt)) + ((overflErr/overflCt)*(overflErr/overflCt)) ) );
  
  else if(lastBinCt==0 && overflCt!=0)
    lastBinErr = overflErr;
  else if(lastBinCt!=0 && overflCt==0)
    lastBinErr = lastBinErr;
  else lastBinErr=0;

  lastBinCt = lastBinCt+overflCt;
  h_hist->SetBinContent(h_hist->GetNbinsX(),lastBinCt);
  h_hist->SetBinError(h_hist->GetNbinsX(),lastBinErr);
    
}
TString getXlabel(TString n){
  if(n=="MET") return "MET (GeV)";
  if(n=="MHT") return "MHT (GeV)";
  if(n=="mT") return "mT(MET,AK8) [GeV]";
  if(n=="mT2J") return "mT(MET,2nd AK8) [GeV]";
  if(n=="NJets") return "N_{jets}";
  if(n=="BTags") return "N_{b jets}";
  if(n=="AK8Pt") return "pT^{AK8} (GeV)";
  if(n=="AK8J2Pt") return "pT^{2nd AK8} (GeV)";
  if(n=="AK8Eta") return "#eta^{AK8}";
  if(n=="AK8J2Eta") return "#eta^{2nd AK8}";
  if(n=="AK8Tau21") return "#tau_{21}^{AK8}";
  if(n=="AK8J2Tau21") return "#tau_{21}^{2nd AK8}";
  if(n=="AK8Mass") return "M^{AK8} (GeV)";
  if(n=="AK8J2Mass") return "M^{2nd AK8} (GeV)";
  if(n=="AK8Pt") return "pT^{AK8} (GeV)";
  if(n=="DeltaPhi1") return "#Delta#Phi_{1}";
  if(n=="DeltaPhi2") return "#Delta#Phi_{2}";
  if(n=="METvBin") return "MET (GeV)";
  if(n=="mTvBin") return "mT(MET,AK8) [GeV]";
  if(n=="mT2JvBin") return "mT(MET,2nd AK8) [GeV]";
  if(n=="mTSumvBin") return "mT^{AK8} + mT^{2AK8} [GeV]";

  return n;
}
