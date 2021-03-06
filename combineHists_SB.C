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

const int nfiles=6;    //Specify no. of files
TFile *f[nfiles];
//int col[10]={kBlack,kPink-2,kOrange,kYellow,kGreen,kTeal+9,kPink+1,kCyan,kBlue,kRed};  //Specify Colors
int col[11]={kRed,kPink+1,kTeal+9,kGreen,kYellow,kOrange,kBlue,kCyan,kMagenta,kBlue+2,kPink+1};
//int col[10]={kBlack,kBlue,kOrange,kYellow,kGreen,kCyan,kPink+1,kTeal+9,kBlue,kRed};  //Specify Colors
//int col[10]={kBlack,kPink-2,kGreen,kYellow,kBlue,kCyan,kPink+1,kCyan,kBlue,kRed};  //Specify Colors

TString name;
bool saveCanvas=1;
void setLastBinAsOverFlow(TH1D*);
void setMyRange(TH1D*,double,double);
void setMyRange(THStack*,double,double);
TString getLegName(TString);
TString getXaxisName(TString);
TH1D* sqrtHist(TH1D*);
double data_Integral=0,mc_Integral=0;
TLatex textOnTop,intLumiE;
double intLumi=137;
void combineHists_SB(){
  TH1::SetDefaultSumw2(1);
  gStyle->SetOptStat(0);
  //  gStyle->SetOptStat("nemri");
  f[0] = new TFile("TChiWH_800_100_MC2018.root");
  f[1] = new TFile("ST__MC2018.root");
  f[2] = new TFile("Rare_MC2018MC2017MC2016.root");
  //f[2] = new TFile("Rare_MC2018.root");
  //  f[3] = new TFile("QCD_HT_MC2018.root");
  f[3] = new TFile("TTJets_MC2018.root");
  f[4] = new TFile("WJetsToLNu_HT_MC2018.root");
  f[5] = new TFile("ZJetsToNuNu_HT_MC2018.root");

  //f[6] = new TFile("TChiWZ_800_100_MC2018.root");
  //  f[7] = new TFile("TChiWZ_1000_100_MC2018.root");
  // f[9] = new TFile("TChiWZ_1000_1_MC2018.root");
 
  
  vector<TString> name1;
  vector<int> rebin;
  vector<double> xLow,xHigh;
  
  //name1.push_back("WH/METvBin_1Wt1Wm1Ht1Hm");  rebin.push_back(1); xLow.push_back(-2000000); xHigh.push_back(1000000);//WH
  //------ Fully W-tagged only
  // name1.push_back("WH/METvBin_1Wt1Wm1Ht0Hm");  rebin.push_back(1); xLow.push_back(-2000000); xHigh.push_back(1000000);//W
  // name1.push_back("WH/METvBin_1Wt1Wm0Ht1Hm");  rebin.push_back(1); xLow.push_back(-2000000); xHigh.push_back(1000000);//W
  // name1.push_back("WH/METvBin_1Wt1Wm0Ht0Hm");  rebin.push_back(1); xLow.push_back(-100000); xHigh.push_back(1000000);//W

  //------- No AK8 passing W/H tag AND nass. Use AK4 pair for W and H.
  // name1.push_back("WH/METvBin_1Wt0wm1Ht0hm");  rebin.push_back(1); xLow.push_back(-100000); xHigh.push_back(1000000);
  // name1.push_back("WH/METvBin_1Wt0wm0Ht0hm");  rebin.push_back(1); xLow.push_back(-100000); xHigh.push_back(1000000);
  // name1.push_back("WH/METvBin_0Wt1wm0Ht0hm");  rebin.push_back(1); xLow.push_back(-100000); xHigh.push_back(1000000);
  // name1.push_back("WH/METvBin_0Wt0wm1Ht0hm");  rebin.push_back(1); xLow.push_back(-100000); xHigh.push_back(1000000);

  // name1.push_back("WH/METvBin_0Wt1wm1Ht0hm");  rebin.push_back(1); xLow.push_back(-100000); xHigh.push_back(1000000);
  // name1.push_back("WH/METvBin_0Wt0wm0Ht0hm");  rebin.push_back(1); xLow.push_back(-100000); xHigh.push_back(1000000);
  
  // name1.push_back("WH/METvBin_0Wt1Wm0Ht1Hm");  rebin.push_back(1); xLow.push_back(-100000); xHigh.push_back(1000000);
  // name1.push_back("WH/METvBin_1Wt0wm0Ht1hm");  rebin.push_back(1); xLow.push_back(-100000); xHigh.push_back(1000000);
  // name1.push_back("WH/METvBin_0Wt0wm0Ht1hm");  rebin.push_back(1); xLow.push_back(-100000); xHigh.push_back(1000000);

  //  deepWdiscr deepDoubleBdiscr
  //------ Fully H-tagged only
  name1.push_back("WH/METvBin_1Wt0Wm1Ht1Hm");  rebin.push_back(1); xLow.push_back(-100000); xHigh.push_back(1000000);//H
  name1.push_back("WH/METvBin_0Wt1Wm1Ht1Hm");  rebin.push_back(1); xLow.push_back(-100000); xHigh.push_back(1000000);//H
  name1.push_back("WH/METvBin_0Wt0Wm1Ht1Hm");  rebin.push_back(1); xLow.push_back(-100000); xHigh.push_back(1000000);//H

  //----- using AK4 jets only. No AK8 info is used.
  //name1.push_back("WH/METvBin_0-1Wm1Hm");  rebin.push_back(1); xLow.push_back(-100000); xHigh.push_back(1000000);//hm
  //  name1.push_back("WH/METvBin_0-1Wm0Hm");  rebin.push_back(1); xLow.push_back(-100000); xHigh.push_back(1000000);//w
  // name1.push_back("WH/METvBin_0-0Wm1Hm");  rebin.push_back(1); xLow.push_back(-100000); xHigh.push_back(1000000);//h

  TLegend *legend;//=new TLegend(0.6, 0.90,  0.98, 0.45);
  TCanvas *c_cA;
  TPad *p_top;
  TPad *p_bot;
  THStack *hs_hist;
  TH1D *h_sum;
  
  name="Combined_Hists";//+name2[i];
  c_cA=new TCanvas(name,name,1500,800);//c_cA->Divide(4,2);
  p_top=new TPad(name+"_top",name+"_top",0,0.4,1,1);
  p_bot=new TPad(name+"_bot",name+"_bot",0,0.0,1,0.4);
  p_top->Draw();p_top->SetGridx();p_top->SetGridy();p_top->SetLogy();
  p_top->SetBottomMargin(0);
  p_bot->SetTopMargin(0);
  p_bot->SetBottomMargin(0.3);
  p_bot->Draw();p_bot->SetGridx();p_bot->SetGridy();    
  name=name1[0]+"_Stack";
  hs_hist = new THStack(name,name);
  legend=new TLegend(0.55, 0.88,  0.85, 0.63);
  //legend->SetBorderSize(0);
  //cout<<getLegName(f[0]->GetName());

  TH1D *h_histG,*h_histE,*h_histi;
  name=name1[0];
  f[0]->GetObject(name,h_histG);
  h_histG->Rebin(rebin[0]);
  setLastBinAsOverFlow(h_histG);
  setMyRange(h_histG,xLow[0],xHigh[0]);
  for(int i=1;i<name1.size();i++){
    name=name1[i];
    f[0]->GetObject(name,h_histi);
    h_histi->Rebin(rebin[i]);
    setLastBinAsOverFlow(h_histi);
    setMyRange(h_histi,xLow[i],xHigh[i]);
    h_histG->Add(h_histi);
  } 

  c_cA->cd();p_top->cd();
  p_top->SetTickx();p_top->SetTicky();
  h_histG->SetLineColor(col[0]);
  h_histG->SetLineWidth(2);
  //    h_histG->SetMarkerStyle(20);
  h_histG->SetMarkerColor(h_histG->GetLineColor());
  legend->AddEntry(h_histG,getLegName(f[0]->GetName()),"l");
  //gPad->Update();
  // gPad->BuildLegend(0.75,0.75,0.95,0.95,"");
  
  for(int p=1;p<nfiles;p++){
    f[p]->GetObject(name1[0],h_histE);
    h_histE->Rebin(rebin[0]);
    setLastBinAsOverFlow(h_histE);
    setMyRange(h_histE,xLow[0],xHigh[0]);
    for(int i=1;i<name1.size();i++){
      f[p]->GetObject(name1[i],h_histi);
      h_histi->Rebin(rebin[i]);
      setLastBinAsOverFlow(h_histi);
      setMyRange(h_histi,xLow[i],xHigh[i]);
      h_histE->Add(h_histi);
    }
    h_histE->SetLineColor(kBlack);
    h_histE->SetLineWidth(1);
    //      h_histE->SetMarkerStyle(21);
    h_histE->SetMarkerColor(col[p]);
    h_histE->SetFillColor(col[p]);
    hs_hist->Add(h_histE);
    
    if(p==1){
      name=name1[0]+"_Sum";
      h_sum = (TH1D*)h_histE->Clone(name);
    }
    else h_sum->Add(h_histE);
    legend->AddEntry(h_histE,getLegName(f[p]->GetName()),"f");
  }
  c_cA->cd();    p_top->cd();
  hs_hist->SetMinimum(0.08);
  hs_hist->SetMaximum(5*hs_hist->GetMaximum());
  hs_hist->Draw("HIST");
    
  h_sum->SetFillColor(1);
  h_sum->SetFillStyle(3013);
  h_sum->Draw("e2 same");

  name = h_histG->GetName();
  if(xLow[0] > -5000. && xHigh[0] < 5000.) hs_hist->GetXaxis()->SetRangeUser(xLow[0],xHigh[0]);
  hs_hist->GetYaxis()->SetNdivisions(5);
  h_histG->Draw("histe same");
  hs_hist->SetTitle(";;Events");
  hs_hist->GetYaxis()->SetTitleOffset(0.50);    
  hs_hist->GetYaxis()->SetTitleSize(0.09);
  hs_hist->GetYaxis()->SetLabelSize(0.09);
    
  c_cA->Modified();
  c_cA->Update();
  legend->SetNColumns(3);
  legend->Draw();
  
  TH1D *h_numr=(TH1D*)h_sum->Clone();
  h_numr = sqrtHist(h_numr);
  h_numr->Scale(2);
  // for(int i=1;i<=h_histG->GetNbinsX();i++){
  //   if(h_histG->GetBinContent(i) < 1e-6)
  // 	h_histG->SetBinContent(i,1e-6);
  // }
  h_numr->Divide(h_histG);
  for(int i=1;i<=h_numr->GetNbinsX();i++){
    if(h_numr->GetBinContent(i) < 1e-10)
      h_numr->SetBinContent(i,1000);
  }
  h_numr->SetLineWidth(2);
  h_numr->SetFillStyle(0);
  h_numr->SetLineColor(col[0]);
  h_numr->SetMarkerColor(col[0]);
  h_numr->SetMarkerStyle(21);
  h_numr->SetMarkerSize(3);
  h_numr->SetTitle(0); name=name1[0];
  h_numr->GetXaxis()->SetTitle(getXaxisName(name));
  h_numr->GetXaxis()->SetTitleOffset(0.96);
  h_numr->GetXaxis()->SetTitleSize(0.13);
  h_numr->GetXaxis()->SetLabelSize(0.13);
  
  h_numr->GetYaxis()->SetTitle("#frac{2#sqrt{BG}}{Signal}");
  h_numr->GetYaxis()->SetTitleOffset(0.35);
  h_numr->GetYaxis()->SetTitleSize(0.13);
  h_numr->GetYaxis()->SetLabelSize(0.14);
  h_numr->GetYaxis()->SetNdivisions(505);
  h_numr->SetMaximum(4.6);
  h_numr->SetMinimum(0.);
  c_cA->cd();    p_bot->cd();
  p_bot->SetTickx();p_bot->SetTicky();
  //    c_cB->cd(i+1);    p_bot->cd();
  //    h_numr->Draw("e1 histe");
  h_numr->Draw("hist");
  
  c_cA->cd();    p_top->cd(); gPad->RedrawAxis();
  char name2[100];
  textOnTop.SetTextSize(0.08);
  intLumiE.SetTextSize(0.08);
  textOnTop.DrawLatexNDC(0.12,0.91,"CMS #it{#bf{Simulation Preliminary}}");
  intLumiE.SetTextSize(0.08);
  sprintf(name2,"#bf{%0.1f fb^{-1} (13 TeV)}",intLumi);
  intLumiE.DrawLatexNDC(0.7,0.91,name2);
  if(saveCanvas){name="cobminedHists.png";c_cA->SaveAs(name);}
  
  cout<<"Integral in Signal "<<data_Integral<<endl
      <<"Integral in MC "<<mc_Integral<<endl;
}

TString getLegName(TString lName){
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
  else if(lName.Contains("TChiWW_800_100")){lName="TChiWW_800_100";}
  else if(lName.Contains("TChiWH_800_100")){lName="TChiWH_800_100";}
  else if(lName.Contains("TChiWZ_1000_1")){lName="TChiWZ_1000_1";}
  else if(lName.Contains("TChiWZ_800_1")){lName="TChiWZ_800_1";}
  else if(lName.Contains("TChiWZ_600")){lName="TChiWZ(600)";}
  else if(lName.Contains("_MC2018.root")){lName.ReplaceAll("_MC2018.root","");}
  return lName;
}
TString getXaxisName(TString n){
  if(n.Contains("mTvBin")) return "mT(MET,AK8) [GeV]";
  if(n.Contains("mT2JvBin")) return "mT(MET,2nd AK8) [GeV]";
  if(n.Contains("mTSumvBin")) return "mT^{AK8} + mT^{2AK8} [GeV]";
  if(n.Contains("MHT")) return "MHT (GeV)";
  if(n.Contains("NJets")) return "N_{jets}";
  if(n.Contains("BTags")) return "N_{b jets}";
  if(n.Contains("AK8Pt")) return "pT^{AK8} (GeV)";
  if(n.Contains("AK8J2Pt")) return "pT^{2nd AK8} (GeV)";
  if(n.Contains("AK8Eta")) return "#eta^{AK8}";
  if(n.Contains("AK8J2Eta")) return "#eta^{2nd AK8}";
  if(n.Contains("AK8Tau21")) return "#tau_{21}^{AK8}";
  if(n.Contains("AK8J2Tau21")) return "#tau_{21}^{2nd AK8}";
  if(n.Contains("AK8Mass")) return "M^{AK8} (GeV)";
  if(n.Contains("AK8J2Mass")) return "M^{2nd AK8} (GeV)";
  if(n.Contains("AK8Pt")) return "pT^{AK8} (GeV)";
  if(n.Contains("DeltaPhi1")) return "#Delta#Phi_{1}";
  if(n.Contains("DeltaPhi2")) return "#Delta#Phi_{2}";
  if(n.Contains("RA2bBins")) return "Bin";
  if(n.Contains("MET")) return "MET (GeV)";
  if(n.Contains("METvBin")) return "MET (GeV)";

  if(n.Contains("WH/")){
    n.ReplaceAll("WH/","");
    n.ReplaceAll("_"," for ");
    return getXaxisName(n);
  }
  return n;
}

void setMyRange(TH1D *h1,double xLow,double xHigh){
  double err=0;
  if(xHigh > 13000) return;
  if(xLow < 13000) return;
  int nMax=h1->FindBin(xHigh);
  h1->SetBinContent(nMax,h1->IntegralAndError(nMax,h1->GetNbinsX(),err));
  h1->SetBinError(nMax,err);
  for(int i=nMax+1;i<=h1->GetNbinsX()+1;i++){
    h1->SetBinContent(i,0);
    h1->SetBinError(i,0);
  }
  //  h1->GetXaxis()->SetRangeUser(xLow,xHigh);
}
void setMyRange(THStack *h1,double xLow,double xHigh){
  // double err=0;
  // if(xHigh > 13000) return;
  // int nMax=h1->FindBin(xHigh);
  // h1->SetBinContent(nMax,h1->IntegralAndError(nMax,h1->GetNbinsX(),err));
  // h1->SetBinError(nMax,err);
  // for(int i=nMax+1;i<=h1->GetNbinsX()+1;i++){
  //   h1->SetBinContent(i,0);
  //   h1->SetBinError(i,0);
  // }
  //h1->GetXaxis()->SetRangeUser(xLow,xHigh);
}

void setLastBinAsOverFlow(TH1D* h_hist){
  double lastBinCt =h_hist->GetBinContent(h_hist->GetNbinsX()),overflCt =h_hist->GetBinContent(h_hist->GetNbinsX()+1);
  double lastBinErr=h_hist->GetBinError(h_hist->GetNbinsX()),  overflErr=h_hist->GetBinError(h_hist->GetNbinsX()+1);
  
  if(lastBinCt!=0 && overflCt!=0)
    lastBinErr = sqrt( (lastBinErr*lastBinErr) + (overflErr*overflErr) );
  
  else if(lastBinCt==0 && overflCt!=0)
    lastBinErr = overflErr;
  else if(lastBinCt!=0 && overflCt==0)
    lastBinErr = lastBinErr;
  else lastBinErr=0;

  lastBinCt = lastBinCt+overflCt;
  h_hist->SetBinContent(h_hist->GetNbinsX(),lastBinCt);
  h_hist->SetBinError(h_hist->GetNbinsX(),lastBinErr);
    
}

TH1D* sqrtHist(TH1D *h1){
  for(int i=0;i<=h1->GetNbinsX()+1;i++){
    if(h1->GetBinContent(i) >=0){
      h1->SetBinContent(i,sqrt(h1->GetBinContent(i)));
      h1->SetBinError(i,0.5*(h1->GetBinError(i)));
    }
    else{
      h1->SetBinContent(i,sqrt(-1*h1->GetBinContent(i)));
      h1->SetBinError(i,0.5*(h1->GetBinError(i)));
    }
  }
  return h1;
}
