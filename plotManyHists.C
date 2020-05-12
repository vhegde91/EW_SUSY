#include<iostream>
#include<iomanip>
#include"TH1.h"
#include"TROOT.h"
#include"TH2.h"
#include"TFile.h"
#include"TDirectory.h"
#include"TF1.h"
#include<string>
#include<vector>
#include"TLegend.h"
#include"TLatex.h"
#include"TCanvas.h"
#include"THStack.h"
#include"TStyle.h"

const int nfiles=2;
TFile *f[nfiles];
TH1D *h1[nfiles];
TLegend *lg1;
int col[10]={kRed,kTeal+9,kBlue,kGray+1,kCyan,kMagenta+2,kYellow+2,kRed,kMagenta,kOrange-9};
vector<TString> fNames,histName,legName;
vector<int> rebin;
TString xAxisName="MET (GeV)";
TString yAxisName="Events";
TString ratioName = "#frac{Fail}{Pass}";
TLatex textOnTop,intLumiE;
double intLumi=137;
char name2[100];

void makeRatioPlot();
void plotManyHists(){
  //-------------------------------------
  gStyle->SetOptStat(0);

  //-------------------------------------
  fNames.push_back("VJets_MC2018_2ndTagF.root"); 
  histName.push_back("WZW/MET_Wmc_WWZ");
  legName.push_back("W+Z+t: 2ndAK8 fail W-corr");
  rebin.push_back(5);
  //-------------------------------------
  fNames.push_back("VJets_MC2018_2ndTagP.root"); 
  histName.push_back("WZW/MET_Wmc_WWZ");
  legName.push_back("W+Z: 2ndAK8 pass W-corr");
  rebin.push_back(5);
  //-------------------------------------
  // fNames.push_back("");
  // histName.push_back("");
  // legName.push_back("");
  // rebin.push_back(1);
  //-------------------------------------
  if(fNames.size()!=nfiles){cout<<"check no. of files and filenames"<<endl; return;}

  TCanvas *c1=new TCanvas("canvasA",yAxisName+" vs "+xAxisName,1500,850);
  //  lg1=new TLegend(0.65, 0.88,  0.83, 0.65);
  lg1=new TLegend(0.65, 0.88,  0.83, 0.75);
  for(int i=0;i<nfiles;i++){
    f[i] = new TFile(fNames[i]);
    f[i]->GetObject(histName[i],h1[i]);
    h1[i]->Rebin(rebin[i]);
    //    h1[i]->Scale(1.0/h1[i]->Integral());
    h1[i]->SetTitle(0);
    h1[i]->SetLineWidth(2);
    h1[i]->SetLineColor(col[i]);
    lg1->AddEntry(h1[i],legName[i],"lep");
    if(i==0){
      h1[i]->GetXaxis()->SetTitle(xAxisName);
      h1[i]->GetXaxis()->SetTitleSize(0.05);
      h1[i]->GetXaxis()->SetLabelSize(0.05);
      h1[i]->GetYaxis()->SetTitle(yAxisName);
      h1[i]->GetYaxis()->SetTitleSize(0.05);
      h1[i]->GetYaxis()->SetLabelSize(0.05);
      h1[i]->Draw();
    }
    else h1[i]->Draw("sames");
  }
  
  lg1->Draw();
  textOnTop.SetTextSize(0.04);
  intLumiE.SetTextSize(0.04);
  textOnTop.DrawLatexNDC(0.12,0.91,"CMS #it{#bf{Simulation Preliminary}}");
  sprintf(name2,"#bf{%0.1f fb^{-1} (13 TeV)}",intLumi);
  intLumiE.DrawLatexNDC(0.7,0.91,name2);

  if(nfiles==2) makeRatioPlot();
  // c1->SaveAs("a.png");
}

void makeRatioPlot(){
  TString name = "c2";
  TCanvas *c2=new TCanvas(name,name,1500,800);//c_cA->Divide(4,2);
  TPad *p_top=new TPad(name+"_top",name+"_top",0,0.4,1,1);
  TPad *p_bot=new TPad(name+"_bot",name+"_bot",0,0.0,1,0.4);
  p_top->Draw();p_top->SetGridx();p_top->SetGridy();p_top->SetLogy();
  p_top->SetBottomMargin(0);
  p_bot->SetTopMargin(0);
  p_bot->SetBottomMargin(0.3);
  p_bot->Draw();p_bot->SetGridx();p_bot->SetGridy();  

  TLegend *lg=new TLegend(0.65, 0.88,  0.85, 0.65);
  // p_top->SetBottomMargin(0);
  // p_bot->SetTopMargin(0);
  // p_bot->SetBottomMargin(0.3);
  p_bot->Draw();p_bot->SetGridx();p_bot->SetGridy();

  p_top->cd();
  p_top->SetTickx();p_top->SetTicky();

  TH1D *h_rCan = (TH1D*)h1[0]->Clone("forRatioCanvas");
  
  h_rCan->GetXaxis()->SetTitle(xAxisName);
  h_rCan->GetXaxis()->SetTitleSize(0.07);
  h_rCan->GetXaxis()->SetLabelSize(0.07);
  h_rCan->GetYaxis()->SetTitle(yAxisName);
  h_rCan->GetYaxis()->SetTitleSize(0.08);
  h_rCan->GetYaxis()->SetTitleOffset(0.51);
  h_rCan->GetYaxis()->SetLabelSize(0.08);
  h_rCan->Draw("histe");
  h1[1]->Draw("sames histe");
  lg->AddEntry(h_rCan,legName[0],"lep");
  lg->AddEntry(h1[1],legName[1],"lep");
  lg->Draw();

  p_bot->cd();
  TH1D *h_numr=(TH1D*)h_rCan->Clone("h1cp");
  h_numr->Divide(h1[1]);
  h_numr->SetLineWidth(2);
  h_numr->SetFillStyle(0);
  h_numr->SetLineColor(col[0]);
  h_numr->SetMarkerColor(col[0]);
  h_numr->SetMarkerStyle(21);
  //  h_numr->SetMarkerSize();
  h_numr->SetTitle(0);
  h_numr->GetXaxis()->SetTitle(xAxisName);
  h_numr->GetXaxis()->SetTitleOffset(0.96);
  h_numr->GetXaxis()->SetTitleSize(0.13);
  h_numr->GetXaxis()->SetLabelSize(0.14);

  //name = "#frac{"+legName[0]+"}{"+legName[1]+"}";
  h_numr->GetYaxis()->SetTitle(ratioName);
  h_numr->GetYaxis()->SetTitleOffset(0.35);
  h_numr->GetYaxis()->SetTitleSize(0.13);
  h_numr->GetYaxis()->SetLabelSize(0.14);
  h_numr->GetYaxis()->SetNdivisions(505);
  h_numr->SetMaximum(1.99);
  h_numr->SetMinimum(0.);
  p_bot->cd();
  p_bot->SetTickx();p_bot->SetTicky();
    //    c_cB->cd(i+1);    p_bot[i]->cd();
    //    h_numr->Draw("e1 histe");
  h_numr->Draw("hist text");
  h_numr->Draw();

  p_top->cd();
  TLatex textOnTop,intLumiE;
  textOnTop.SetTextSize(0.05);
  intLumiE.SetTextSize(0.05);
  textOnTop.DrawLatexNDC(0.12,0.91,"CMS #it{#bf{Simulation Preliminary}}");
  sprintf(name2,"#bf{%0.1f fb^{-1} (13 TeV)}",intLumi);
  intLumiE.DrawLatexNDC(0.75,0.91,name2);
}
