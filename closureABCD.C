#include "headerFunctions.h"
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

TLatex textOnTop,intLumiE;
double intLumi=137;
TString name;
bool saveCanvas=0;

void closureABCD(){
  vector<TString> fNames,histName,legName;
  gStyle->SetOptStat(0);
  int rebin;
  int col[2] = {kOrange+8,kBlue};
  TString ratioName = "#frac{Pred }{Exp}";
  TString xAxisName="MET (GeV)";
  TString yAxisName="Events";
 
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

  rebin = 1;
  bool changeXrange = 0;
  double xMin = 200, xMax = 1500;
  //-------------------------------------
  fNames.push_back("VJets_MC2018_2ndTagF.root");
  //  fNames.push_back("TotalBG_2ndTagF.root");
  legName.push_back("W+Z: Prediction");
  histName.push_back("WZW/METvBin_Tau21_WWZ");
  //-------------------------------------
  fNames.push_back("VJets_MC2018_2ndTagP.root");
  //fNames.push_back("TotalBG_2ndTagP.root");
  histName.push_back(histName[0]);
  legName.push_back("W+Z: Expected");
    
  TH1D *h1[2];
    
  TFile *f[2];
  for(int i=0;i<2;i++){
    f[i] = new TFile(fNames[i]);
    f[i]->GetObject(histName[i],h1[i]);
    h1[i]->Rebin(rebin);
    h1[i]->SetLineColor(col[i]);
    if(changeXrange) setMyRange(h1[i],xMin,xMax);
    setLastBinAsOverFlow(h1[i]);
  }
  //  h1[0]->Scale(12.3/69.12);
  //  h1[0]->Scale(15.28/80.15);
  float scaleHist = h1[1]->Integral(h1[1]->FindBin(200.001),h1[1]->FindBin(299.99))/h1[0]->Integral(h1[0]->FindBin(200.001),h1[0]->FindBin(299.99));
  h1[0]->Scale(scaleHist);

  TH1D *h_rCan = (TH1D*)h1[0]->Clone("forRatioCanvas");

  h_rCan->SetTitle(0);
  h_rCan->GetXaxis()->SetTitle(xAxisName);
  h_rCan->GetXaxis()->SetTitleSize(0.07);
  h_rCan->GetXaxis()->SetLabelSize(0.07);
  h_rCan->GetYaxis()->SetTitle(yAxisName);
  h_rCan->GetYaxis()->SetTitleSize(0.08);
  h_rCan->GetYaxis()->SetTitleOffset(0.51);
  h_rCan->GetYaxis()->SetLabelSize(0.08);
  h_rCan->SetMinimum(0.05);
  h_rCan->Draw("e1");
  h1[1]->Draw("same e1");
  lg->AddEntry(h_rCan,legName[0],"lep");
  lg->AddEntry(h1[1],legName[1],"lep");
  lg->Draw();

  p_bot->cd();
  TH1D *h_numr=(TH1D*)h_rCan->Clone("h1cp");
  h_numr->Divide(h1[1]);
  h_numr->SetLineWidth(2);
  h_numr->SetFillStyle(0);
  h_numr->SetLineColor(kBlack);
  h_numr->SetMarkerColor(h_numr->GetLineColor());
  h_numr->SetMarkerStyle(21);
  //  h_numr->SetMarkerSize();
  h_numr->SetTitle(0);
  h_numr->GetXaxis()->SetTitle(xAxisName);
  h_numr->GetXaxis()->SetTitleOffset(0.96);
  h_numr->GetXaxis()->SetTitleSize(0.13);
  h_numr->GetXaxis()->SetLabelSize(0.13);

  //name = "#frac{"+legName[0]+"}{"+legName[1]+"}";
  h_numr->GetYaxis()->SetTitle(ratioName);
  h_numr->GetYaxis()->SetTitleOffset(0.35);
  h_numr->GetYaxis()->SetTitleSize(0.13);
  h_numr->GetYaxis()->SetLabelSize(0.13);
  h_numr->GetYaxis()->SetNdivisions(505);
  h_numr->SetMaximum(1.99);
  h_numr->SetMinimum(0.01);
  p_bot->cd();
  p_bot->SetTickx();p_bot->SetTicky();
  //    c_cB->cd(i+1);    p_bot[i]->cd();
  //    h_numr->Draw("e1 histe");
  //h_numr->Draw("e1");
  h_numr->DrawCopy("e1");
  h_numr->DrawCopy("e0 same");

  TLine *line;
  if(!changeXrange) line = new TLine(h_numr->GetBinLowEdge(1),1.0, h_numr->GetBinLowEdge(h_numr->GetNbinsX()+1),1);
  else line = new TLine(xMin,1.0, xMax,1.0);
  line->SetLineStyle(3);
  line->Draw();
  
  c2->cd();    p_top->cd(); gPad->RedrawAxis();
  char name2[100];
  textOnTop.SetTextSize(0.06);
  intLumiE.SetTextSize(0.06);
  textOnTop.DrawLatexNDC(0.12,0.91,"CMS #it{#bf{Simulation Preliminary}}");
  sprintf(name2,"#bf{%0.1f fb^{-1} (13 TeV)}",intLumi);
  intLumiE.DrawLatexNDC(0.7,0.91,name2);
  if(saveCanvas){name="METclosure.png";c2->SaveAs(name);}

}
