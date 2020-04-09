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
void resetNegativeBins(TH1D*);
void setMyRange(TH1D*,double,double);
void setMyRange(THStack*,double,double);
TString getLegName(TString);
TString getXaxisName(TString);
TH1D* sqrtHist(TH1D*);
double data_Integral=0,mc_Integral=0;
TLatex textOnTop,intLumiE;
double intLumi=137;
void signalVsTotalBG_EvtCat(){
  TH1::SetDefaultSumw2(1);
  gStyle->SetOptStat(0);
  //  gStyle->SetOptStat("nemri");
  f[0] = new TFile("TChiWH_800_100_MC2018.root");
  //f[0] = new TFile("TChiWH_800_100_MC2018.root");
  f[1] = new TFile("ST__MC2018.root");
  f[2] = new TFile("Rare_MC2018MC2017MC2016.root");
  //f[2] = new TFile("Rare_MC2018.root");
  //  f[3] = new TFile("QCD_HT_MC2018.root");
  f[3] = new TFile("TTJets_MC2018.root");
  f[4] = new TFile("WJetsToLNu_HT_MC2018.root");
  f[5] = new TFile("ZJetsToNuNu_HT_MC2018.root");

  //  f[0] = new TFile("TChiWZ_800_400_MC2018.root");
  //  f[7] = new TFile("TChiWZ_1000_100_MC2018.root");
  // f[9] = new TFile("TChiWZ_1000_1_MC2018.root");
 
  
  vector<string> name1;
  vector<int> rebin;
  vector<double> xLow,xHigh;
  //  name1.push_back("mTvBin");  rebin.push_back(1); xLow.push_back(200); xHigh.push_back(1000000);
  //  name1.push_back("MT2vBin");  rebin.push_back(1); xLow.push_back(200); xHigh.push_back(1000000);
  //  name1.push_back("AK8Pt");  rebin.push_back(10); xLow.push_back(200); xHigh.push_back(1000000);
  name1.push_back("EvtTypeWH");  rebin.push_back(1); xLow.push_back(-100000000); xHigh.push_back(1600000);
  name1.push_back("EvtTypeWH_0AK8M");  rebin.push_back(1); xLow.push_back(-100000000); xHigh.push_back(1000000);
  name1.push_back("EvtTypeWH_ak84");  rebin.push_back(1); xLow.push_back(-100000000); xHigh.push_back(1000000);

  TLegend *legend[name1.size()];//=new TLegend(0.6, 0.90,  0.98, 0.45);
  TCanvas *c_cA[name1.size()];
  TPad *p_top[name1.size()];
  TPad *p_bot[name1.size()];
  THStack *hs_hist[name1.size()];
  TH1D *h_sum[name1.size()];
  // gStyle->SetTitleYSize(0.08);
  // gStyle->SetTitleYOffset(0.68);
  for(int i=0;i<name1.size();i++){
    name=name1[i];//+name2[i];
    c_cA[i]=new TCanvas(name,name,1500,800);//c_cA[i]->Divide(4,2);
    p_top[i]=new TPad(name+"_top",name+"_top",0,0.5,1,1);
    p_bot[i]=new TPad(name+"_bot",name+"_bot",0,0.0,1,0.5);
    p_top[i]->Draw();p_top[i]->SetGridx();p_top[i]->SetGridy();p_top[i]->SetLogy();
    p_top[i]->SetBottomMargin(0);
    p_bot[i]->SetTopMargin(0);
    p_bot[i]->SetBottomMargin(0.3);
    p_bot[i]->Draw();p_bot[i]->SetGridx();p_bot[i]->SetGridy();    p_bot[i]->SetLogy();
    name=name1[i]+"_Stack";
    hs_hist[i] = new THStack(name,name);
    legend[i]=new TLegend(0.15, 0.88,  0.45, 0.63);
    //    legend[i]=new TLegend(0.7, 0.88,  0.88, 0.43);
    //legend[i]->SetBorderSize(0);
  }
  //cout<<getLegName(f[0]->GetName());
  TH1D *h_histG,*h_histE,*h_histGcopy;
  for(int i=0;i<name1.size();i++){
    c_cA[i]->cd();
    name=name1[i];
    h_histG=(TH1D*)f[0]->FindObjectAny(name);
    resetNegativeBins(h_histG);
    
    c_cA[i]->cd();p_top[i]->cd();
    p_top[i]->SetTickx();p_top[i]->SetTicky();
    h_histG->SetLineColor(col[0]);
    h_histG->SetLineWidth(2);
    
    h_histG->SetMarkerColor(h_histG->GetLineColor());
    legend[i]->AddEntry(h_histG,getLegName(f[0]->GetName()),"l");

    for(int p=1;p<nfiles;p++){  
      name=name1[i];
      h_histE=(TH1D*)f[p]->FindObjectAny(name);
      resetNegativeBins(h_histE);
      h_histE->SetLineColor(kBlack);
      h_histE->SetLineWidth(1);
    
      h_histE->SetMarkerColor(col[p]);
      h_histE->SetFillColor(col[p]);
     
      hs_hist[i]->Add(h_histE);
      if(p==1){
	name=name1[i]+"_Sum";
	h_sum[i] = (TH1D*)h_histE->Clone(name);
      }
      else h_sum[i]->Add(h_histE);
      legend[i]->AddEntry(h_histE,getLegName(f[p]->GetName()),"f");
    }
    c_cA[i]->cd();    p_top[i]->cd();
    hs_hist[i]->SetMinimum(0.08);
    hs_hist[i]->SetMaximum(10*hs_hist[i]->GetMaximum());
    hs_hist[i]->Draw("HIST");

    h_sum[i]->SetFillColor(1);
    h_sum[i]->SetFillStyle(3013);
    h_sum[i]->Draw("e2 same");

    name = h_histG->GetName();
    
    hs_hist[i]->GetYaxis()->SetNdivisions(5);
    //    if(name.Contains("MuPt")) setMyRange(hs_hist[i],xLow[i],xHigh[i]);
    h_histG->Draw("histe same");
    hs_hist[i]->SetTitle(";;Events");
    hs_hist[i]->GetYaxis()->SetTitleOffset(0.50);    
    hs_hist[i]->GetYaxis()->SetTitleSize(0.09);
    hs_hist[i]->GetYaxis()->SetLabelSize(0.09);

    c_cA[i]->Modified();
    c_cA[i]->Update();
    legend[i]->SetNColumns(2);
    legend[i]->Draw();
    TLine *line1=new TLine( 0,5.0,  h_histG->GetNbinsX(),5.0);
    line1->SetLineWidth(2);
    line1->SetLineColor(kBlue);
    line1->SetLineStyle(2);
    line1->Draw();

    TH1D *h_numr=(TH1D*)h_histG->Clone();
    
    //    h_numr->Scale(10);
    // for(int i=1;i<=h_histG->GetNbinsX();i++){
    //   if(h_histG->GetBinContent(i) < 1e-6)
    // 	h_histG->SetBinContent(i,1e-6);
    // }
    h_numr->Divide(h_sum[i]);
    for(int i=1;i<=h_numr->GetNbinsX();i++){
      if(h_numr->GetBinContent(i) < 1e-10)
	h_numr->SetBinContent(i,1000);
    }
    h_numr->SetLineWidth(2);
    h_numr->SetFillStyle(0);
    h_numr->SetLineColor(col[0]);
    h_numr->SetMarkerColor(col[0]);
    h_numr->SetMarkerStyle(21);
    h_numr->SetMarkerSize(2);
    h_numr->SetTitle(0); name=name1[i];
    //    h_numr->GetXaxis()->SetTitle(getXaxisName(name));
    h_numr->GetXaxis()->SetTitleOffset(0.96);
    h_numr->GetXaxis()->SetTitleSize(0.13);
    h_numr->GetXaxis()->SetLabelSize(0.11);

    h_numr->GetYaxis()->SetTitle("#frac{Signal}{BG}");
    h_numr->GetYaxis()->SetTitleOffset(0.47);
    h_numr->GetYaxis()->SetTitleSize(0.11);
    h_numr->GetYaxis()->SetLabelSize(0.1);
    h_numr->GetYaxis()->SetNdivisions(505);
    h_numr->SetMaximum(1.8);
    h_numr->SetMinimum(0.0001);
    c_cA[i]->cd();    p_bot[i]->cd();
    p_bot[i]->SetTickx();p_bot[i]->SetTicky();
    //    c_cB->cd(i+1);    p_bot[i]->cd();
    //    h_numr->Draw("e1 histe");
    h_numr->Draw("hist text");

    c_cA[i]->cd();    p_top[i]->cd(); gPad->RedrawAxis();
    char name2[100];
    textOnTop.SetTextSize(0.08);
    intLumiE.SetTextSize(0.08);
    textOnTop.DrawLatexNDC(0.12,0.91,"CMS #it{#bf{Simulation Preliminary}}");
    intLumiE.SetTextSize(0.08);
    sprintf(name2,"#bf{%0.1f fb^{-1} (13 TeV)}",intLumi);
    intLumiE.DrawLatexNDC(0.7,0.91,name2);
    if(saveCanvas){name=name1[i]+".png";c_cA[i]->SaveAs(name);}
    
  }
  
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
  return lName;
}
TString getXaxisName(TString n){
  if(n=="MET") return "MET (GeV)";
  if(n=="METvBin") return "MET (GeV)";
  if(n=="mTvBin") return "mT(MET,AK8) [GeV]";
  if(n=="mT2JvBin") return "mT(MET,2nd AK8) [GeV]";
  if(n=="mTSumvBin") return "mT^{AK8} + mT^{2AK8} [GeV]";
  if(n=="MHT") return "MHT (GeV)";
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
  if(n=="RA2bBins") return "Bin";

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

void resetNegativeBins(TH1D *h_hist){
  for(int i=1;i<=h_hist->GetNbinsX();i++){
    if(h_hist->GetBinContent(i) < 0){
      h_hist->SetBinContent(i,0);
      h_hist->SetBinError(i,0);
    }
  }
}
