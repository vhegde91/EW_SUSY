#ifndef Skimming_H
#define Skimming_H

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

class Skimming : public NtupleVariables{

 public:
  Skimming(const TString &inputFileList="foo.txt", const char *outFileName="histo.root",const char *dataset="data");
  ~Skimming();
  Bool_t   FillChain(TChain *chain, const TString &inputFileList);
  Long64_t LoadTree(Long64_t entry);
  void     EventLoop(const char *,const char *);
  void     BookHistogram(const char *);
  void print(Long64_t);
  double wt=0,lumiInfb=137.0;
  
  //Variables defined
  TH1D *h_filters;
  TH1F *h_cutflow;
 
  TFile *oFile;
  
};
#endif

#ifdef Skimming_cxx

void Skimming::BookHistogram(const char *outFileName) {

  //  char hname[200], htit[200];
  //  double xlow = 0.0,  xhigh = 2000.0;
  //  int nbins = 2000;
  TString name,title;
 
  oFile = new TFile(outFileName, "recreate");
  //  oFile->mkdir("TreeMaker2");
  //  oFile->cd("TreeMaker2");
  TH1::SetDefaultSumw2(1);
  h_cutflow = new TH1F("CutFlow","cut flow",15,0,15);
  h_filters = new TH1D("Filters","Filters: Bin1 : all nEvnts, other bins: filter pass/fail",10,0,10);

}

Skimming::Skimming(const TString &inputFileList, const char *outFileName, const char* dataset) {
  string nameData=dataset;
  TString nameData2 = nameData;
  TChain *tree = new TChain("tree");
  if( ! FillChain(tree, inputFileList) ) {
    std::cerr << "Cannot get the tree " << std::endl;
  } else {
    std::cout << "Initiating analysis of dataset " << dataset << std::endl;
  }

  cout<<"Treating the input files as "<<nameData<<" for setting tree branches"<<endl;
  NtupleVariables::Init(tree,nameData);

  BookHistogram(outFileName);
  
}

Bool_t Skimming::FillChain(TChain *chain, const TString &inputFileList) {
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

Long64_t Skimming::LoadTree(Long64_t entry) {
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

Skimming::~Skimming() { 

  if (!fChain) return;
  delete fChain->GetCurrentFile();
  oFile->Write();
  oFile->Close();

}

#endif

