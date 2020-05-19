#define PlotLimits_cxx
#include "PlotLimits.h"
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

  PlotLimits ana(inputFileList, outFileName, data);
  cout << "dataset " << data << " " << endl;

  ana.EventLoop(data,inputFileList);

  return 0;
}

void PlotLimits::EventLoop(const char *data,const char *inputFileList) {
  if (fChain == 0) return;

  Long64_t nentries = fChain->GetEntriesFast();
  cout << "nentries " << nentries << endl;
  cout << "Analyzing dataset " << data << " " << endl;
  Long64_t nbytes = 0, nb = 0;
  int decade = 0;
  int evtSurvived=0;

  TString s_data=data,name;

  vector<TString> fNames;
  ifstream filein(inputFileList);
  if(filein.is_open()){
    string line1;
    while(getline(filein,line1))
      fNames.push_back(line1);    
  }
  else cout<<"Could not open file: "<<inputFileList<<endl;
  filein.close();

  for (Long64_t jentry=0; jentry<nentries;jentry++) {
   
    // ==============print number of events done == == == == == == == =
    double progress = 10.0 * jentry / (1.0 * nentries);
    int k = int (progress);
    // if (k > decade)
    //   cout << 10 * k << " %" <<endl;
    // decade = k;
    // cout<<"j:"<<jentry<<" fcurrent:"<<fCurrent<<endl;
    // ===============read this entry == == == == == == == == == == == 
    Long64_t ientry = LoadTree(jentry);
    if (ientry < 0) break;
    nb = fChain->GetEntry(jentry);   nbytes += nb;
    //    cout<<fCurrent<<" "<<fNames[fCurrent]<<endl;

    //    if(quantileExpected > 0) continue;
    int m1=mh, m2=roundf((mh-m1)*10000);
    double mParent=m1,mKid=round(m2),xsec=-100,xsecUnc=0;
    name = to_string(m1)+"_"+to_string(m2);
    //    cout<<name<<endl;
    if(limit > 100. || limit < 0.0000001) cout<<"Limit:"<<limit<<" (mParent,mKid):"<<mParent<<" "<<mKid<<endl;

    if(quantileExpected < 0){
      h2_mParentmKid_r->Fill(mParent,mKid,limit);
      h2_mParentmKid_r_fb->Fill(mParent,mKid,limit);

      h2_mParentmKid_rUnc->Fill(mParent,mKid,limitErr);
      //      cout<<"("<<mKid<<" , "<<xsecUnc<<")"<<endl;
      //      cout<<"("<<mKid<<" , "<<mKid<<"),"<<endl;
    }
    if(abs(quantileExpected - 0.1599999) <= 0.0001){
      h2_mParentmKid_16pc->Fill(mParent,mKid,limit);
    }
    if(abs(quantileExpected - 0.8399999) <= 0.0001){
      h2_mParentmKid_84pc->Fill(mParent,mKid,limit);
    }
    if(abs(quantileExpected - 0.5) <= 0.0001){
      h2_mParentmKid_median->Fill(mParent,mKid,limit);
      if(jentry < 90) h_expLimitR->Fill(name,limit);
      cout<<fNames[fCurrent]<<" "<<name<<" "<<limit<<endl;
    }
    if(abs(quantileExpected - 0.025) <= 0.0001){
      h2_mParentmKid_2p5pc->Fill(mParent,mKid,limit);
    }
    if(abs(quantileExpected - 0.975) <= 0.0001){
      h2_mParentmKid_97p5pc->Fill(mParent,mKid,limit);
    }
    //    cout<<GluinoMass<<" "<<NLSPMass<<endl;
  }
}
