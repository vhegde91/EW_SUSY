void checkTree(TString name){
  auto *f = new TFile(name);
  auto *t = (TTree*)f->Get("TreeMaker2/PreSelection");
  t->Scan("SusyMotherMass:SusyLSPMass:CrossSection:NumEvents:Weight","","",1,5);
  //  cout<<endl;
}
