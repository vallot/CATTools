void truePileUp(){

  vector<string> file;
  file.push_back("MyDataPileupHistogram.root");

  cout << "adding... " << file[0] << endl;
  TFile* f = new TFile(file[0].c_str());
  TH1 * hData = (TH1*)  f->Get("pileup") ;

  //if it is more than one file
  for(int i=1; i < file.size(); i++){
   cout << "adding... " << file[i] << endl;
   TFile* ftmp = new TFile(file[i].c_str());
   TH1 * tmp = (TH1*) ftmp->Get("pileup") ;
   hData->Add(tmp);
  }
   
  int bins = hData->GetNbinsX();
  cout << "number of bins from data= " << bins << endl;
  for(int i=1; i <= 100 ; i++){
    if( hData->GetBinContent(i) ) cout << "    " << hData->GetBinContent(i) << "," << endl;
    else cout << "    " << "0.0" << "," << endl;
  } 

  TCanvas * c = new TCanvas("c","c",1); 
  hData->Draw();
  c->Print("pileup.pdf");
}
