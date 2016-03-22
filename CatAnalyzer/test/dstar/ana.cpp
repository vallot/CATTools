#include<iostream>
#include<TH1F.h>
#include"OptionParser.h"
#include<TFile.h>
#include"dstarAna.h"
#include<unordered_map>
#include<TLorentzVector.h>



using namespace std;
using optparse::OptionParser;


class DStarAnalyzer {

private :
  unordered_map<std::string, TH1F*> hist;
  std::string inFileURL, outFileURL;

public :
  DStarAnalyzer( std::string inFileURL, std::string outFileURL){
    this->inFileURL=inFileURL;
    this->outFileURL=outFileURL;

    float d0_xbin, d0_xmin, d0_xmax;
    d0_xbin = 100; 
    d0_xmin = 1.6;
    d0_xmax = 2.2;
    auto title = TString::Format("D0 invariant mass ; D0 mass[GeV/c^2] ; Entries/%.1fMeV",(d0_xmax-d0_xmin)/d0_xbin*1000);
    bookHist1D("d0_mass",title.Data(), d0_xbin, d0_xmin, d0_xmax);
    bookHist1D("d0_mass_true",("True "+title).Data(),d0_xbin, d0_xmin, d0_xmax);
    bookHist1D("dstar_mass","D* invariant mass",100,1.6,2.2);
    bookHist1D("dstar_mass_true","True D* invariant mass",100,1.6,2.2);
    bookHist1D("dstar_mass_diff","Mass difference of D* and D0 meson",35,0.135,0.17);
    bookHist1D("dstar_mass_diff_1","Mass difference of D* and D0(p+ k-) meson",100,0.135,0.17);
    bookHist1D("dstar_mass_diff_2","Mass difference of D* and D0(p- k+) meson",100,0.135,0.17);
    bookHist1D("dstar_mass_diff_true","Mass difference of true D* and D0 meson",100,0.135,0.17);
    ana();
  }


  void bookHist1D( std::string name, std::string title, int nbin, double xmin, double xmax) {
    hist[name] = new TH1F(name.c_str(), title.c_str(), nbin, xmin, xmax);
    //hist[name]->Sumw2();
    hist[name]->SetDirectory(0);
  }

  void ana(){
    TFile* fIn = TFile::Open(inFileURL.c_str());

    TTree* tree = (TTree*)fIn->Get("cattree/nom");
    int totalEntries = tree->GetEntriesFast();
   

    dstarAna* ntuple = new dstarAna(tree);

    std::cout<<totalEntries<<std::endl;
    for( int idx = 0 ; idx< totalEntries ; ++idx) {
      ntuple->GetEntry(idx);
      //if ( !ntuple->step5 ) continue;
      for ( int i=0 ; i < ntuple->d0_pt->size() ; ++i) {
        float d0_mass = (*(ntuple->d0_m))[i];
        hist["d0_mass"]->Fill( d0_mass);
        if ( (*(ntuple->d0_true))[i]) hist["d0_mass_true"]->Fill( d0_mass); 
      }
      for( int i =0 ; i< ntuple->dstar_pt->size(); ++i) {
        TLorentzVector dau1, dau2, dau3;
        
        dau1.SetPtEtaPhiM( (*(ntuple->dstar_dau1_pt))[i], (*(ntuple->dstar_dau1_eta))[i], (*(ntuple->dstar_dau1_phi))[i], (*(ntuple->dstar_dau1_m))[i]   );
        dau2.SetPtEtaPhiM( (*(ntuple->dstar_dau2_pt))[i], (*(ntuple->dstar_dau2_eta))[i], (*(ntuple->dstar_dau2_phi))[i], (*(ntuple->dstar_dau2_m))[i]   );
        dau3.SetPtEtaPhiM( (*(ntuple->dstar_dau3_pt))[i], (*(ntuple->dstar_dau3_eta))[i], (*(ntuple->dstar_dau3_phi))[i], (*(ntuple->dstar_dau3_m))[i]   );

        float diff = (dau1+dau2+dau3).M()-(dau1+dau2).M();
        hist["dstar_mass"]->Fill( (*(ntuple->dstar_m))[i] );
        
        hist["dstar_mass_diff"]->Fill(diff);
        if ( (*(ntuple->dstar_dau1_q))[i] >0 && (*(ntuple->dstar_dau2_q))[i] <0 ) hist["dstar_mass_diff_1"]->Fill(diff);
        if ( (*(ntuple->dstar_dau1_q))[i] <0 && (*(ntuple->dstar_dau2_q))[i] >0 ) hist["dstar_mass_diff_2"]->Fill(diff);
        if ( (*(ntuple->dstar_true))[i] ) hist["dstar_mass_diff_true"]->Fill(diff);
      }
    }
    fIn->Close();
    TFile* fOut = TFile::Open(outFileURL.c_str(), "RECREATE");
    for( auto ptr = hist.begin() ; ptr != hist.end() ; ++ptr) hist[ptr->first]->Write();
    fOut->Close();

  }
};


int main()
{
  string inFileURL = "cattree.root";
  string outFileURL = "output.root";
 
  DStarAnalyzer da(inFileURL, outFileURL);
 
   

  return 0;
}
