#include <iostream>
#include <vector>
#include <string>

#include <TROOT.h>
#include <TSystem.h>
#include <TFile.h>
#include <TH1.h>

using namespace std;

void combine(const char* fNameCen, const char* fNameUp, const char* fNameDn,
             const std::vector<std::string> fNames,
             const std::vector<std::string> plotNames,
             const char* combineByCStr)
{
  enum class COMBINE { HESSE, ENVELOPE, GAUSS, } combineBy;
  const string combineByStr(combineByCStr);
  if      ( combineByStr == "hessian"  ) combineBy = COMBINE::HESSE;
  else if ( combineByStr == "envelope" ) combineBy = COMBINE::ENVELOPE;
  else if ( combineByStr == "gaussian" ) combineBy = COMBINE::GAUSS;
  else {
    cerr << "CombineBy parameter is invalid, was " << combineByStr << endl;
    return;
  }

  TFile* fcen = TFile::Open(fNameCen);
  if ( !fcen or fcen->IsZombie() ) return;
  fcen->SetBit(TFile::kDevNull);

  cout << "@@ Processing with C " << fNameCen << endl;

  TFile* fup = TFile::Open(fNameUp, "recreate");
  TFile* fdn = TFile::Open(fNameDn, "recreate");

  // Load all files for the uncertainty combinations
  std::vector<TFile*> files;
  for ( auto fName : fNames ) {
    TFile* f = TFile::Open(fName.c_str());
    f->SetBit(TFile::kDevNull);
    files.push_back(f);
  }

  for ( auto pName : plotNames ) {
    TH1* hcen = (TH1*)fcen->Get(pName.c_str());
    if ( !hcen ) continue;
    const int nbins = hcen->GetNbinsX();
    //cout << "->" << pName << "...";

    std::vector<std::vector<double> > diffs(nbins+2);
    for ( int i=0, n=files.size(); i<n; ++i ) {
      TH1* h = (TH1*)files[i]->Get(pName.c_str());
      for ( int b = 0; b <= nbins +1; ++ b ) {
        diffs[b].push_back(h->GetBinContent(b)-hcen->GetBinContent(b));
      }
    }

    auto rdir = gSystem->DirName(pName.c_str());

    TDirectory* dup = fup->GetDirectory(rdir);
    if ( !dup ) { fup->mkdir(rdir); dup = fup->GetDirectory(rdir); }
    dup->cd();
    TH1* hup = (TH1*)hcen->Clone();

    TDirectory* ddn = fdn->GetDirectory(rdir);
    if ( !ddn ) { fdn->mkdir(rdir); ddn = fdn->GetDirectory(rdir); }
    ddn->cd();
    TH1* hdn = (TH1*)hcen->Clone();

    if ( combineBy == COMBINE::HESSE ) {
      for ( int b = 0; b <= nbins+1; ++b ) {
        // Combination for the Hessian set http://arxiv.org/pdf/1510.03865v1.pdf p.49, eqn.20
        double dysqr = 0;
        for ( auto dyi : diffs ) { dysqr += dyi[b]*dyi[b]; }
        hup->AddBinContent(b,  sqrt(dysqr));
        hdn->AddBinContent(b, -sqrt(dysqr));
      }
    }
    else if ( combineBy == COMBINE::ENVELOPE ) {
      for ( int b = 0; b <= nbins+1; ++b ) {
        // Combination by envelope, take the maximum/minimum
        double dymax = 0, dymin = 0;
        for ( auto dyi : diffs ) {
          dymax = max(dyi[b], dymax);
          dymin = min(dyi[b], dymin);
        }
        hup->AddBinContent(b, dymax);
        hdn->AddBinContent(b, dymin);
      }
    }
    else if ( combineBy == COMBINE::GAUSS ) {
      cout << "!!!! COMBINE BY GAUSSIAN may be incorrect !!!!" << endl;
      for ( int b = 0; b <= nbins+1; ++b ) {
        // Combination by Gaussian
        // FIME : To be verified!!!!!
        const int n = diffs.size();
        double dysqr = 0;
        for ( auto dyi : diffs ) { dysqr += dyi[b]*dyi[b]; }
        hup->AddBinContent(b,  sqrt(dysqr)/n);
        hdn->AddBinContent(b, -sqrt(dysqr)/n);
      }
    }

    dup->cd();
    hup->Write("", TObject::kOverwrite);

    ddn->cd();
    hdn->Write("", TObject::kOverwrite);
  }

  cout << "@@ Clean up opened files" << endl;
  for ( auto f : files ) {
    // Close input files. Traditional Close() is very slow due to the large number of histograms.
    // Use the trick by Phillipe : https://root.cern.ch/phpBB3/viewtopic.php?t=14450
    // This is OK for this read-only files.
    gROOT->GetListOfFiles()->Remove(f);
  }
  fcen->Close();
  fup->Close();
  fdn->Close();

  cout << "@@ Finished " << fNameCen << endl;
}

