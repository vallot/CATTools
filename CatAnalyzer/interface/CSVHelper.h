#ifndef CSVHelper_h__
#define CSVHelper_h__

// original code : https://github.com/cms-ttH/MiniAOD/blob/master/MiniAODHelper/interface/CSVHelper.h

#include "TFile.h"
#include "TH1D.h"
#include "CATTools/DataFormats/interface/Jet.h"

class CSVHelper
{
  public:
    // nHFptBins specifies how many of these pt bins are used:
    // (jetPt >= 19.99 && jetPt < 30), (jetPt >= 30 && jetPt < 40), (jetPt >= 40 && jetPt < 60), 
    // (jetPt >= 60 && jetPt < 100), (jetPt >= 100 && jetPt < 160), (jetPt >= 160 && jetPt < 10000).
    // If nHFptBins < 6, the last on is inclusive (eg jetPt >=100 && jetPt < 10000 for nHFptBins=5).
    // The SFs from data have 5 bins, the pseudo data scale factors 6 bins.
    CSVHelper(std::string hf="", std::string lf="", int nHFptBins=6);

    //double getCSVWeight(std::vector<double> jetPts, std::vector<double> jetEtas, std::vector<double> jetCSVs,
    //float getCSVWeight(std::vector<double> jetPts, std::vector<double> jetEtas, std::vector<double> jetCSVs,
    //                   std::vector<int> jetFlavors, int iSys, double &csvWgtHF, double &csvWgtLF, double &csvWgtCF);
    double getCSVWeight(cat::JetCollection jets, int iSys);//, double &csvWgtHF, double &csvWgtLF, double &csvWgtCF);


  private:
    void fillCSVHistos(TFile *fileHF, TFile *fileLF);

    // CSV reweighting
    TH1D *h_csv_wgt_hf[9][6];
    TH1D *c_csv_wgt_hf[9][6];
    TH1D *h_csv_wgt_lf[9][4][3];
    const int nHFptBins;
};

#endif
