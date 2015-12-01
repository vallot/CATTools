#ifndef CSVHelper_h__
#define CSVHelper_h__

// original code : https://github.com/cms-ttH/MiniAOD/blob/master/MiniAODHelper/interface/CSVHelper.h

#include "TFile.h"
#include "TH1D.h"
#include "CATTools/DataFormats/interface/Jet.h"

class CSVHelper
{
  public:
    CSVHelper(std::string hf="", std::string lf="");

    //double getCSVWeight(std::vector<double> jetPts, std::vector<double> jetEtas, std::vector<double> jetCSVs,
    //float getCSVWeight(std::vector<double> jetPts, std::vector<double> jetEtas, std::vector<double> jetCSVs,
    //                   std::vector<int> jetFlavors, int iSys, double &csvWgtHF, double &csvWgtLF, double &csvWgtCF);
    float getCSVWeight(cat::JetCollection jets, int iSys);//, double &csvWgtHF, double &csvWgtLF, double &csvWgtCF);


  private:
    void fillCSVHistos(TFile *fileHF, TFile *fileLF);

    // CSV reweighting
    TH1D *h_csv_wgt_hf[9][6];
    TH1D *c_csv_wgt_hf[9][6];
    TH1D *h_csv_wgt_lf[9][4][3];
};

#endif
