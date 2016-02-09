#include "CATTools/CatAnalyzer/interface/BTagWeightEvaluator.h"
#include "FWCore/ParameterSet/interface/FileInPath.h"

using namespace cat;
using namespace std;

BTagWeightEvaluator::BTagWeightEvaluator(const edm::ParameterSet& pset)
{
  method_ = pset.getParameter<string>("method");
  const string btagName = pset.getParameter<string>("btagName"); // csvv2 or ...
  if ( btagName == "csvv2" ) btagAlgo_ = BTAG_CSVv2;
  else if ( btagName == "cmvav2" ) btagAlgo_ = BTAG_cMVAv2;
  else if ( btagName == "jp" ) btagAlgo_ = BTAG_JP;
  else btagAlgo_ = "undefined"; // FIXME: Eventually raise error somewhere?

  if ( method_ == "CSVWeight" ) {
    const string rootFileNameHF = "csv_rwt_fit_hf_2016_01_28.root";
    const string rootFileNameLF = "csv_rwt_fit_lf_2016_01_28.root";
    const auto inputFileHF = edm::FileInPath("CATTools/CatAnalyzer/data/scaleFactors/"+rootFileNameHF).fullPath();
    const auto inputFileLF = edm::FileInPath("CATTools/CatAnalyzer/data/scaleFactors/"+rootFileNameLF).fullPath();

    TFile f_CSVwgt_HF(inputFileHF.c_str());
    TFile f_CSVwgt_LF(inputFileLF.c_str());

    fillCSVHistos(&f_CSVwgt_HF, &f_CSVwgt_LF, 6);

    return;
  }
  else if ( method_ == "iterativefit" ) {
    //const string csvFileName = "ttH_BTV_CSVv2_13TeV_2015D_20151120.csv";
    const string csvFileName = pset.getParameter<string>("csvFileName"); // CSVv2_prelim.csv or cMVAv2_prelim.csv
    const auto csvFile = edm::FileInPath("CATTools/CatAnalyzer/data/scaleFactors/"+csvFileName).fullPath();
    BTagCalibration calib(btagName, csvFile);
    BTagEntry::OperatingPoint op = BTagEntry::OP_RESHAPING;
    uncNames_ = {
      "central", "up_jes", "down_jes",
      "up_lf", "down_lf", "up_hf", "down_hf", 
      "up_hfstats1", "down_hfstats1", "up_hfstats2", "down_hfstats2",
      "up_lfstats1", "down_lfstats1", "up_lfstats2", "down_lfstats2",
      "up_cferr1", "down_cferr1", "up_cferr2", "down_cferr2"
    };
    for ( unsigned int i=0; i<uncNames_.size(); ++i ) {
      readers_[i] = BTagCalibrationReader(&calib, op, method_, uncNames_[i]); 
    }
  }
  else {
    const string csvFileName = pset.getParameter<string>("csvFileName"); // CSVv2_prelim.csv or cMVAv2_prelim.csv
    const auto csvFile = edm::FileInPath("CATTools/CatAnalyzer/data/scaleFactors/"+csvFileName).fullPath();
    BTagCalibration calib(btagName, csvFile);

    const string opName = pset.getParameter<string>("operatingPoint");
    BTagEntry::OperatingPoint op = BTagEntry::OP_LOOSE;
    if      ( opName == "loose"  ) op = BTagEntry::OP_LOOSE;
    else if ( opName == "medium" ) op = BTagEntry::OP_MEDIUM;
    else if ( opName == "tight"  ) op = BTagEntry::OP_TIGHT;
    uncNames_ = {"central", "up", "down"};

    readers_[0] = BTagCalibrationReader(&calib, op, method_, "central");
    readers_[1] = BTagCalibrationReader(&calib, op, method_, "up"     );
    readers_[2] = BTagCalibrationReader(&calib, op, method_, "down"   );
  }
}

double BTagWeightEvaluator::computeWeight(const cat::JetCollection& jets, const int unc) const
{
  // Do the simplest weighting methods
  if ( method_ == "CSVWeight" ) {
    double weight = 1.0;
    for ( auto& jet : jets ) weight *= computeCSVWeightFromROOT(jet, unc);
    return weight;
  }
  else if ( method_ == "iterativefit" ) {
    double weight = 1.0;
    for ( auto& jet : jets ) weight *= computeCSVWeightFromCSV(jet, unc);
    return weight;
  }

  // Other standard weights
  //auto& reader = readers_[unc];

  return 1.0;

}

double BTagWeightEvaluator::computeCSVWeightFromCSV(const cat::Jet& jet, const int unc) const
{
  const double pt = std::min(jet.pt(), 999.);
  const double aeta = std::abs(jet.eta());
  if ( pt <= 20 or aeta >= 2.4 ) return 1;

  double csv = jet.bDiscriminator(BTAG_CSVv2);
  if ( csv < 0.0 ) csv = -0.05;
  else if ( csv > 1.0 ) csv = 1.0;

  const int flav = std::abs(jet.hadronFlavour());
  BTagEntry::JetFlavor jf = BTagEntry::FLAV_UDSG;
  if ( flav == 5 ) jf = BTagEntry::FLAV_B;
  else if ( flav == 4 ) jf = BTagEntry::FLAV_C;

  int uncKey = unc;
  // Special care for the flavour dependent SFs
  if ( flav == 5 ) {
    if ( unc != LF_UP and unc != LF_DN and
        unc != HFSTAT1_UP and unc != HFSTAT1_DN and
        unc != HFSTAT2_UP and unc != HFSTAT2_DN ) uncKey = CENTRAL;
  }
  else if ( flav == 4 ) {
    if ( unc != CFERR1_UP and unc != CFERR1_DN and
        unc != CFERR2_UP and unc != CFERR2_DN ) uncKey = CENTRAL;
  }
  else {
    if ( unc != HF_UP and unc != HF_DN and
        unc != LFSTAT1_UP and unc != LFSTAT1_DN and
        unc != LFSTAT2_UP and unc != LFSTAT2_DN ) uncKey = CENTRAL;
  }

  const auto reader = readers_.find(uncKey);
  if ( reader == readers_.end() ) return 1;

  return reader->second.eval(jf, aeta, pt, csv);
}

double BTagWeightEvaluator::computeCSVWeightFromROOT(const cat::Jet& jet, const int unc) const
{
  // original code : https://github.com/cms-ttH/MiniAOD/blob/master/MiniAODHelper/interface/CSVHelper.h
  const double pt = std::min(jet.pt(), 999.);
  const double aeta = std::abs(jet.eta());
  if ( pt <= 20 or aeta >= 2.4 ) return 1;

  double csv = jet.bDiscriminator(BTAG_CSVv2);
  if ( csv < 0.0 ) csv = -0.05;
  else if ( csv > 1.0 ) csv = 1.0;

  const int flav = std::abs(jet.hadronFlavour());

  int iSysHF = 0;
  switch (unc) {
    case JES_UP: iSysHF = 1; break; // JESUp
    case JES_DN: iSysHF = 2; break; // JESDown
    case LF_UP: iSysHF = 3; break; // LFUp
    case LF_DN: iSysHF = 4; break; // LFDown
    case HFSTAT1_UP: iSysHF = 5; break; // Stats1Up
    case HFSTAT1_DN: iSysHF = 6; break; // Stats1Down
    case HFSTAT2_UP: iSysHF = 7; break; // Stats2Up
    case HFSTAT2_DN: iSysHF = 8; break; // Stats2Down
    default: iSysHF = 0;  break; // NoSys
  }

  int iSysC = 0;
  switch (unc) {
    case CFERR1_UP: iSysC = 1; break;
    case CFERR1_DN: iSysC = 2; break;
    case CFERR2_UP: iSysC = 3; break;
    case CFERR2_DN: iSysC = 4; break;
    default:  iSysC = 0; break;
  }

  int iSysLF = 0;
  switch (unc) {
    case JES_UP: iSysLF = 1; break; // JESUp
    case JES_DN: iSysLF = 2; break; // JESDown
    case HF_UP: iSysLF = 3; break; // HFUp
    case HF_DN: iSysLF = 4; break; // HFDown
    case LFSTAT1_UP: iSysLF = 5; break; // Stats1Up
    case LFSTAT1_DN: iSysLF = 6; break; // Stats1Down
    case LFSTAT2_UP: iSysLF = 7; break; // Stats2Up
    case LFSTAT2_DN: iSysLF = 8; break; // Stats2Down
    default: iSysLF = 0; break; // NoSys
  }

  int iPt = -1, iEta = -1;
  if (pt >=19.99 && pt<30) iPt = 0;
  else if (pt >=30 && pt<40) iPt = 1;
  else if (pt >=40 && pt<60) iPt = 2;
  else if (pt >=60 && pt<100) iPt = 3;
  else if (pt >=100)          iPt = 4;

  if (aeta >= 0 && aeta < 0.8)          iEta = 0;
  else if (aeta >= 0.8 && aeta < 1.6)   iEta = 1;
  else if (aeta >= 1.6 && aeta < 2.41)  iEta = 2;
  else  iEta = 2;

  if (iPt < 0 || iEta < 0) return (float) 1.0;

  if (flav == 5) {
    if(iPt>=nHFptBins_) iPt=nHFptBins_-1;
    int useCSVBin = (csv >= 0.) ? h_csv_wgt_hf[iSysHF][iPt]->FindBin(csv) : 1;
    return (float) h_csv_wgt_hf[iSysHF][iPt]->GetBinContent(useCSVBin);
  } else if (flav == 4) {
    if(iPt>=nHFptBins_) iPt=nHFptBins_-1;
    int useCSVBin = (csv >= 0.) ? c_csv_wgt_hf[iSysC][iPt]->FindBin(csv) : 1;
    return (float) c_csv_wgt_hf[iSysC][iPt]->GetBinContent(useCSVBin);
  } else {
    if (iPt >= 3)  iPt = 3; /// [30-40], [40-60] and [60-10000] only 3 Pt bins for lf
    int useCSVBin = (csv >= 0.) ? h_csv_wgt_lf[iSysLF][iPt][iEta]->FindBin(csv) : 1;
    return  (float) h_csv_wgt_lf[iSysLF][iPt][iEta]->GetBinContent(useCSVBin);
  }
}

// fill the histograms (done once)
void BTagWeightEvaluator::fillCSVHistos(TFile *fileHF, TFile *fileLF, int nHFptBins)
{
  nHFptBins_ = nHFptBins;
  for (int iSys = 0; iSys < 9; iSys++) {
    for (int iPt = 0; iPt < 5; iPt++) h_csv_wgt_hf[iSys][iPt] = 0;
    for (int iPt = 0; iPt < 3; iPt++) {
      for (int iEta = 0; iEta < 3; iEta++) h_csv_wgt_lf[iSys][iPt][iEta] = 0;
    }
  }
  for (int iSys = 0; iSys < 5; iSys++) {
    for (int iPt = 0; iPt < 5; iPt++) c_csv_wgt_hf[iSys][iPt] = 0;
  }

  // CSV reweighting /// only care about the nominal ones
  for (int iSys = 0; iSys < 9; iSys++) {
    TString syst_csv_suffix_hf = "final";
    TString syst_csv_suffix_c  = "final";
    TString syst_csv_suffix_lf = "final";

    switch (iSys) {
      case 0:
        // this is the nominal case
        break;
      case 1:
        // JESUp
        syst_csv_suffix_hf = "final_JESUp";
        syst_csv_suffix_lf = "final_JESUp";
        syst_csv_suffix_c = "final_cErr1Up";
        break;
      case 2:
        // JESDown
        syst_csv_suffix_hf = "final_JESDown";
        syst_csv_suffix_lf = "final_JESDown";
        syst_csv_suffix_c = "final_cErr1Down";
        break;
      case 3:
        // purity up
        syst_csv_suffix_hf = "final_LFUp";
        syst_csv_suffix_lf = "final_HFUp";
        syst_csv_suffix_c = "final_cErr2Up";
        break;
      case 4:
        // purity down
        syst_csv_suffix_hf = "final_LFDown";
        syst_csv_suffix_lf = "final_HFDown";
        syst_csv_suffix_c = "final_cErr2Down";
        break;
      case 5:
        // stats1 up
        syst_csv_suffix_hf = "final_Stats1Up";
        syst_csv_suffix_lf = "final_Stats1Up";
        break;
      case 6:
        // stats1 down
        syst_csv_suffix_hf = "final_Stats1Down";
        syst_csv_suffix_lf = "final_Stats1Down";
        break;
      case 7:
        // stats2 up
        syst_csv_suffix_hf = "final_Stats2Up";
        syst_csv_suffix_lf = "final_Stats2Up";
        break;
      case 8:
        // stats2 down
        syst_csv_suffix_hf = "final_Stats2Down";
        syst_csv_suffix_lf = "final_Stats2Down";
        break;
    }

    for (int iPt = 0; iPt < nHFptBins_; iPt++) {
      TH1D* h = (TH1D *)fileHF->Get(Form("csv_ratio_Pt%i_Eta0_%s", iPt, syst_csv_suffix_hf.Data()));
      if ( !h ) continue;
      h->SetDirectory(0);
      h_csv_wgt_hf[iSys][iPt] = h;
    }

    if (iSys < 5) {
      for (int iPt = 0; iPt < nHFptBins_; iPt++) {
        TH1D* h = (TH1D *)fileHF->Get(Form("c_csv_ratio_Pt%i_Eta0_%s", iPt, syst_csv_suffix_c.Data()));
        if ( !h ) continue;
        h->SetDirectory(0);
        c_csv_wgt_hf[iSys][iPt] = h;
      }
    }

    for (int iPt = 0; iPt < 4; iPt++) {
      for (int iEta = 0; iEta < 3; iEta++) {
        TH1D* h = (TH1D *)fileLF->Get(Form("csv_ratio_Pt%i_Eta%i_%s", iPt, iEta, syst_csv_suffix_lf.Data()));
        if ( !h ) continue;
        h->SetDirectory(0);
        h_csv_wgt_lf[iSys][iPt][iEta] = h;
      }
    }
  }
}

