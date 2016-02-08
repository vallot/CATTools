#include "CATTools/CatAnalyzer/interface/BTagScaleFactorEvaluators.h"
#include "FWCore/ParameterSet/interface/FileInPath.h"

using namespace cat;
using namespace std;

CSVWeightEvaluator::CSVWeightEvaluator(const int inputType)
{
  const string csvFileName = "ttH_BTV_CSVv2_13TeV_2015D_20151120.csv";
  const string rootFileNameHF = "csv_rwt_fit_hf_2016_01_28.root";
  const string rootFileNameLF = "csv_rwt_fit_lf_2016_01_28.root";

  if ( inputType == CSV ) {
    const char* method = "iterativefit";
    // setup calibration readers (once)
    const auto csvFile = edm::FileInPath("CATTools/CatAnalyzer/data/scaleFactors/"+csvFileName).fullPath();
    BTagCalibration calib_csvv2("csvv2", csvFile);
    readers_[CENTRAL] = BTagCalibrationReader(&calib_csvv2, BTagEntry::OP_RESHAPING, method, "central");

    readers_[JES_UP] = BTagCalibrationReader(&calib_csvv2, BTagEntry::OP_RESHAPING, method, "up_jes");
    readers_[JES_DN] = BTagCalibrationReader(&calib_csvv2, BTagEntry::OP_RESHAPING, method, "down_jes");

    readers_[LF_UP] = BTagCalibrationReader(&calib_csvv2, BTagEntry::OP_RESHAPING, method, "up_lf");
    readers_[LF_DN] = BTagCalibrationReader(&calib_csvv2, BTagEntry::OP_RESHAPING, method, "down_lf");

    readers_[HF_UP] = BTagCalibrationReader(&calib_csvv2, BTagEntry::OP_RESHAPING, method, "up_hf");
    readers_[HF_DN] = BTagCalibrationReader(&calib_csvv2, BTagEntry::OP_RESHAPING, method, "down_hf");

    readers_[HFSTAT1_UP] = BTagCalibrationReader(&calib_csvv2, BTagEntry::OP_RESHAPING, method, "up_hfstats1");
    readers_[HFSTAT1_DN] = BTagCalibrationReader(&calib_csvv2, BTagEntry::OP_RESHAPING, method, "down_hfstats1");

    readers_[HFSTAT2_UP] = BTagCalibrationReader(&calib_csvv2, BTagEntry::OP_RESHAPING, method, "up_hfstats2");
    readers_[HFSTAT2_DN] = BTagCalibrationReader(&calib_csvv2, BTagEntry::OP_RESHAPING, method, "down_hfstats2");

    readers_[LFSTAT1_UP] = BTagCalibrationReader(&calib_csvv2, BTagEntry::OP_RESHAPING, method, "up_lfstats1");
    readers_[LFSTAT1_DN] = BTagCalibrationReader(&calib_csvv2, BTagEntry::OP_RESHAPING, method, "down_lfstats1");

    readers_[LFSTAT2_UP] = BTagCalibrationReader(&calib_csvv2, BTagEntry::OP_RESHAPING, method, "up_lfstats2");
    readers_[LFSTAT2_DN] = BTagCalibrationReader(&calib_csvv2, BTagEntry::OP_RESHAPING, method, "down_lfstats2");

    readers_[CFERR1_UP] = BTagCalibrationReader(&calib_csvv2, BTagEntry::OP_RESHAPING, method, "up_cferr1");
    readers_[CFERR1_DN] = BTagCalibrationReader(&calib_csvv2, BTagEntry::OP_RESHAPING, method, "down_cferr1");

    readers_[CFERR2_UP] = BTagCalibrationReader(&calib_csvv2, BTagEntry::OP_RESHAPING, method, "up_cferr2");
    readers_[CFERR2_DN] = BTagCalibrationReader(&calib_csvv2, BTagEntry::OP_RESHAPING, method, "down_cferr2");
  }
  else {
    const auto inputFileHF = edm::FileInPath("CATTools/CatAnalyzer/data/scaleFactors/"+rootFileNameHF).fullPath();
    const auto inputFileLF = edm::FileInPath("CATTools/CatAnalyzer/data/scaleFactors/"+rootFileNameLF).fullPath();

    TFile f_CSVwgt_HF(inputFileHF.c_str());
    TFile f_CSVwgt_LF(inputFileLF.c_str());

    fillCSVHistos(&f_CSVwgt_HF, &f_CSVwgt_LF, 6);
  }
}

double CSVWeightEvaluator::operator()(const cat::Jet& jet, const int unc) const
{
  if ( readers_.empty() ) return computeWeightFromROOT(jet, unc);
  return computeWeightFromCSV(jet, unc);
}

double CSVWeightEvaluator::computeWeightFromROOT(const cat::Jet& jet, const int unc) const
{
  const double pt = std::min(jet.pt(), 999.);
  const double aeta = std::abs(jet.eta());
  if ( pt <= 20 or aeta >= 2.4 ) return 1;

  double csv = jet.bDiscriminator(BTAG_CSVv2);
  if ( csv < 0.0 ) csv = -0.05;
  else if ( csv > 1.0 ) csv = 1.0;

  const int flav = std::abs(jet.hadronFlavour());

  int iSysHF = 0;
  switch (unc) {
    case CSVWeightEvaluator::JES_UP: iSysHF = 1; break; // JESUp
    case CSVWeightEvaluator::JES_DN: iSysHF = 2; break; // JESDown
    case CSVWeightEvaluator::LF_UP: iSysHF = 3; break; // LFUp
    case CSVWeightEvaluator::LF_DN: iSysHF = 4; break; // LFDown
    case CSVWeightEvaluator::HFSTAT1_UP: iSysHF = 5; break; // Stats1Up
    case CSVWeightEvaluator::HFSTAT1_DN: iSysHF = 6; break; // Stats1Down
    case CSVWeightEvaluator::HFSTAT2_UP: iSysHF = 7; break; // Stats2Up
    case CSVWeightEvaluator::HFSTAT2_DN: iSysHF = 8; break; // Stats2Down
    default: iSysHF = 0;  break; // NoSys
  }

  int iSysC = 0;
  switch (unc) {
    case CSVWeightEvaluator::CFERR1_UP: iSysC = 1; break;
    case CSVWeightEvaluator::CFERR1_DN: iSysC = 2; break;
    case CSVWeightEvaluator::CFERR2_UP: iSysC = 3; break;
    case CSVWeightEvaluator::CFERR2_DN: iSysC = 4; break;
    default:  iSysC = 0; break;
  }

  int iSysLF = 0;
  switch (unc) {
    case CSVWeightEvaluator::JES_UP: iSysLF = 1; break; // JESUp
    case CSVWeightEvaluator::JES_DN: iSysLF = 2; break; // JESDown
    case CSVWeightEvaluator::HF_UP: iSysLF = 3; break; // HFUp
    case CSVWeightEvaluator::HF_DN: iSysLF = 4; break; // HFDown
    case CSVWeightEvaluator::LFSTAT1_UP: iSysLF = 5; break; // Stats1Up
    case CSVWeightEvaluator::LFSTAT1_DN: iSysLF = 6; break; // Stats1Down
    case CSVWeightEvaluator::LFSTAT2_UP: iSysLF = 7; break; // Stats2Up
    case CSVWeightEvaluator::LFSTAT2_DN: iSysLF = 8; break; // Stats2Down
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

double CSVWeightEvaluator::computeWeightFromCSV(const cat::Jet& jet, const int unc) const
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

// fill the histograms (done once)
void CSVWeightEvaluator::fillCSVHistos(TFile *fileHF, TFile *fileLF, int nHFptBins)
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

