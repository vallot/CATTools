#include "CATTools/CatAnalyzer/interface/BTagWeightEvaluator.h"
#include "FWCore/ParameterSet/interface/FileInPath.h"

using namespace cat;
using namespace std;

double BTagWeightEvaluator::eventWeight(const cat::JetCollection& jets, const int unc) const
{
  if ( method_ == 3 ) {
    // Method 1c) using scale factors only
    const int njet = jets.size();
    if ( njet == 0 ) return 1;

    double w0n = 1; // w(0|n) : weight of 0-tag among n-jets
    double w1n = 0; // w(1|n) : weight of 1-tag among n-jets
    for ( int i=0; i<njet; ++i ) {
      w0n *= 1-getSF(jets[i], unc);
      double prodw1n = 1;
      for ( int j=0; j<njet; ++j ) {
        if ( i == j ) prodw1n *= getSF(jets[j], unc);
        else prodw1n *= 1-getSF(jets[j], unc);
      }
      w1n += prodw1n;
    }

    if      ( minNbjet_ == 0 ) return w0n;
    else if ( minNbjet_ == 1 ) return 1-w0n;
    else if ( minNbjet_ == 2 ) return 1-w0n-w1n;
  }
  else if ( method_ == 4 ) {
    // Method 1d) using discriminator-dependent scale factors;
    double weight = 1.0;
    for ( auto& jet : jets ) weight *= getSF(jet, unc);
    return weight;
  }
  return 1.0;
}

void BTagWeightEvaluator::initCSVWeight(const bool isFromROOT, const string btagName)
{
  method_ = 4;
  if ( isFromROOT ) {
    type_ = CSVWEIGHT;
//    const string rootFileNameHF = "csv_rwt_fit_hf_2016_01_28.root";
//    const string rootFileNameLF = "csv_rwt_fit_lf_2016_01_28.root";
    const string rootFileNameHF = "csv_rwt_fit_hf_76x_2016_02_08.root";
    const string rootFileNameLF = "csv_rwt_fit_lf_76x_2016_02_08.root";
    const auto inputFileHF = edm::FileInPath("CATTools/CatAnalyzer/data/scaleFactors/"+rootFileNameHF).fullPath();
    const auto inputFileLF = edm::FileInPath("CATTools/CatAnalyzer/data/scaleFactors/"+rootFileNameLF).fullPath();

    TFile f_CSVwgt_HF(inputFileHF.c_str());
    TFile f_CSVwgt_LF(inputFileLF.c_str());

    fillCSVHistos(&f_CSVwgt_HF, &f_CSVwgt_LF, 6);
  }
  else {
    type_ = ITERATIVEFIT;

    string csvFileName;
    if      ( btagName == "csvv2" ) {
      btagAlgo_ = BTAG_CSVv2 ;
      csvFileName = "CSVv2_prelim.csv";
      //csvFileName = "ttH_BTV_CSVv2_13TeV_2015D_20151120.csv";
    }
    else if ( btagName == "mva" ) {
      btagAlgo_ = BTAG_cMVAv2;
      csvFileName = "cMVAv2_prelim.csv";
    }
    //else if ( btagName == "jp"    ) btagAlgo_ = BTAG_JP    ;
    else btagAlgo_ = "undefined"; // FIXME: Eventually raise error somewhere?

    const auto csvFile = edm::FileInPath("CATTools/CatAnalyzer/data/scaleFactors/"+csvFileName).fullPath();
    BTagCalibration calib(btagName, csvFile);
    uncNames_ = {
      "central", "up_jes", "down_jes",
      "up_lf", "down_lf", "up_hf", "down_hf", 
      "up_hfstats1", "down_hfstats1", "up_hfstats2", "down_hfstats2",
      "up_lfstats1", "down_lfstats1", "up_lfstats2", "down_lfstats2",
      "up_cferr1", "down_cferr1", "up_cferr2", "down_cferr2"
    };
    for ( unsigned int i=0; i<uncNames_.size(); ++i ) {
      readers_[i] = BTagCalibrationReader(&calib, BTagEntry::OP_RESHAPING, "iterativefit", uncNames_[i]); 
    }
  }
}

void BTagWeightEvaluator::init(const int method,
                               const string btagName, BTagEntry::OperatingPoint operationPoint, int minNbjet)
{
  type_ = STANDARD;

  minNbjet_ = minNbjet;
  method_ = method;

  string csvFileName;
  if      ( btagName == "csvv2" ) {
    btagAlgo_ = BTAG_CSVv2 ;
    csvFileName = "CSVv2_prelim.csv";
  }
  else if ( btagName == "mva" ) {
    btagAlgo_ = BTAG_cMVAv2;
    csvFileName = "cMVAv2_prelim.csv";
  }
  //else if ( btagName == "jp"    ) btagAlgo_ = BTAG_JP    ;
  else btagAlgo_ = "undefined"; // FIXME: Eventually raise error somewhere?

  const auto csvFile = edm::FileInPath("CATTools/CatAnalyzer/data/scaleFactors/"+csvFileName).fullPath();
  BTagCalibration calib(btagName, csvFile);
  uncNames_ = {"central", "up", "down"};

  readers_[0] = BTagCalibrationReader(&calib, operationPoint, "incl", "central");
  readers_[1] = BTagCalibrationReader(&calib, operationPoint, "incl", "up"     );
  readers_[2] = BTagCalibrationReader(&calib, operationPoint, "incl", "down"   );

  readers_[3] = BTagCalibrationReader(&calib, operationPoint, "mujets", "central");
  readers_[4] = BTagCalibrationReader(&calib, operationPoint, "mujets", "up"     );
  readers_[5] = BTagCalibrationReader(&calib, operationPoint, "mujets", "down"   );
}

double BTagWeightEvaluator::getSF(const cat::Jet& jet, const int unc) const
{
  const double pt = std::min(jet.pt(), 999.);
  const double eta = jet.eta();
  const double aeta = std::abs(eta);
  if ( pt <= 20 or aeta >= 2.4 ) return 1.0;

  double discr = jet.bDiscriminator(btagAlgo_);
  const int flav = std::abs(jet.hadronFlavour());

  if ( type_ == CSVWEIGHT or type_ == ITERATIVEFIT ) {
    if      ( discr < -1.0 ) discr = -0.05;
    else if ( discr >  1.0 ) discr = 1.0;

    if ( type_ == CSVWEIGHT ) return getCSVWeightSFFromROOT(pt, aeta, flav, discr, unc);

    // Special care for the flavour dependent SFs
    int uncKey = unc;
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

    auto readerItr = readers_.find(uncKey);
    if ( readerItr == readers_.end() ) return 1.0;

    BTagEntry::JetFlavor jf = BTagEntry::FLAV_UDSG;
    if      ( flav == 5 ) jf = BTagEntry::FLAV_B;
    else if ( flav == 4 ) jf = BTagEntry::FLAV_C;

    return readerItr->second.eval(jf, aeta, pt, discr);
  }
  else {
    int uncKey = unc;
    if ( flav == 5 or flav == 4 ) uncKey += 3; // Use mujets reader for b flavour and c-flavour

    auto readerItr = readers_.find(uncKey);
    if ( readerItr == readers_.end() ) return 1.0;

    if      ( flav == 5 ) return readerItr->second.eval(BTagEntry::FLAV_B, eta, pt, discr);
    else if ( flav == 4 ) return readerItr->second.eval(BTagEntry::FLAV_C, eta, pt, discr);

    return readerItr->second.eval(BTagEntry::FLAV_UDSG, aeta, pt, discr);
  }

  return 1.;
};

double BTagWeightEvaluator::getCSVWeightSFFromROOT(const double pt, const double aeta,
                                                   const int flav, const double discr,
                                                   const int unc) const
{
  // original code : https://github.com/cms-ttH/MiniAOD/blob/master/MiniAODHelper/interface/CSVHelper.h
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
    int useCSVBin = (discr >= 0.) ? h_csv_wgt_hf[iSysHF][iPt]->FindBin(discr) : 1;
    return (float) h_csv_wgt_hf[iSysHF][iPt]->GetBinContent(useCSVBin);
  } else if (flav == 4) {
    if(iPt>=nHFptBins_) iPt=nHFptBins_-1;
    int useCSVBin = (discr >= 0.) ? c_csv_wgt_hf[iSysC][iPt]->FindBin(discr) : 1;
    return (float) c_csv_wgt_hf[iSysC][iPt]->GetBinContent(useCSVBin);
  } else {
    if (iPt >= 3)  iPt = 3; /// [30-40], [40-60] and [60-10000] only 3 Pt bins for lf
    int useCSVBin = (discr >= 0.) ? h_csv_wgt_lf[iSysLF][iPt][iEta]->FindBin(discr) : 1;
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

