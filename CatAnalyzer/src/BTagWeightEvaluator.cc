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
    if ( type_ == ITERATIVEFIT ) {
      // Method 1d) using discriminator-dependent scale factors;
      double weight = 1.0;
      for ( auto& jet : jets ) weight *= getSF(jet, unc);
      return weight;
    }
    else if ( type_ == CSVWEIGHT ) {
      csvHelper_->getCSVWeight(jets, unc == 0 ? 0 : unc+7);
    }
  }
  return 1.0;
}

void BTagWeightEvaluator::initCSVWeight(const bool useCSVHelper, const string btagName)
{
  method_ = 4;

  if ( useCSVHelper ) {
    type_ = CSVWEIGHT;
    csvHelper_.reset(new CSVHelper());
  }
  else {
    type_ = ITERATIVEFIT;

    string csvFileName;
    if      ( btagName == "csvv2" ) {
      btagAlgo_ = BTAG_CSVv2 ;
      csvFileName = "CSVv2_ichep.csv";
      //csvFileName = "ttH_BTV_CSVv2_13TeV_2015D_20151120.csv";
    }
    else if ( btagName == "mva" ) {
      btagAlgo_ = BTAG_cMVAv2;
      csvFileName = "cMVAv2_ichep.csv";
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
    csvFileName = "CSVv2_ichep.csv";
  }
  else if ( btagName == "mva" ) {
    btagAlgo_ = BTAG_cMVAv2;
    csvFileName = "cMVAv2_ichep.csv";
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

  if ( type_ == ITERATIVEFIT ) {
    if      ( discr < -1.0 ) discr = -0.05;
    else if ( discr >  1.0 ) discr = 1.0;

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

