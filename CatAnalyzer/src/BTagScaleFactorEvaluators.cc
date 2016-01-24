#include "CATTools/CatAnalyzer/interface/BTagScaleFactorEvaluators.h"
#include "FWCore/ParameterSet/interface/FileInPath.h"

using namespace cat;

using namespace std;
CSVWeightEvaluator::CSVWeightEvaluator()
{
  const char* method = "iterativefit";
  // setup calibration readers (once)
  const auto csvFile = edm::FileInPath("CATTools/CatAnalyzer/data/scaleFactors/ttH_BTV_CSVv2_13TeV_2015D_20151120.csv").fullPath();
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

double CSVWeightEvaluator::operator()(const cat::Jet& jet, const int unc) const
{
  const double pt = std::min(jet.pt(), 999.);
  const double aeta = std::abs(jet.eta());
  if ( pt <= 20 or aeta >= 2.4 ) return 1;

  double csv = jet.bDiscriminator("pfCombinedInclusiveSecondaryVertexV2BJetTags");
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

