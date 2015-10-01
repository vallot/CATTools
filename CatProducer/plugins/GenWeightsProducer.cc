#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/stream/EDProducer.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"

#include "SimDataFormats/GeneratorProducts/interface/LHEEventProduct.h"
#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"
#include "DataFormats/Common/interface/View.h"

#include <LHAPDF/LHAPDF.h>

#include <memory>
#include <vector>
#include <string>

using namespace std;

class GenWeightsProducer : public edm::stream::EDProducer<>
{
public:
  GenWeightsProducer(const edm::ParameterSet& pset);
  void produce(edm::Event& event, const edm::EventSetup& eventSetup) override;

private:
  const bool enforceUnitGenWeight_;
  bool doReweightPdf_;
  const std::string pdfName_;
  std::string generatedPdfName_;
  const edm::EDGetTokenT<LHEEventProduct> lheToken_;
  const edm::EDGetTokenT<GenEventInfoProduct> genInfoToken_;
  const unsigned int lheWeightIndex_, genWeightIndex_;
};

GenWeightsProducer::GenWeightsProducer(const edm::ParameterSet& pset):
  enforceUnitGenWeight_(pset.getParameter<bool>("enforceUnitGenWeight")),
  pdfName_(pset.getParameter<std::string>("pdfName")),
  lheToken_(consumes<LHEEventProduct>(pset.getParameter<edm::InputTag>("lheEvent"))),
  genInfoToken_(consumes<GenEventInfoProduct>(pset.getParameter<edm::InputTag>("genEventInfo"))),
  lheWeightIndex_(pset.getParameter<int>("lheWeightIndex")),
  genWeightIndex_(pset.getParameter<int>("genWeightIndex"))
{
  doReweightPdf_ = false;
  if ( pset.existsAs<std::string>("generatedPdfName") )
  {
    generatedPdfName_ = pset.getParameter<std::string>("generatedPdfName");
    if ( generatedPdfName_ != pdfName_ ) doReweightPdf_ = true;
  }

  produces<float>("genWeight");
  produces<float>("lheWeight");
  produces<std::vector<float> >("pdfWeights");
  produces<int>("id1");
  produces<int>("id2");
  produces<float>("x1");
  produces<float>("x2");
  produces<float>("Q");

  LHAPDF::initPDFSet(1, pdfName_.c_str());
  if ( doReweightPdf_ ) LHAPDF::initPDFSet(2, generatedPdfName_.c_str());
}

void GenWeightsProducer::produce(edm::Event& event, const edm::EventSetup& eventSetup)
{
  float lheWeight = 1, genWeight = 1;
  std::auto_ptr<std::vector<float> > pdfWeights(new std::vector<float>);
  if ( event.isRealData() )
  {
    pdfWeights->push_back(1); // no reweighting
    event.put(std::auto_ptr<float>(new float(lheWeight)), "lheWeight");
    event.put(std::auto_ptr<float>(new float(genWeight)), "genWeight");
    event.put(pdfWeights, "pdfWeights");
  }

  edm::Handle<LHEEventProduct> lheHandle;
  event.getByToken(lheToken_, lheHandle);
  edm::Handle<GenEventInfoProduct> genInfoHandle;
  event.getByToken(genInfoToken_, genInfoHandle);

  // Generator weights
  //   enforceUnitGenWeight == true  : genWeights are scaled to +1 or -1 by themselves
  //   enforceUnitGenWeight == false : genWeights are scaled by 1/originalWeight.
  //                                   do not scale weight if LHE is not available.
  double originalWeight = 1;
  if ( lheHandle.isValid() and lheHandle->weights().size() > lheWeightIndex_ ) {
    lheWeight = lheHandle->weights().at(lheWeightIndex_).wgt;
    originalWeight = std::abs(lheHandle->originalXWGTUP());
  }
  if ( genInfoHandle->weights().size() > genWeightIndex_ ) genWeight = genInfoHandle->weights().at(genWeightIndex_);
  if ( enforceUnitGenWeight_ ) {
    lheWeight = lheWeight == 0 ? 0 : lheWeight/std::abs(lheWeight);
    genWeight = genWeight == 0 ? 0 : genWeight/std::abs(genWeight);
  }
  else {
    lheWeight = originalWeight == 0 ? 0 : lheWeight/originalWeight;
    genWeight = originalWeight == 0 ? 0 : genWeight/originalWeight;
  }

  const float q = genInfoHandle->pdf()->scalePDF;
  const int id1 = genInfoHandle->pdf()->id.first;
  const int id2 = genInfoHandle->pdf()->id.second;
  const float x1 = genInfoHandle->pdf()->x.first;
  const float x2 = genInfoHandle->pdf()->x.second;

  const int generatedPdfIdx = doReweightPdf_ ? 2 : 1;
  const float xpdf1 = LHAPDF::xfx(generatedPdfIdx, x1, q, id1);
  const float xpdf2 = LHAPDF::xfx(generatedPdfIdx, x2, q, id2);
  const float w0 = xpdf1*xpdf2;

  if ( !doReweightPdf_ ) pdfWeights->push_back(1);
  else
  {
    const float xpdf1_new = LHAPDF::xfx(1, x1, q, id1);
    const float xpdf2_new = LHAPDF::xfx(1, x2, q, id2);
    const float w_new = xpdf1_new*xpdf2_new;
    pdfWeights->push_back(w_new/w0);
  }

  for ( unsigned int i=1, n=LHAPDF::numberPDF(1); i<=n; ++i )
  {
    LHAPDF::usePDFMember(1, i);
    const float xpdf1_syst = LHAPDF::xfx(1, x1, q, id1);
    const float xpdf2_syst = LHAPDF::xfx(1, x2, q, id2);
    pdfWeights->push_back(xpdf1_syst*xpdf2_syst/w0);
  }

  event.put(std::auto_ptr<float>(new float(lheWeight)), "lheWeight");
  event.put(std::auto_ptr<float>(new float(genWeight)), "genWeight");
  event.put(pdfWeights, "pdfWeights");
  event.put(std::auto_ptr<int>(new int(id1)), "id1");
  event.put(std::auto_ptr<int>(new int(id2)), "id2");
  event.put(std::auto_ptr<float>(new float(x1)), "x1");
  event.put(std::auto_ptr<float>(new float(x2)), "x2");
  event.put(std::auto_ptr<float>(new float(q)), "Q");

}

DEFINE_FWK_MODULE(GenWeightsProducer);

