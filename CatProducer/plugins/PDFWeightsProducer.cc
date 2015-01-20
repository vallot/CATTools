#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"

#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"
#include "DataFormats/Common/interface/View.h"

#include <LHAPDF/LHAPDF.h>

#include <memory>
#include <vector>
#include <string>

using namespace std;

class PDFWeightsProducer : public edm::EDProducer
{
public:
  PDFWeightsProducer(const edm::ParameterSet& pset);
  void beginJob() override;
  void produce(edm::Event& event, const edm::EventSetup& eventSetup) override;

private:
  bool doReweightPdf_;
  std::string pdfName_;
  std::string generatedPdfName_;
  edm::InputTag genInfoToken_;
};

PDFWeightsProducer::PDFWeightsProducer(const edm::ParameterSet& pset)
{
  genInfoToken_ = consumes<GenEventInfoProduct>(pset.getParameter<edm::InputTag>("genEventInfo"));
  pdfName_ = pset.getParameter<std::string>("pdfName");

  doReweightPdf_ = false;
  if ( pset.existsAs<std::string>("generatedPdfName") )
  {
    generatedPdfName_ = pset.getParameter<std::string>("generatedPdfName");
    if ( generatedPdfName_ != pdfName_ ) doReweightPdf_ = true;
  }

  produces<std::vector<double> >();
  produces<int>("id1");
  produces<int>("id2");
  produces<double>("x1");
  produces<double>("x2");
  produces<double>("Q");
}

void PDFWeightsProducer::beginJob()
{
  LHAPDF::initPDFSet(1, pdfName_.c_str());
  if ( doReweightPdf_ ) LHAPDF::initPDFSet(2, generatedPdfName_.c_str());
}

void PDFWeightsProducer::produce(edm::Event& event, const edm::EventSetup& eventSetup)
{
  std::auto_ptr<std::vector<double> > weights(new std::vector<double>);
  if ( event.isRealData() )
  {
    weights->push_back(1); // no reweighting
    event.put(weights);
    return;
  }

  edm::Handle<GenEventInfoProduct> genInfoHandle;
  event.getByLabel(genInfoToken_, genInfoHandle);

  const float q = genInfoHandle->pdf()->scalePDF;
  const int id1 = genInfoHandle->pdf()->id.first;
  const int id2 = genInfoHandle->pdf()->id.second;
  const double x1 = genInfoHandle->pdf()->x.first;
  const double x2 = genInfoHandle->pdf()->x.second;

  const int generatedPdfIdx = doReweightPdf_ ? 2 : 1;
  const double xpdf1 = LHAPDF::xfx(generatedPdfIdx, x1, q, id1);
  const double xpdf2 = LHAPDF::xfx(generatedPdfIdx, x2, q, id2);
  const double w0 = xpdf1*xpdf2;

  if ( !doReweightPdf_ ) weights->push_back(1);
  else
  {
    const double xpdf1_new = LHAPDF::xfx(1, x1, q, id1);
    const double xpdf2_new = LHAPDF::xfx(1, x2, q, id2);
    const double w_new = xpdf1_new*xpdf2_new;
    weights->push_back(w_new/w0);
  }

  for ( unsigned int i=1, n=LHAPDF::numberPDF(1); i<=n; ++i )
  {
    LHAPDF::usePDFMember(1, i);
    const double xpdf1_syst = LHAPDF::xfx(1, x1, q, id1);
    const double xpdf2_syst = LHAPDF::xfx(1, x2, q, id2);
    weights->push_back(xpdf1_syst*xpdf2_syst/w0);
  }

  event.put(weights);
  event.put(std::auto_ptr<int>(new int(id1)), "id1");
  event.put(std::auto_ptr<int>(new int(id2)), "id2");
  event.put(std::auto_ptr<double>(new double(x1)), "x1");
  event.put(std::auto_ptr<double>(new double(x2)), "x2");
  event.put(std::auto_ptr<double>(new double(q)), "Q");

}

DEFINE_FWK_MODULE(PDFWeightsProducer);

