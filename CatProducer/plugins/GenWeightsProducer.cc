#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/one/EDProducer.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"

#include "SimDataFormats/GeneratorProducts/interface/LHERunInfoProduct.h"
#include "SimDataFormats/GeneratorProducts/interface/LHEEventProduct.h"
#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"
#include "DataFormats/Common/interface/View.h"

#include <LHAPDF/LHAPDF.h>
#include <TDOMParser.h>
#include <TXMLNode.h>
#include <TXMLAttr.h>

#include <memory>
#include <vector>
#include <string>

using namespace std;

class GenWeightsProducer : public edm::one::EDProducer<edm::one::SharedResources, edm::BeginRunProducer>
{
public:
  GenWeightsProducer(const edm::ParameterSet& pset);
  void beginRunProduce(edm::Run& run, const edm::EventSetup&) override;
  void produce(edm::Event& event, const edm::EventSetup& eventSetup) override;

  typedef std::vector<int> vint;
  typedef std::vector<float> vfloat;
  typedef std::vector<std::string> vstring;
  typedef std::vector<vstring> vvstring;

private:
  const bool enforceUnitGenWeight_;
  const bool doLOPDFReweight_;
  bool reweightToNewPDF_;
  const std::string pdfName_;
  std::string generatedPdfName_;
  const edm::InputTag lheLabel_;
  const edm::EDGetTokenT<LHEEventProduct> lheToken_;
  const edm::EDGetTokenT<GenEventInfoProduct> genInfoToken_;

  vint weightIdxUB_;
};

GenWeightsProducer::GenWeightsProducer(const edm::ParameterSet& pset):
  enforceUnitGenWeight_(pset.getParameter<bool>("enforceUnitGenWeight")),
  doLOPDFReweight_(pset.getParameter<bool>("doLOPDFReweight")),
  pdfName_(pset.getParameter<std::string>("pdfName")),
  lheLabel_(pset.getParameter<edm::InputTag>("lheEvent")),
  lheToken_(consumes<LHEEventProduct>(pset.getParameter<edm::InputTag>("lheEvent"))),
  genInfoToken_(consumes<GenEventInfoProduct>(pset.getParameter<edm::InputTag>("genEventInfo")))
{
  if ( doLOPDFReweight_ and pset.existsAs<std::string>("generatedPdfName") )
  {
    generatedPdfName_ = pset.getParameter<std::string>("generatedPdfName");
    if ( generatedPdfName_ != pdfName_ ) reweightToNewPDF_ = true;
  }

  produces<vstring, edm::InRun>("weightTypes");
  produces<vint, edm::InRun>("weightIdxUB");
  produces<vvstring, edm::InRun>("weightParams");

  produces<float>("genWeight");
  produces<float>("lheWeight");
  produces<vfloat>("pdfWeights");
  produces<int>("id1");
  produces<int>("id2");
  produces<float>("x1");
  produces<float>("x2");
  produces<float>("Q");

  if ( doLOPDFReweight_ )
  {
    LHAPDF::initPDFSet(1, pdfName_.c_str());
    if ( reweightToNewPDF_ ) LHAPDF::initPDFSet(2, generatedPdfName_.c_str());

    usesResource(); // FIXME What is the resource name of LHAPDF?
  }
}

void GenWeightsProducer::beginRunProduce(edm::Run& run, const edm::EventSetup&)
{
  std::auto_ptr<vstring> weightTypes(new vstring);
  std::auto_ptr<vint> weightIdxUB(new vint);
  std::auto_ptr<vvstring> weightParams(new vvstring);

  do {
    edm::Handle<LHERunInfoProduct> lheHandle;
    run.getByLabel(lheLabel_, lheHandle);
    if ( !lheHandle.isValid() ) break;

    // Find interested LHE header
    auto weightHeader = lheHandle->headers_end();
    for ( auto itr = lheHandle->headers_begin(); itr != lheHandle->headers_end(); ++ itr )
    {
      // "initrwgt" is the header for the weights
      if ( string(itr->tag()) == "initrwgt" )
      {
        weightHeader = itr;
        break;
      }
    }
    if ( weightHeader == lheHandle->headers_end() ) break;

    // Read LHE using the XML parser
    string contents = "<lhe>\n"; // Need root node
    contents.reserve(10000); // ~50 char per line, >100 weights
    for ( auto line : weightHeader->lines() )
    {
      line = line.substr(line.find_first_not_of(" \t\r\t"), line.find_last_not_of(" \n\r\t"));
      if ( line.empty() or line[0] != '<' or line[line.size()-1] != '>' ) continue;
      contents += line + "\n";
    }
    contents += "\n</lhe>\n"; // Close the root node
    TDOMParser xmlParser; xmlParser.SetValidate(false);
    xmlParser.ParseBuffer(contents.c_str(), contents.size());
    if ( !xmlParser.GetXMLDocument() ) break;

    // XML is ready. Browser the xmldoc and find nodes with weight information
    TXMLNode* topNode = xmlParser.GetXMLDocument()->GetRootNode();
    for ( TXMLNode* grpNode = topNode->GetChildren(); grpNode != 0; grpNode = grpNode->GetNextNode() )
    {
      if ( string(grpNode->GetNodeName()) != "weightgroup" ) continue;
      auto weightTypeObj = (TXMLAttr*)grpNode->GetAttributes()->FindObject("type");
      if ( !weightTypeObj ) continue;

      weightTypes->push_back(weightTypeObj->GetValue());
      int weightSize = 0;
      weightParams->push_back(vstring());
      for ( TXMLNode* weightNode = grpNode->GetChildren(); weightNode != 0; weightNode = weightNode->GetNextNode() )
      {
        if ( string(weightNode->GetNodeName()) != "weight" ) continue;
        weightParams->back().push_back(weightNode->GetText());
        ++weightSize;
      }
      if ( weightIdxUB->empty() ) weightIdxUB->push_back(weightSize);
      else weightIdxUB->push_back(weightIdxUB->back()+weightSize);
    }

  } while ( false );

  weightIdxUB_.clear();
  std::copy(weightIdxUB->begin(), weightIdxUB->end(), std::back_inserter(weightIdxUB_));

  run.put(weightTypes, "weightTypes");
  run.put(weightIdxUB, "weightIdxUB");
  run.put(weightParams, "weightParams");
}

void GenWeightsProducer::produce(edm::Event& event, const edm::EventSetup& eventSetup)
{
  float lheWeight = 1, genWeight = 1;
  std::auto_ptr<vfloat> pdfWeights(new vfloat);
  if ( event.isRealData() )
  {
    pdfWeights->push_back(1); // no reweighting
    event.put(std::auto_ptr<float>(new float(lheWeight)), "lheWeight");
    event.put(std::auto_ptr<float>(new float(genWeight)), "genWeight");
    event.put(pdfWeights, "pdfWeights");
    event.put(std::auto_ptr<int>(new int(0)), "id1");
    event.put(std::auto_ptr<int>(new int(0)), "id2");
    event.put(std::auto_ptr<float>(new float(0)), "x1");
    event.put(std::auto_ptr<float>(new float(0)), "x2");
    event.put(std::auto_ptr<float>(new float(0)), "Q");
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
  if ( lheHandle.isValid() and !lheHandle->weights().empty() ) {
    lheWeight = lheHandle->weights().at(0).wgt;
    originalWeight = std::abs(lheHandle->originalXWGTUP());
  }
  genWeight = genInfoHandle->weight();
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

  if ( doLOPDFReweight_ )
  {
    const int generatedPdfIdx = reweightToNewPDF_ ? 2 : 1;
    const float xpdf1 = LHAPDF::xfx(generatedPdfIdx, x1, q, id1);
    const float xpdf2 = LHAPDF::xfx(generatedPdfIdx, x2, q, id2);
    const float w0 = xpdf1*xpdf2;

    const float xpdf1_new = LHAPDF::xfx(1, x1, q, id1);
    const float xpdf2_new = LHAPDF::xfx(1, x2, q, id2);
    const float w_new = xpdf1_new*xpdf2_new;
    pdfWeights->push_back(w_new/w0);

    for ( unsigned int i=1, n=LHAPDF::numberPDF(1); i<=n; ++i )
    {
      LHAPDF::usePDFMember(1, i);
      const float xpdf1_syst = LHAPDF::xfx(1, x1, q, id1);
      const float xpdf2_syst = LHAPDF::xfx(1, x2, q, id2);
      pdfWeights->push_back(xpdf1_syst*xpdf2_syst/w0);
    }
  }
  else
  {
    if ( lheHandle.isValid() )
    {
      for ( size_t i=1; i<lheHandle->weights().size(); ++i )
      {
        const double w = lheHandle->weights().at(i).wgt;
        if ( enforceUnitGenWeight_ ) pdfWeights->push_back(w/std::abs(lheWeight));
        else pdfWeights->push_back(w/originalWeight);
      }
    }
    else
    {
      for ( size_t i=1; i<genInfoHandle->weights().size(); ++i )
      {
        const double w = genInfoHandle->weights().at(i);
        if ( enforceUnitGenWeight_ ) pdfWeights->push_back(w/std::abs(genWeight));
        else pdfWeights->push_back(w/originalWeight);
      }
    }
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

