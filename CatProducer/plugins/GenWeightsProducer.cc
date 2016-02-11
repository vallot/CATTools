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
  const edm::InputTag lheLabel_;
  const edm::EDGetTokenT<LHERunInfoProduct> lheRunToken_;
  const edm::EDGetTokenT<LHEEventProduct> lheToken_;
  const edm::EDGetTokenT<GenEventInfoProduct> genInfoToken_;

  std::set<int> scaleWeightIdxs_, pdfWeightIdxs_;
};

GenWeightsProducer::GenWeightsProducer(const edm::ParameterSet& pset):
  enforceUnitGenWeight_(pset.getParameter<bool>("enforceUnitGenWeight")),
  doLOPDFReweight_(pset.getParameter<bool>("doLOPDFReweight")),
  lheLabel_(pset.getParameter<edm::InputTag>("lheEvent")),
  lheRunToken_(consumes<LHERunInfoProduct, edm::InRun>(pset.getParameter<edm::InputTag>("lheEvent"))),
  lheToken_(consumes<LHEEventProduct>(pset.getParameter<edm::InputTag>("lheEvent"))),
  genInfoToken_(consumes<GenEventInfoProduct>(pset.getParameter<edm::InputTag>("genEventInfo")))
{
  std::string pdfName, generatedPdfName;
  if ( doLOPDFReweight_ and pset.existsAs<std::string>("generatedPdfName") )
  {
    generatedPdfName = pset.getParameter<std::string>("generatedPdfName");
    pdfName = pset.getParameter<std::string>("pdfName");
    if ( generatedPdfName != pdfName ) reweightToNewPDF_ = true;
  }

  produces<string, edm::InRun>("combineScaleBy");
  produces<string, edm::InRun>("combinePDFBy");

  produces<float>("genWeight");
  produces<float>("lheWeight");
  produces<vfloat>("scaleWeights");
  produces<vfloat>("pdfWeights");
  produces<vfloat>("otherWeights"); // Should be empty mostly
  produces<int>("id1");
  produces<int>("id2");
  produces<float>("x1");
  produces<float>("x2");
  produces<float>("Q");

  if ( doLOPDFReweight_ )
  {
    LHAPDF::initPDFSet(1, pdfName.c_str());
    if ( reweightToNewPDF_ ) LHAPDF::initPDFSet(2, generatedPdfName.c_str());

    usesResource(); // FIXME What is the resource name of LHAPDF?
  }
}

void GenWeightsProducer::beginRunProduce(edm::Run& run, const edm::EventSetup&)
{
  vstring weightTypes;
  //vvstring weightParams;
  scaleWeightIdxs_.clear();
  pdfWeightIdxs_.clear();

  std::auto_ptr<string> combineScaleBy(new string);
  std::auto_ptr<string> combinePDFBy(new string);

  do {
    // Workaround found in HN, physicstools #3437
    edm::Handle<LHERunInfoProduct> lheHandle;
    //run.getByToken(lheRunToken_, lheHandle);
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
    for ( string line : weightHeader->lines() )
    {
      const auto s0 = line.find_first_of("<");
      const auto s1 = line.find_last_of(">");
      if ( s0 == string::npos or s1 == string::npos ) continue;
      line = line.substr(s0, s1-s0+1);
      contents += line + "\n";
    }
    contents += "\n</lhe>\n"; // Close the root node
    TDOMParser xmlParser; xmlParser.SetValidate(false);
    xmlParser.ParseBuffer(contents.c_str(), contents.size());
    if ( !xmlParser.GetXMLDocument() ) break;

    // XML is ready. Browser the xmldoc and find nodes with weight information
    TXMLNode* topNode = xmlParser.GetXMLDocument()->GetRootNode();
    int weightTotalSize = 0;
    for ( TXMLNode* grpNode = topNode->GetChildren(); grpNode != 0; grpNode = grpNode->GetNextNode() )
    {
      if ( string(grpNode->GetNodeName()) != "weightgroup" ) continue;
      auto weightTypeObj = (TXMLAttr*)grpNode->GetAttributes()->FindObject("name");
      if ( !weightTypeObj ) weightTypeObj = (TXMLAttr*)grpNode->GetAttributes()->FindObject("type"); // FIXME: this may not needed - double check LHE header doc.
      if ( !weightTypeObj ) continue;

      weightTypes.push_back(weightTypeObj->GetValue());
      int weightType = 0;
      if ( weightTypes.back().substr(0, 5) == "scale" ) weightType = 1;
      else if ( weightTypes.back().substr(0, 3) == "PDF" ) weightType = 2;
      int weightSize = 0;
      //weightParams.push_back(vstring());
      for ( TXMLNode* weightNode = grpNode->GetChildren(); weightNode != 0; weightNode = weightNode->GetNextNode() )
      {
        if ( string(weightNode->GetNodeName()) != "weight" ) continue;
        ++weightSize;
        if ( weightSize == 1 ) continue; // Skip the first one of the weight group since it is the nominal value.
        //weightParams.back().push_back(weightNode->GetText());
        if ( weightType == 1 ) scaleWeightIdxs_.insert(weightTotalSize);
        else if ( weightType == 2 ) pdfWeightIdxs_.insert(weightTotalSize);
        ++weightTotalSize;
      }

      auto weightCombineByObj = (TXMLAttr*)grpNode->GetAttributes()->FindObject("combine");
      if ( weightCombineByObj ) {
        if ( weightType == 1 ) *combineScaleBy = weightCombineByObj->GetValue();
        else if ( weightType == 2 ) *combinePDFBy = weightCombineByObj->GetValue();
      }
    }

  } while ( false );

  run.put(combineScaleBy, "combineScaleBy");
  run.put(combinePDFBy, "combinePDFBy");
}

void GenWeightsProducer::produce(edm::Event& event, const edm::EventSetup& eventSetup)
{
  float lheWeight = 1, genWeight = 1;
  std::auto_ptr<vfloat> scaleWeights(new vfloat);
  std::auto_ptr<vfloat> pdfWeights(new vfloat);
  std::auto_ptr<vfloat> otherWeights(new vfloat);
  if ( event.isRealData() )
  {
    pdfWeights->push_back(1); // no reweighting
    event.put(std::auto_ptr<float>(new float(lheWeight)), "lheWeight");
    event.put(std::auto_ptr<float>(new float(genWeight)), "genWeight");
    event.put(scaleWeights, "scaleWeights");
    event.put(pdfWeights, "pdfWeights");
    event.put(otherWeights, "otherWeights");
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
      for ( size_t i=0; i<lheHandle->weights().size(); ++i )
      {
        const double w0 = lheHandle->weights().at(i).wgt;
        const double w = w0/(enforceUnitGenWeight_ ? std::abs(lheWeight) : originalWeight);
        if ( scaleWeightIdxs_.find(i) != scaleWeightIdxs_.end() ) scaleWeights->push_back(w);
        else if ( pdfWeightIdxs_.find(i) != pdfWeightIdxs_.end() ) pdfWeights->push_back(w);
        else otherWeights->push_back(w);
      }
    }
    else
    {
      for ( size_t i=0; i<genInfoHandle->weights().size(); ++i )
      {
        const double w0 = genInfoHandle->weights().at(i);
        const double w = w0/(enforceUnitGenWeight_ ? std::abs(genWeight) : originalWeight);
        if ( scaleWeightIdxs_.find(i) != scaleWeightIdxs_.end() ) scaleWeights->push_back(w);
        else if ( pdfWeightIdxs_.find(i) != pdfWeightIdxs_.end() ) pdfWeights->push_back(w);
        else otherWeights->push_back(w);
      }
    }
  }

  event.put(std::auto_ptr<float>(new float(lheWeight)), "lheWeight");
  event.put(std::auto_ptr<float>(new float(genWeight)), "genWeight");
  event.put(scaleWeights, "scaleWeights");
  event.put(pdfWeights, "pdfWeights");
  event.put(otherWeights, "otherWeights");
  event.put(std::auto_ptr<int>(new int(id1)), "id1");
  event.put(std::auto_ptr<int>(new int(id2)), "id2");
  event.put(std::auto_ptr<float>(new float(x1)), "x1");
  event.put(std::auto_ptr<float>(new float(x2)), "x2");
  event.put(std::auto_ptr<float>(new float(q)), "Q");

}

DEFINE_FWK_MODULE(GenWeightsProducer);

