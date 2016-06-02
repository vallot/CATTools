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

#include "CATTools/DataFormats/interface/GenWeights.h"

#include <LHAPDF/LHAPDF.h>
#include <TDOMParser.h>
#include <TXMLNode.h>
#include <TXMLAttr.h>
#include <regex>
#include <boost/lexical_cast.hpp>

#include <memory>
#include <vector>
#include <string>

using namespace std;

class GenWeightsProducer;
class GenWeightsToFlatWeights;

class GenWeightsToFlatWeights : public edm::stream::EDProducer<>
{
public:
  GenWeightsToFlatWeights(const edm::ParameterSet& pset);
  void beginRun(const edm::Run& run, const edm::EventSetup&) override;
  void produce(edm::Event& event, const edm::EventSetup&) override;

  typedef std::vector<float> vfloat;

private:
  const edm::InputTag srcLabel_;
  edm::EDGetTokenT<cat::GenWeights> srcToken_;
  const bool doSaveOthers_, doKeepFirstOnly_;

  std::set<size_t> key_sup_, key_sdn_, key_pdf_, key_oth_;
};

GenWeightsToFlatWeights::GenWeightsToFlatWeights(const edm::ParameterSet& pset):
  srcLabel_(pset.getParameter<edm::InputTag>("src")),
  doSaveOthers_(pset.getParameter<bool>("saveOthers")),
  doKeepFirstOnly_(pset.getParameter<bool>("keepFirstOnly"))
{
  srcToken_ = consumes<cat::GenWeights>(srcLabel_);
  consumes<cat::GenWeightInfo, edm::InRun>(srcLabel_);

  produces<float>();
  produces<vfloat>("scaleup");
  produces<vfloat>("scaledown");
  produces<vfloat>("pdf");
  if ( doSaveOthers_ ) produces<vfloat>("others");
}

void GenWeightsToFlatWeights::beginRun(const edm::Run& run, const edm::EventSetup&)
{
  edm::Handle<cat::GenWeightInfo> srcHandle;
  if ( !run.getByLabel(srcLabel_, srcHandle) ) return;

  for ( int i=0, n=srcHandle->nGroups(); i<n; ++i ) {
    const auto& keys = srcHandle->keys(i);

    string name = srcHandle->name(i);
    std::transform(name.begin(), name.end(), name.begin(), ::toupper);
    if ( name.find("SCALE") != string::npos ) {
      if ( !key_sup_.empty() and !key_sdn_.empty() ) {
        cout << "!!! Duplicated scale variation from " << name << "!!!\n";
        continue;
      }
      const auto& params = srcHandle->params(i);
      for ( int j=1, m=keys.size(); j<m; ++j ) {
        string par = params[j];
        // Skip unphysical combinations
        // By fixed index: up=(1002, 1004, 1005), down=(1003, 1007, 1009), unphysical=(1006, 1008)
        // Or by matching to the muR and muF parameters
        // First, cleanup the string
        std::transform(par.begin(), par.end(), par.begin(), ::toupper);
        par.erase(std::remove(par.begin(), par.end(), '_'), par.end());
        std::replace(par.begin(), par.end(), '=', ' ');
        // Then tokenize parameter string
        std::vector<std::string> tokens;
        std::regex ws_re("\\s+"); // whitespace
        std::copy(std::sregex_token_iterator(par.begin(), par.end(), ws_re, -1),
                  std::sregex_token_iterator(), std::back_inserter(tokens));
        double muR = 0, muF = 0;
        for ( int iToken=0, nToken=tokens.size(); iToken<nToken; ++iToken) {
          if ( tokens[iToken] == "MUR" and iToken+1 < nToken ) {
            const double orig = muR;
            try { muR = boost::lexical_cast<double>(tokens[++iToken]); }
            catch ( boost::bad_lexical_cast ) { muR = orig; }
          }
          if ( tokens[iToken] == "MUF" and iToken+1 < nToken ) {
            const double orig = muF;
            try { muF = boost::lexical_cast<double>(tokens[++iToken]); }
            catch ( boost::bad_lexical_cast ) { muF = orig; }
          }
        }
        if ( (muR == 0 or muF == 0) or (muR > 1.5 and muF < 0.7) or (muR < 0.7 and muF > 1.5 ) ) {
          cout << "@@@ Skipping unphysical parameter" << params[j] << endl;
          continue;
        }

        const size_t key = keys[j];
        if      ( muR > 1.5 or muF > 1.5 ) {
          cout << "@@@ Inserting into scaleup weight[" << (key_sup_.size()) << "]" << params[j] << endl;
          key_sup_.insert(key);
        }
        else if ( muR < 0.7 or muF < 0.7 ) {
          cout << "@@@ Inserting into scaledown weight[" << (key_sdn_.size()) << "]" << params[j] << endl;
          key_sdn_.insert(key);
        }
      }
    }
    else if ( name.find("PDF") != string::npos ) {
      if ( doKeepFirstOnly_ and !key_pdf_.empty() ) {
        cout << "@@@ Skipping PDF weight " << name << " since the first weight group is already set " << endl;
        continue;
      }
      key_pdf_.insert(keys.begin(), keys.end());
    }
    else if ( doSaveOthers_ ) {
      key_oth_.insert(keys.begin(), keys.end());
    }
  }
}

void GenWeightsToFlatWeights::produce(edm::Event& event, const edm::EventSetup&)
{
  std::auto_ptr<float> out_weight(new float(1));
  std::auto_ptr<vfloat> out_sup(new vfloat);
  std::auto_ptr<vfloat> out_sdn(new vfloat);
  std::auto_ptr<vfloat> out_pdf(new vfloat);
  std::auto_ptr<vfloat> out_oth(new vfloat);

  edm::Handle<cat::GenWeights> srcHandle;
  event.getByToken(srcToken_, srcHandle);

  if ( srcHandle.isValid() ) {
    *out_weight = srcHandle->genWeight();
    const auto weights = srcHandle->weights();
    for ( int i=0, n=weights.size(); i<n; ++i ) {
      const auto& w = weights[i];
      if      ( key_sup_.find(i) != key_sup_.end() ) out_sup->push_back(w);
      else if ( key_sdn_.find(i) != key_sdn_.end() ) out_sdn->push_back(w);
      else if ( key_pdf_.find(i) != key_pdf_.end() ) out_pdf->push_back(w);
      else if ( key_oth_.find(i) != key_oth_.end() ) out_oth->push_back(w);
    }
  }

  event.put(out_weight);
  event.put(out_sup, "scaleup");
  event.put(out_sdn, "scaledown");
  event.put(out_pdf, "pdf");
  if ( doSaveOthers_ ) event.put(out_oth, "others");
}

class GenWeightsProducer : public edm::one::EDProducer<edm::one::SharedResources, edm::BeginRunProducer>
{
public:
  GenWeightsProducer(const edm::ParameterSet& pset);
  void beginRunProduce(edm::Run& run, const edm::EventSetup&) override;
  void produce(edm::Event& event, const edm::EventSetup& eventSetup) override;

  typedef std::vector<float> vfloat;
  typedef std::vector<std::string> vstring;
  typedef std::vector<unsigned short> vushort;

private:
  const edm::InputTag lheLabel_;

  const bool enforceUnitGenWeight_;
  const bool doLOPDFReweight_;
  bool reweightToNewPDF_;
  const edm::EDGetTokenT<LHEEventProduct> lheToken_;
  const edm::EDGetTokenT<GenEventInfoProduct> genInfoToken_;

};

GenWeightsProducer::GenWeightsProducer(const edm::ParameterSet& pset):
  lheLabel_(pset.getParameter<edm::InputTag>("lheEvent")),
  enforceUnitGenWeight_(pset.getParameter<bool>("enforceUnitGenWeight")),
  doLOPDFReweight_(pset.getParameter<bool>("doLOPDFReweight")),
  lheToken_(consumes<LHEEventProduct>(pset.getParameter<edm::InputTag>("lheEvent"))),
  genInfoToken_(consumes<GenEventInfoProduct>(pset.getParameter<edm::InputTag>("genEventInfo")))
{
  consumes<LHERunInfoProduct, edm::InRun>(lheLabel_);
  produces<cat::GenWeightInfo, edm::InRun>();
  produces<cat::GenWeights>();

  std::string pdfName, generatedPdfName;
  if ( doLOPDFReweight_ and pset.existsAs<std::string>("generatedPdfName") )
  {
    generatedPdfName = pset.getParameter<std::string>("generatedPdfName");
    pdfName = pset.getParameter<std::string>("pdfName");
    if ( generatedPdfName != pdfName ) reweightToNewPDF_ = true;
  }

  if ( doLOPDFReweight_ )
  {
    LHAPDF::initPDFSet(1, pdfName.c_str());
    if ( reweightToNewPDF_ ) LHAPDF::initPDFSet(2, generatedPdfName.c_str());

    usesResource(); // FIXME What is the resource name of LHAPDF?
  }
}

void GenWeightsProducer::beginRunProduce(edm::Run& run, const edm::EventSetup&)
{
  std::auto_ptr<cat::GenWeightInfo> out_genWeightInfo(new cat::GenWeightInfo);

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

      const string weightTypeStr = weightTypeObj->GetValue();
      vushort keys;
      vstring params;
      for ( TXMLNode* weightNode = grpNode->GetChildren(); weightNode != 0; weightNode = weightNode->GetNextNode() )
      {
        if ( string(weightNode->GetNodeName()) != "weight" ) continue;
        ++weightTotalSize;
        params.push_back(weightNode->GetText());
        keys.push_back(weightTotalSize-1);
      }

      auto weightCombineByObj = (TXMLAttr*)grpNode->GetAttributes()->FindObject("combine");
      string combineBy = weightCombineByObj ? weightCombineByObj->GetValue() : "";
      out_genWeightInfo->addWeightGroup(weightTypeStr, combineBy, params, keys);
    }

  } while ( false );

  run.put(out_genWeightInfo);
}

void GenWeightsProducer::produce(edm::Event& event, const edm::EventSetup& eventSetup)
{
  float lheWeight = 1, genWeight = 1;
  std::auto_ptr<cat::GenWeights> out_genWeights(new cat::GenWeights);

  if ( event.isRealData() ) {
    out_genWeights->setLHEWeight(lheWeight);
    out_genWeights->setGenWeight(genWeight);
    event.put(out_genWeights);
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
  out_genWeights->setInfo(id1, id2, x1, x2, q);
  out_genWeights->setLHEWeight(lheWeight);
  out_genWeights->setGenWeight(genWeight);

  if ( doLOPDFReweight_ )
  {
    const int generatedPdfIdx = reweightToNewPDF_ ? 2 : 1;
    const float xpdf1 = LHAPDF::xfx(generatedPdfIdx, x1, q, id1);
    const float xpdf2 = LHAPDF::xfx(generatedPdfIdx, x2, q, id2);
    const float w0 = xpdf1*xpdf2;

    const float xpdf1_new = LHAPDF::xfx(1, x1, q, id1);
    const float xpdf2_new = LHAPDF::xfx(1, x2, q, id2);
    const float w_new = xpdf1_new*xpdf2_new;
    out_genWeights->addWeight(w_new/w0);

    for ( unsigned int i=1, n=LHAPDF::numberPDF(1); i<=n; ++i )
    {
      LHAPDF::usePDFMember(1, i);
      const float xpdf1_syst = LHAPDF::xfx(1, x1, q, id1);
      const float xpdf2_syst = LHAPDF::xfx(1, x2, q, id2);
      out_genWeights->addWeight(xpdf1_syst*xpdf2_syst/w0);
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
        out_genWeights->addWeight(w);
      }
    }
    else
    {
      for ( size_t i=0; i<genInfoHandle->weights().size(); ++i )
      {
        const double w0 = genInfoHandle->weights().at(i);
        const double w = w0/(enforceUnitGenWeight_ ? std::abs(genWeight) : originalWeight);
        out_genWeights->addWeight(w);
      }
    }
  }

  event.put(out_genWeights);

}

DEFINE_FWK_MODULE(GenWeightsProducer);
DEFINE_FWK_MODULE(GenWeightsToFlatWeights);

