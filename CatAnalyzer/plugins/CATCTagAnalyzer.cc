#include "FWCore/ServiceRegistry/interface/Service.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "CATTools/DataFormats/interface/Jet.h"

#include "TTree.h"
#include <memory>
using namespace std;
using namespace cat;

class CATCTagAnalyzer : public edm::one::EDAnalyzer<edm::one::SharedResources,edm::one::WatchLuminosityBlocks> {
public:
  virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
  explicit CATCTagAnalyzer(const edm::ParameterSet&);
  void beginLuminosityBlock(const edm::LuminosityBlock& lumi, const edm::EventSetup&) override {};
  void endLuminosityBlock(const edm::LuminosityBlock&, const edm::EventSetup&) override {};
  ~CATCTagAnalyzer() {};


private:

  edm::EDGetTokenT<cat::JetCollection>      jetToken_;
  

  void resetBr();

  std::vector<float> b_ctag_CVL, b_ctag_CVB;

  std::vector<std::vector<int> > cutflow_;
  TTree* ttree_;

  

};
//
// constructors and destructor
//
CATCTagAnalyzer::CATCTagAnalyzer(const edm::ParameterSet& iConfig)
{
  jetToken_  = consumes<cat::JetCollection>(iConfig.getParameter<edm::InputTag>("jets"));

  usesResource("TFileService");
  edm::Service<TFileService> fs;
  ttree_ = fs->make<TTree>("ctag","ctag");
  TTree* tr = ttree_;
  tr->Branch("CvsL","std::vector<float>",&b_ctag_CVL);
  tr->Branch("CvsB","std::vector<float>",&b_ctag_CVB);


}


void CATCTagAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  //const bool runOnMC = !iEvent.isRealData();

  resetBr();

  edm::Handle<cat::JetCollection> jets;    iEvent.getByToken(jetToken_,jets);
  for( const auto& jet : *jets ) {
    float CVL = jet.bDiscriminator("pfCombinedCvsLJetTags");
    float CVB = jet.bDiscriminator("pfCombinedCvsBJetTags");
    //std::cout<<CVL<<"   "<<CVB<<std::endl;
    b_ctag_CVL.push_back(CVL);
    b_ctag_CVB.push_back(CVB);
  }
  ttree_->Fill();

}


void CATCTagAnalyzer::resetBr()
{
  b_ctag_CVL.clear(); b_ctag_CVB.clear();
}

//define this as a plug-in
DEFINE_FWK_MODULE(CATCTagAnalyzer);
