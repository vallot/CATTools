/**
  \class    PFMuonProducer PFMuonProducer.h "CATTools/CatProducer/interface/PFMuonProducer.h"
  \brief    PF Muon 
*/


#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "DataFormats/Common/interface/Association.h"
#include "DataFormats/Common/interface/RefToPtr.h"
#include "DataFormats/MuonReco/interface/Muon.h"
#include "DataFormats/MuonReco/interface/MuonFwd.h"

#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidate.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidateFwd.h"
#include "CommonTools/UtilAlgos/interface/StringCutObjectSelector.h"
#include "RecoEcal/EgammaCoreTools/interface/EcalClusterLazyTools.h"

#include "RecoParticleFlow/PFProducer/interface/PFMuonAlgo.h"  //PFMuons

using namespace edm;
using namespace std;

namespace reco {

  class PFMuonProducer : public edm::EDProducer {
    public:
      explicit PFMuonProducer(const edm::ParameterSet & iConfig);
      virtual ~PFMuonProducer() { }

      virtual void produce(edm::Event & iEvent, const edm::EventSetup & iSetup);

    private:
    edm::InputTag src_;
    
      PFMuonAlgo *pfmu_;
  };

} // namespace

reco::PFMuonProducer::PFMuonProducer(const edm::ParameterSet & iConfig) :
  src_(iConfig.getParameter<edm::InputTag>( "src" ))
{
    pfmu_ = new PFMuonAlgo();
    pfmu_->setParameters(iConfig);
    produces<std::vector<reco::Muon> >();
}

void 
reco::PFMuonProducer::produce(edm::Event & iEvent, const edm::EventSetup & iSetup) {


    edm::Handle<reco::MuonCollection> src;
    iEvent.getByLabel(src_, src);

    auto_ptr<vector<reco::Muon> >  out(new vector<reco::Muon>());

    std::cout << "Total Number of Muons = " << src->size() << std::endl; 

    for (reco::MuonCollection::size_type i = 0; i < src->size(); ++i) {
     
      const reco::Muon aRecoMuon = src->at(i);

      reco::MuonRef aRecoMuonRef(src,i);

      //bool pfmuon = pfmu_->isMuon(aRecoMuonRef) || pfmu_->isLooseMuon(aRecoMuonRef);
      bool pfmuon = pfmu_->isMuon(aRecoMuonRef);
      bool defaultPFMuon = aRecoMuon.isPFMuon(); 

      if(defaultPFMuon) std::cout << "Default PF muon at " << i ;

      if(pfmuon){
        std::cout << " ---> new pf muon at " << i ;
        out->push_back(aRecoMuon);
      }

      std::cout << endl;

    }

    iEvent.put(out);
}

#include "FWCore/Framework/interface/MakerMacros.h"
using namespace reco;
DEFINE_FWK_MODULE(PFMuonProducer);
