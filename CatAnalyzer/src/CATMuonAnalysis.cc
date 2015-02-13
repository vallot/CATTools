// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ServiceRegistry/interface/Service.h"

#include "DataFormats/Common/interface/Association.h"
#include "DataFormats/Common/interface/RefToPtr.h"

#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "DataFormats/PatCandidates/interface/Muon.h"
#include "CATTools/DataFormats/interface/Muon.h"
#include "DataFormats/Math/interface/LorentzVector.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "DataFormats/Candidate/interface/Candidate.h"
#include "CommonTools/Utils/interface/StringObjectFunction.h"
#include "DataFormats/PatCandidates/interface/LookupTableRecord.h"

#include "TH1F.h"
#include "TFile.h"
#include "TTree.h"

using namespace std;
//
// class decleration
//
class CATMuonAnalysis : public edm::EDAnalyzer {
  public:
    explicit CATMuonAnalysis(const edm::ParameterSet&);
    ~CATMuonAnalysis();

  private:
    virtual void beginJob() ;
    virtual void analyze(const edm::Event&, const edm::EventSetup&);
    virtual void endJob() ;

    // ----------member data ---------------------------

    edm::EDGetTokenT<edm::View<cat::Muon> > src_;

    TH1F* phi;
    TH1F* eta;
    TH1F* pt;
    TH1F* chIso;
    TH1F* nhIso;
    TH1F* phIso;
    TH1F* puIso;
    
};

CATMuonAnalysis::CATMuonAnalysis(const edm::ParameterSet& iConfig):
  src_(consumes<edm::View<cat::Muon> >(iConfig.getParameter<edm::InputTag>("src")))
{
  edm::Service<TFileService> fs;
  
  phi   = fs->make<TH1F>("phi","phi",400,0,4);
  eta   = fs->make<TH1F>("eta","eta",400,0,4);
  pt    = fs->make<TH1F>("pt","pt",400,0,4);
  chIso = fs->make<TH1F>("chIso","chIso",400,0,4);
  nhIso = fs->make<TH1F>("nhIso","nhIso",400,0,4);
  phIso = fs->make<TH1F>("phIso","phIso",400,0,4);
  puIso = fs->make<TH1F>("puIso","puIso",400,0,4);

}

CATMuonAnalysis::~CATMuonAnalysis()
{
  // do anything here that needs to be done at desctruction time
  // (e.g. close files, deallocate resources etc.)

}

void CATMuonAnalysis::beginJob(){
  //Add event and RUN BRANCHING         
}

void CATMuonAnalysis::endJob(){

}

void CATMuonAnalysis::analyze( const edm::Event& iEvent, const edm::EventSetup& iSetup){

  using namespace edm;
  using namespace std;
  using namespace reco;

  Handle<View<cat::Muon> > src;
  iEvent.getByToken(src_, src);

  for (unsigned int i = 0; i < src->size() ; i++) {
    const cat::Muon & muon = src->at(i);

    phi  ->Fill( muon.phi() );
    eta  ->Fill( muon.eta() );
    pt   ->Fill( muon.pt() );
    chIso->Fill( muon.chargedHadronIso() );
    nhIso->Fill( muon.neutralHadronIso() );
    phIso->Fill( muon.photonIso() );
    puIso->Fill( muon.puChargedHadronIso() );
  }

}

//define this as a plug-in
DEFINE_FWK_MODULE(CATMuonAnalysis);

