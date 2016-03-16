#include "FWCore/Framework/interface/one/EDAnalyzer.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ServiceRegistry/interface/Service.h"

#include "DataFormats/Common/interface/Association.h"
#include "DataFormats/Common/interface/RefToPtr.h"

#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "DataFormats/PatCandidates/interface/Muon.h"
#include "CATTools/DataFormats/interface/SecVertex.h"
#include "DataFormats/Math/interface/LorentzVector.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "DataFormats/Candidate/interface/Candidate.h"
#include "CommonTools/Utils/interface/StringObjectFunction.h"
#include "DataFormats/PatCandidates/interface/LookupTableRecord.h"

#include "TH1F.h"
#include "TTree.h"
#include "TNtuple.h"
using namespace std;
//
// class decleration
//
class CATDStarAnalysis : public edm::one::EDAnalyzer<edm::one::SharedResources> {
  public:
    explicit CATDStarAnalysis(const edm::ParameterSet&);
    ~CATDStarAnalysis();

  private:
    virtual void analyze(const edm::Event&, const edm::EventSetup&) override;

    // ----------member data ---------------------------

    edm::EDGetTokenT< vector<cat::SecVertex> > D0Src_, DstarSrc_;

    TNtuple *nt, *nt2 ;

};

CATDStarAnalysis::CATDStarAnalysis(const edm::ParameterSet& iConfig):
  D0Src_(consumes<cat::SecVertexCollection>(iConfig.getParameter<edm::InputTag>("D0Src"))),
  DstarSrc_(consumes<cat::SecVertexCollection>(iConfig.getParameter<edm::InputTag>("DstarSrc")))
{
  usesResource("TFileService");
  edm::Service<TFileService> fs;

  nt = fs->make<TNtuple>("D0","D0","D0_Mass");
  nt2 = fs->make<TNtuple>("DStar","DStar","DStar_Mass");
}

CATDStarAnalysis::~CATDStarAnalysis()
{
  // do anything here that needs to be done at desctruction time
  // (e.g. close files, deallocate resources etc.)

}

void CATDStarAnalysis::analyze( const edm::Event& iEvent, const edm::EventSetup& iSetup){

  using namespace edm;
  using namespace std;
  using namespace reco;

  Handle<cat::SecVertexCollection> D0Src;
  iEvent.getByToken(D0Src_, D0Src);
  Handle<cat::SecVertexCollection> DstarSrc;
  iEvent.getByToken(DstarSrc_, DstarSrc);

  for (unsigned int i = 0; i < D0Src->size() ; i++) {
    const cat::SecVertex & D0Cand = D0Src->at(i);
    nt->Fill( D0Cand.mass());
  }
  
  for(unsigned int i=0 ; i< DstarSrc->size() ; i++) {
    const cat::SecVertex & DstarCand = DstarSrc->at(i);
    nt2->Fill(DstarCand.mass());

  }

}

//define this as a plug-in
DEFINE_FWK_MODULE(CATDStarAnalysis);

