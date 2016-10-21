#include "FWCore/Framework/interface/one/EDAnalyzer.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ServiceRegistry/interface/Service.h"

#include "DataFormats/Common/interface/Association.h"
#include "DataFormats/Common/interface/RefToPtr.h"

#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "DataFormats/PatCandidates/interface/Muon.h"
#include "CATTools/DataFormats/interface/Electron.h"
#include "CATTools/DataFormats/interface/Muon.h"
#include "CATTools/DataFormats/interface/Jet.h"
#include "CATTools/DataFormats/interface/FatJet.h"
#include "CATTools/DataFormats/interface/MET.h"
#include "DataFormats/Math/interface/LorentzVector.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "DataFormats/Candidate/interface/Candidate.h"
#include "CommonTools/Utils/interface/StringObjectFunction.h"
#include "DataFormats/PatCandidates/interface/LookupTableRecord.h"

#include "TH1F.h"
#include "TH2F.h"
#include "TTree.h"

using namespace std;

typedef struct 
{
        float pt; 
        float eta;
        float tau1;
        float tau2;
        float tau3;
        float sdmass;
        float pmass;
        float mass;
} Fatjetstruct;

//
// class decleration
//
class CATFatjetAnalysis : public edm::one::EDAnalyzer<edm::one::SharedResources> {
  public:
    explicit CATFatjetAnalysis(const edm::ParameterSet&);
    ~CATFatjetAnalysis();

  private:
    virtual void analyze(const edm::Event&, const edm::EventSetup&) override;

    // ----------member data ---------------------------

    edm::EDGetTokenT<cat::MuonCollection>          muonToken_;
    edm::EDGetTokenT<cat::ElectronCollection>      electronToken_;
    edm::EDGetTokenT<cat::JetCollection>           jetToken_;
    edm::EDGetTokenT<cat::FatJetCollection>        fatjetToken_;
    edm::EDGetTokenT<cat::METCollection>           metToken_;

    TH1F* phi;
    TH1F* eta;
    TH1F* pt;
    TH2F* htau123;
    TH1F* hpmass;


    Fatjetstruct afatjet;

    TTree *tree_;

};

CATFatjetAnalysis::CATFatjetAnalysis(const edm::ParameterSet& iConfig)
{
    muonToken_         = consumes<cat::MuonCollection>          (iConfig.getParameter<edm::InputTag>("muonLabel"));
    electronToken_     = consumes<cat::ElectronCollection>      (iConfig.getParameter<edm::InputTag>("electronLabel"));
    jetToken_          = consumes<cat::JetCollection>           (iConfig.getParameter<edm::InputTag>("jetLabel"));
    fatjetToken_          = consumes<cat::FatJetCollection>           (iConfig.getParameter<edm::InputTag>("fatjetLabel"));
    metToken_          = consumes<cat::METCollection>           (iConfig.getParameter<edm::InputTag>("metLabel"));
  usesResource("TFileService");
  edm::Service<TFileService> fs;

  phi   = fs->make<TH1F>("phi","phi",100,0,4);
  eta   = fs->make<TH1F>("eta","eta",100,0,4);
  pt    = fs->make<TH1F>("pt","pt",100,0,1000);
  htau123 = fs->make<TH2F>("htau123","#tau_{23} vs #tau_{1}",30,0,1, 30, 0, 1);
  hpmass = fs->make<TH1F>("pmass", "Pruned mass", 50, 100, 300);
  tree_ = fs->make<TTree>("event", "event");
  tree_->Branch("fatjet", &afatjet, "pt/F:eta:tau1:tau2:tau3:sdmass:pmass:mass");
}
  


CATFatjetAnalysis::~CATFatjetAnalysis()
{
  // do anything here that needs to be done at desctruction time
  // (e.g. close files, deallocate resources etc.)

}

void CATFatjetAnalysis::analyze( const edm::Event& iEvent, const edm::EventSetup& iSetup){

  using namespace edm;
  using namespace std;
  using namespace reco;

    Handle<cat::MuonCollection> muons;
    iEvent.getByToken(muonToken_, muons);
    Handle<cat::ElectronCollection> electrons;
    iEvent.getByToken(electronToken_, electrons);
    Handle<cat::METCollection> MET;
    iEvent.getByToken(metToken_, MET);
    Handle<cat::JetCollection> jets;
    iEvent.getByToken(jetToken_, jets);
    Handle<cat::FatJetCollection> fatjets;
    iEvent.getByToken(fatjetToken_, fatjets);

    //cout << fatjets->size() << endl;
    for (unsigned int i = 0; i < fatjets->size() ; i++) {
        const cat::FatJet & j = fatjets->at(i);

        phi  ->Fill( j.phi() );
        eta  ->Fill( j.eta() );
        pt   ->Fill( j.pt() );
        afatjet.pt = j.pt();
        afatjet.eta = j.eta();
        afatjet.tau1 = j.tau1();
        afatjet.tau2 = j.tau2();
        afatjet.tau3 = j.tau3();
        afatjet.sdmass = j.softdropmass();
        afatjet.pmass = j.prunedmass();
        afatjet.mass = j.mass();
        tree_->Fill();

        if (j.pt()>500.0 && j.mass()>150.0 && j.mass()<200.0)
        {
            htau123->Fill( j.tau1(), j.tau3()/j.tau2());
            hpmass->Fill(j.prunedmass());
        }
    }
}

//define this as a plug-in
DEFINE_FWK_MODULE(CATFatjetAnalysis);

