#include "FWCore/Framework/interface/one/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ServiceRegistry/interface/Service.h"

#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "CATTools/DataFormats/interface/GenJet.h"
#include "CATTools/DataFormats/interface/GenTop.h"
#include "CATTools/DataFormats/interface/Muon.h"
#include "CATTools/DataFormats/interface/Electron.h"
#include "CATTools/DataFormats/interface/Jet.h"
#include "CATTools/DataFormats/interface/MET.h"
#include "CATTools/DataFormats/interface/GenWeights.h"
#include "CATTools/CommonTools/interface/AnalysisHelper.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "TH1.h"
#include "TTree.h"
#include "TLorentzVector.h"

using namespace cat;

class testAnalyzer : public edm::one::EDAnalyzer<edm::one::SharedResources> {
public:
  explicit testAnalyzer(const edm::ParameterSet&);
  ~testAnalyzer();

private:
  virtual void analyze(const edm::Event&, const edm::EventSetup&);

  //----------------------------------------------------------------
  bool IsSelectMuon    (const cat::Muon     & i_muon_candidate);
  bool IsVetoMuon      (const cat::Muon     & i_muon_candidate);
  bool IsSelectElectron(const cat::Electron & i_electron_candidate);
  bool IsVetoElectron  (const cat::Electron & i_electron_candidate);
  float weight_valid   (float weight);
  float jer_valid      (float weight);

  bool isMC_;
  int TTbarMC_; // 0->No ttbar, 1->ttbar Signal

  edm::EDGetTokenT<float>                   genWeightToken_;
  edm::EDGetTokenT<std::vector<float>>      scaleUpWeightToken_, scaleDownWeightToken_;
  edm::EDGetTokenT<double>                  prefweightToken_, prefweightupToken_, prefweightdownToken_;
  edm::EDGetTokenT<int>                     genttbarHiggsCatToken_;
  edm::EDGetTokenT<cat::GenTopCollection>   genttbarCatToken_;
  edm::EDGetTokenT<cat::MuonCollection>     muonToken_;
  edm::EDGetTokenT<cat::ElectronCollection> electronToken_;
  edm::EDGetTokenT<cat::JetCollection>      jetToken_;
  edm::EDGetTokenT<cat::METCollection>      metToken_;
  edm::EDGetTokenT<int>                     pvToken_, nTrueVertToken_, recoFiltersToken_, recoFiltersMCToken_;
  edm::EDGetTokenT<int>                     trigMuFiltersToken_, trigElFiltersToken_, trigElHTFiltersToken_;

  TTree *tree;

  int b_Event, b_Run, b_Lumi_Number;
  float b_GenWeight;
  std::vector<float> *b_ScaleWeight;

  int b_nGoodPV, b_nTruePV;
  std::vector<double> *b_PrefireWeight;

  int b_Channel, b_GenHiggsCatID;
  std::vector<float> *b_GenCone_gJet_pt;
  float b_GenLepton_pt, b_GenNu_pt, b_GenTop1_pt;

  float b_Lepton_pt, b_MET;
  float b_addbjet1_pt, b_Hbjet1_pt, b_Jet_pt;
  int b_Jet_partonFlavour, b_Jet_hadronFlavour;
  float b_Jet_JES_Up, b_Jet_JES_Down, b_Jet_JER_Up, b_Jet_JER_Nom, b_Jet_JER_Down;
  float b_Jet_deepCSV, b_Jet_deepJet;
  float b_Jet_deepCvsL, b_Jet_deepCvsB, b_Jet_deepJetCvsL, b_Jet_deepJetCvsB;

  // Histograms: Number of Events and Weights
  TH1D *EventInfo;

};

testAnalyzer::testAnalyzer(const edm::ParameterSet& iConfig):
  TTbarMC_    (iConfig.getUntrackedParameter<int>("TTbarSampleLabel", 0))
{
  // Weights
  auto genWeightLabel = iConfig.getParameter<edm::InputTag>("genWeightLabel");
  genWeightToken_       = consumes<float>             (edm::InputTag(genWeightLabel.label()));
  scaleUpWeightToken_   = consumes<std::vector<float>>(edm::InputTag(genWeightLabel.label(), "scaleup"));
  scaleDownWeightToken_ = consumes<std::vector<float>>(edm::InputTag(genWeightLabel.label(), "scaledown"));
  // Prefire (16 and 17 only)
  prefweightToken_      = consumes<double>(edm::InputTag("prefiringweight:nonPrefiringProb"));
  prefweightupToken_    = consumes<double>(edm::InputTag("prefiringweight:nonPrefiringProbUp"));
  prefweightdownToken_  = consumes<double>(edm::InputTag("prefiringweight:nonPrefiringProbDown"));
  // GEN and ttbar Categorization
  genttbarCatToken_      = consumes<cat::GenTopCollection>(iConfig.getParameter<edm::InputTag>("genttbarCatLabel"));
  genttbarHiggsCatToken_ = consumes<int>                  (iConfig.getParameter<edm::InputTag>("genHiggsCatLabel"));
  // Object Collections
  muonToken_        = consumes<cat::MuonCollection>    (iConfig.getParameter<edm::InputTag>("muonLabel"));
  electronToken_    = consumes<cat::ElectronCollection>(iConfig.getParameter<edm::InputTag>("electronLabel"));
  jetToken_         = consumes<cat::JetCollection>     (iConfig.getParameter<edm::InputTag>("jetLabel"));
  metToken_         = consumes<cat::METCollection>     (iConfig.getParameter<edm::InputTag>("metLabel"));
  pvToken_          = consumes<int>                    (iConfig.getParameter<edm::InputTag>("pvLabel"));
  // Trigger
  recoFiltersMCToken_ = consumes<int>(iConfig.getParameter<edm::InputTag>("recoFiltersMC"));
  recoFiltersToken_ = consumes<int>(iConfig.getParameter<edm::InputTag>("recoFilters"));
  trigMuFiltersToken_ = consumes<int>(iConfig.getParameter<edm::InputTag>("trigMuFilters"));
  trigElFiltersToken_ = consumes<int>(iConfig.getParameter<edm::InputTag>("trigElFilters"));
  trigElHTFiltersToken_ = consumes<int>(iConfig.getParameter<edm::InputTag>("trigElHTFilters"));

  nTrueVertToken_    = consumes<int>(iConfig.getParameter<edm::InputTag>("nTrueVertLabel"));
 
  b_PrefireWeight  = new std::vector<double>;
  b_ScaleWeight    = new std::vector<float>;
  b_GenCone_gJet_pt   = new std::vector<float>;

  usesResource("TFileService");
  edm::Service<TFileService> fs;
  tree = fs->make<TTree>("tree", "TopTree");

  tree->Branch("event",      &b_Event,       "event/I");
  tree->Branch("run",        &b_Run,         "run/I");
  tree->Branch("luminumber", &b_Lumi_Number, "luminumber/I");
  tree->Branch("genweight",  &b_GenWeight,   "genweight/F");
  tree->Branch("TruePV",     &b_nTruePV,     "TruePV/I");
  tree->Branch("GoodPV",     &b_nGoodPV,     "GoodPV/I");
  tree->Branch("channel",    &b_Channel,     "channel/I");

  tree->Branch("prefireweight", "std::vector<double>",&b_PrefireWeight);
  tree->Branch("scaleweight",   "std::vector<float>", &b_ScaleWeight );

  tree->Branch("MET",       &b_MET,       "MET/F");
  tree->Branch("lepton_pt", &b_Lepton_pt, "lepton_pt/F");

  tree->Branch("jet_pt",          &b_Jet_pt,          "jet_pt/F");
  tree->Branch("jet_deepCSV",     &b_Jet_deepCSV,     "jet_deepCSV/F" );
  tree->Branch("jet_deepJet",     &b_Jet_deepJet,     "jet_deepJet/F" );
  tree->Branch("jet_deepCvsL",    &b_Jet_deepCvsL,    "jet_deepCvsL/F" );
  tree->Branch("jet_deepCvsB",    &b_Jet_deepCvsB,    "jet_deepCvsB/F" );
  tree->Branch("jet_deepJetCvsL", &b_Jet_deepJetCvsL, "jet_deepJetCvsL/F" );
  tree->Branch("jet_deepJetCvsB", &b_Jet_deepJetCvsB, "jet_deepJetCvsB/F" );

  tree->Branch("jet_partonFlavour", &b_Jet_partonFlavour, "jet_partonFlavour/I" );
  tree->Branch("jet_hadronFlavour", &b_Jet_hadronFlavour, "jet_hadronFlavour/I" );

  tree->Branch("jet_JES_Up",   &b_Jet_JES_Up,   "jet_JES_Up/F" );
  tree->Branch("jet_JES_Down", &b_Jet_JES_Down, "jet_JES_Down/F" );
  tree->Branch("jet_JER_Up",   &b_Jet_JER_Up,   "jet_JER_Up/F" );
  tree->Branch("jet_JER_Nom",  &b_Jet_JER_Nom,  "jet_JER_Nom/F" );
  tree->Branch("jet_JER_Down", &b_Jet_JER_Down, "jet_JER_Down/F" );
  tree->Branch("Hbjet1_pt",    &b_Hbjet1_pt,    "Hbjet1_pt/F");

  // GEN Variables (only ttbarSignal)
  if(TTbarMC_ == 1){
    tree->Branch("gencone_gjet_pt", "std::vector<float>", &b_GenCone_gJet_pt);
    tree->Branch("genhiggscatid",   &b_GenHiggsCatID,     "genhiggscatid/I");
    tree->Branch("genlepton_pt",    &b_GenLepton_pt,      "genlepton_pt/F");
    tree->Branch("gennu_pt",        &b_GenNu_pt,          "gennu_pt/F");
    tree->Branch("gentop1_pt",      &b_GenTop1_pt,        "gentop1_pt/F");
    tree->Branch("addbjet1_pt",     &b_addbjet1_pt,       "addbjet1_pt/F");
  }

  EventInfo = fs->make<TH1D>("EventInfo","Event Information",3,0,3);
  EventInfo->GetXaxis()->SetBinLabel(1,"Number of Events");
  EventInfo->GetXaxis()->SetBinLabel(2,"Sum of Weights");
}


testAnalyzer::~testAnalyzer()
{
  delete b_PrefireWeight;
  delete b_ScaleWeight;
  delete b_GenCone_gJet_pt;
}

void testAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  using namespace edm;

  b_PrefireWeight->clear();
  b_ScaleWeight  ->clear();
  b_GenCone_gJet_pt  ->clear();
  b_addbjet1_pt = b_Hbjet1_pt = -1.0 ;

  // Event Info
  isMC_ = !iEvent.isRealData();

  b_Event        = iEvent.id().event();
  b_Run          = iEvent.id().run();
  b_Lumi_Number  = iEvent.luminosityBlock();
  EventInfo->Fill(0.5, 1.0);         // Number of Events

  bool METfiltered = false;
  if( isMC_ ) {

    edm::Handle<int> nTrueVertHandle;
    iEvent.getByToken( nTrueVertToken_, nTrueVertHandle );
    b_nTruePV = *nTrueVertHandle;

    edm::Handle<float> genWeight;
    iEvent.getByToken(genWeightToken_, genWeight);
    b_GenWeight = *genWeight;
    EventInfo->Fill(1.5, b_GenWeight); // Sum of aMC@NLO Weights

    edm::Handle<double> prefweight;
    iEvent.getByToken(prefweightToken_, prefweight );
    b_PrefireWeight->push_back(*prefweight);

    edm::Handle<double> prefweight_up;
    iEvent.getByToken(prefweightupToken_, prefweight_up );
    b_PrefireWeight->push_back(*prefweight_up);

    edm::Handle<double> prefweight_down;
    iEvent.getByToken(prefweightdownToken_, prefweight_down );
    b_PrefireWeight->push_back(*prefweight_down);

    edm::Handle<int> recoFiltersMCHandle;
    iEvent.getByToken(recoFiltersMCToken_, recoFiltersMCHandle);
    METfiltered = *recoFiltersMCHandle == 0 ? true : false;

  }
  else{
    b_ScaleWeight->push_back(1.0);
    b_GenWeight = 1.0;
    b_nTruePV = 0;

    edm::Handle<int> recoFiltersHandle;
    iEvent.getByToken(recoFiltersToken_, recoFiltersHandle);
    METfiltered = *recoFiltersHandle == 0 ? true : false;
  }

  // Gen ttbar info
  if( TTbarMC_ == 1 ) {
    edm::Handle<std::vector<float>> scaleUpWeightsHandle, scaleDownWeightsHandle;
    iEvent.getByToken(scaleUpWeightToken_,   scaleUpWeightsHandle);
    iEvent.getByToken(scaleDownWeightToken_, scaleDownWeightsHandle);

    b_ScaleWeight->push_back( weight_valid(scaleUpWeightsHandle  ->at(0)) ); // muR=Nom  muF=Up
    b_ScaleWeight->push_back( weight_valid(scaleDownWeightsHandle->at(0)) ); // muR=Nom  muF=Down
  }

  if( TTbarMC_ > 0 ) {
    edm::Handle<int> genttbarHiggsCatHandle;
    iEvent.getByToken( genttbarHiggsCatToken_, genttbarHiggsCatHandle );
    b_GenHiggsCatID = *genttbarHiggsCatHandle;

    edm::Handle<cat::GenTopCollection> genttbarConeCat;
    iEvent.getByToken( genttbarCatToken_, genttbarConeCat );

    b_addbjet1_pt  = genttbarConeCat->begin()->addbJets1().Pt();
    b_Hbjet1_pt  = genttbarConeCat->begin()->HbJets1().Pt();
    b_GenTop1_pt = genttbarConeCat->begin()->topquark1().Pt();

    math::XYZTLorentzVector gJetGenCone[7];
    gJetGenCone[0] = genttbarConeCat->begin()->bJetsFromTop1();
    gJetGenCone[1] = genttbarConeCat->begin()->bJetsFromTop2();
    gJetGenCone[2] = genttbarConeCat->begin()->JetsFromW1();
    gJetGenCone[3] = genttbarConeCat->begin()->JetsFromW2();
    gJetGenCone[4] = genttbarConeCat->begin()->addJets1();
    gJetGenCone[5] = genttbarConeCat->begin()->addJets2();
    gJetGenCone[6] = genttbarConeCat->begin()->FCNCupJet();

    for( int ijGT=0; ijGT<7; ijGT++ )
      b_GenCone_gJet_pt  ->push_back(gJetGenCone[ijGT].Pt());
    
    if( isMC_ && TTbarMC_== 1 ){
      if( genttbarConeCat->begin()->lepton1().pt() != 0. ){
        b_GenLepton_pt  = genttbarConeCat->begin()->lepton1().Pt();
        b_GenNu_pt      = genttbarConeCat->begin()->nu1().Pt();
      }
      else{
        b_GenLepton_pt  = genttbarConeCat->begin()->lepton2().Pt();
        b_GenNu_pt      = genttbarConeCat->begin()->nu2().Pt();
      }

    }// if(GENTTbarMCTree_)
  } // if(TTbarMC==0)

  // Primary Vertex Info
  edm::Handle<int> pvHandle;
  iEvent.getByToken( pvToken_, pvHandle );
  b_nGoodPV = *pvHandle;

  // Missing E_T
  Handle<cat::METCollection> MET;
  iEvent.getByToken(metToken_, MET);
  b_MET = MET->begin()->pt();

  // Electrons
  cat::ElectronCollection selectedElectrons;
  cat::ElectronCollection vetoElectrons;
  Handle<cat::ElectronCollection> electrons;
  iEvent.getByToken(electronToken_, electrons);
  for( unsigned int i = 0; i < electrons->size() ; i++ ){
    const cat::Electron & electron = electrons->at(i);
    if     ( IsSelectElectron( electron ) ) selectedElectrons.push_back( electron );
    else if( IsVetoElectron( electron ) ) vetoElectrons.push_back( electron ); // does not Include selected electrons
  }

  // Muons
  cat::MuonCollection selectedMuons;
  cat::MuonCollection vetoMuons;
  Handle<cat::MuonCollection> muons;
  iEvent.getByToken(muonToken_, muons);
  for( unsigned int i = 0; i < muons->size() ; i++ ){
    const cat::Muon & muon = muons->at(i);
    if     ( IsSelectMuon(muon) ) selectedMuons.push_back(muon);
    else if( IsVetoMuon(muon) )   vetoMuons.push_back(muon); // does not Include selected muons
  }

  // Channel Selection
  TLorentzVector lepton;

  int ch_tag  = 999;
  if( selectedMuons.size()     == 1 &&
      vetoMuons.size()         == 0 &&
      selectedElectrons.size() == 0 &&
      vetoElectrons.size()     == 0){
    lepton.SetPtEtaPhiE(selectedMuons[0].pt(), selectedMuons[0].eta(), selectedMuons[0].phi(), selectedMuons[0].energy());
    ch_tag = 0; //muon + jets
  }

  if( selectedMuons.size()     == 0 &&
      vetoMuons.size()         == 0 &&
      selectedElectrons.size() == 1 &&
      vetoElectrons.size()     == 0){
    ch_tag = 1; //electron + jets
    lepton.SetPtEtaPhiE(selectedElectrons[0].pt(), selectedElectrons[0].eta(), selectedElectrons[0].phi(), selectedElectrons[0].energy());
  }

  // HLTrigger
  bool EvTrigger = false;
  bool IsTriggerMu = false;
  bool IsTriggerEl = false;
  bool IsTriggerElHT = false;

  edm::Handle<int> trigMuFiltersHandle;
  iEvent.getByToken(trigMuFiltersToken_, trigMuFiltersHandle);
  IsTriggerMu = *trigMuFiltersHandle == 0 ? false : true;

  edm::Handle<int> trigElFiltersHandle;
  iEvent.getByToken(trigElFiltersToken_, trigElFiltersHandle);
  IsTriggerEl = *trigElFiltersHandle == 0 ? false : true;

  edm::Handle<int> trigElHTFiltersHandle;
  iEvent.getByToken(trigElHTFiltersToken_, trigElHTFiltersHandle);
  IsTriggerElHT = *trigElHTFiltersHandle == 0 ? false : true;

  if( (ch_tag == 0 && IsTriggerMu) ||
      (ch_tag == 1 && (IsTriggerEl || IsTriggerElHT) ) ){
    EvTrigger = true;
  }

  // Fill Tree with events that have ONLY one lepton
  if( b_nGoodPV < 1 ) ch_tag = 999; 

  if( ch_tag<2 && EvTrigger && !METfiltered ){ // Single lepton event
    b_Channel  = ch_tag;
    b_Lepton_pt  = lepton.Pt();

    // Jets
    Handle<cat::JetCollection> jets;
    iEvent.getByToken(jetToken_, jets);

    for( unsigned int i = 0; i < jets->size() ; i++ ){
      const cat::Jet & jet = jets->at(i);

      bool goodJet  = false;
      bool cleanJet = false;

      // Jet Selection (pT>20GeV to take into account SYST Variations)
      if( std::abs(jet.eta()) < 2.4 && jet.pt() > 20. && jet.tightLepVetoJetID() ) goodJet = true;

      // Jet Cleaning
      TLorentzVector vjet;
      vjet.SetPtEtaPhiE(jet.pt(), jet.eta(), jet.phi(), jet.energy());
      double dr_LepJet = vjet.DeltaR(lepton);
      if( dr_LepJet > 0.4 ) cleanJet = true;

      if( goodJet && cleanJet ){
        b_Jet_pt = static_cast<float>(jet.pt());
        b_Jet_partonFlavour = jet.partonFlavour();
        b_Jet_hadronFlavour = jet.hadronFlavour();

        b_Jet_deepCSV = jet.bDiscriminator(BTAG_DeepCSVb)+jet.bDiscriminator(BTAG_DeepCSVbb);
        b_Jet_deepJet = jet.bDiscriminator(BTAG_DeepJetb)+jet.bDiscriminator(BTAG_DeepJetbb)+jet.bDiscriminator(BTAG_DeepJetlepb);
        b_Jet_deepCvsL = jet.bDiscriminator(BTAG_DeepCSVc)/(jet.bDiscriminator(BTAG_DeepCSVc)+jet.bDiscriminator(BTAG_DeepCSVudsg));
        b_Jet_deepCvsB = jet.bDiscriminator(BTAG_DeepCSVc)/(jet.bDiscriminator(BTAG_DeepCSVc)+jet.bDiscriminator(BTAG_DeepCSVb)+jet.bDiscriminator(BTAG_DeepCSVbb));
        b_Jet_deepJetCvsL = jet.bDiscriminator(BTAG_DeepJetc)/(jet.bDiscriminator(BTAG_DeepJetc)+jet.bDiscriminator(BTAG_DeepJetuds)+jet.bDiscriminator(BTAG_DeepJetg));
        b_Jet_deepJetCvsB = jet.bDiscriminator(BTAG_DeepJetc)/(jet.bDiscriminator(BTAG_DeepJetc)+jet.bDiscriminator(BTAG_DeepJetb)+jet.bDiscriminator(BTAG_DeepJetbb)+jet.bDiscriminator(BTAG_DeepJetlepb));

        if( isMC_ ){
          b_Jet_JES_Up   = jet.shiftedEnUp();
          b_Jet_JES_Down = jet.shiftedEnDown();

          // Ref: https://twiki.cern.ch/twiki/bin/view/CMS/JetResolution
          b_Jet_JER_Up   = jer_valid(jet.smearedResUp());
          b_Jet_JER_Nom  = jer_valid(jet.smearedRes());
          b_Jet_JER_Down = jer_valid(jet.smearedResDown());

        } // if(isMC_)
        break; //Check only first jet

      }// if(GoodJets)
    }// for(AllJets)

    tree->Fill();
        
  }// if(ch_tag)
}

//------------- Good Muon Selection -----------------------
bool testAnalyzer::IsSelectMuon(const cat::Muon & i_muon_candidate){

  bool GoodMuon=true;
  GoodMuon &= (i_muon_candidate.passed(reco::Muon::CutBasedIdTight|reco::Muon::PFIsoTight));
  GoodMuon &= (i_muon_candidate.isPFMuon());           // PF
  GoodMuon &= (i_muon_candidate.pt()> 20);             // pT
  GoodMuon &= (std::abs(i_muon_candidate.eta())< 2.4); // eta

  return GoodMuon;
}

//------------- Loose Muon Selection -----------------------
bool testAnalyzer::IsVetoMuon(const cat::Muon & i_muon_candidate){

  bool GoodMuon=true;
  GoodMuon &= (i_muon_candidate.passed(reco::Muon::CutBasedIdLoose|reco::Muon::PFIsoLoose));
  GoodMuon &= (i_muon_candidate.isPFMuon());           // PF
  GoodMuon &= (i_muon_candidate.pt()> 15);             // pT
  GoodMuon &= (std::abs(i_muon_candidate.eta())< 2.4); // eta

  return GoodMuon;
}

//------------- Good Electron Selection -----------------------
bool testAnalyzer::IsSelectElectron(const cat::Electron & i_electron_candidate){

  bool GoodElectron=true;
  GoodElectron &= (i_electron_candidate.isPF() );                // PF
  GoodElectron &= (i_electron_candidate.pt() > 20);              // pT
  GoodElectron &= (std::abs(i_electron_candidate.eta()) < 2.4);  // eta
  GoodElectron &= (std::abs(i_electron_candidate.scEta()) < 1.4442 || // eta Super-Cluster
                   std::abs(i_electron_candidate.scEta()) > 1.566);
  GoodElectron &= i_electron_candidate.electronID("cutBasedElectronID-Fall17-94X-V2-tight") > 0.0;

  return GoodElectron;

}

//------------- Loose Electron Selection -----------------------
bool testAnalyzer::IsVetoElectron(const cat::Electron & i_electron_candidate){

  bool GoodElectron=true;
  GoodElectron &= (i_electron_candidate.isPF() );                // PF
  GoodElectron &= (i_electron_candidate.pt() > 15);              // pT
  GoodElectron &= (std::abs(i_electron_candidate.eta()) < 2.4);  // eta
  GoodElectron &= (std::abs(i_electron_candidate.scEta()) < 1.4442 || // eta Super-Cluster
                   std::abs(i_electron_candidate.scEta()) > 1.566);

  GoodElectron &= i_electron_candidate.electronID("cutBasedElectronID-Fall17-94X-V2-veto") > 0.0;

  return GoodElectron;
}

float testAnalyzer::weight_valid(float weight){
  if( std::isnan(weight) or std::isinf(weight) ) return 1.0;
  if( (weight < 0.01) or (weight > 100) )        return 1.0;
  return weight;
};

float testAnalyzer::jer_valid(float weight){
  if( std::isnan(weight) or std::isinf(weight) ) return 0.0;
  if( weight < 0.0 )                             return 0.0;
  return weight;
};

//define this as a plug-in
DEFINE_FWK_MODULE(testAnalyzer);
