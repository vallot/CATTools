// -*- C++ -*-
//
// Package:    TtbarSingleLeptonAnalyzer
// Class:      TtbarSingleLeptonAnalyzer
//
/**\class TtbarSingleLeptonAnalyzer TtbarSingleLeptonAnalyzer.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Javier Brochero Cifuentes,512 1-001,+41789428316
//         Created:  Tue Feb  3 09:52:55 CET 2015
// $Id$
//
//

#include "FWCore/Framework/interface/one/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ServiceRegistry/interface/Service.h"

#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "CATTools/DataFormats/interface/GenTop.h"
#include "CATTools/DataFormats/interface/Muon.h"
#include "CATTools/DataFormats/interface/Electron.h"
#include "CATTools/DataFormats/interface/Jet.h"
#include "CATTools/DataFormats/interface/MET.h"

#include "CATTools/CommonTools/interface/AnalysisHelper.h"

#include "CommonTools/UtilAlgos/interface/TFileService.h"

// b-tagging Eff
//#include "CATTools/CatAnalyzer/interface/BTagSFUtil.h"

#include "TH1.h"
#include "TTree.h"
#include "TLorentzVector.h"
//
// class declaration
//

using namespace cat;
using namespace std;

template<class T>
struct bigger_second
: std::binary_function<T,T,bool>
{
   inline bool operator()(const T& lhs, const T& rhs)
   {
      return lhs.second > rhs.second;
   }
};
typedef std::pair<int,double> data_t;

class TtbarSingleLeptonAnalyzer : public edm::one::EDAnalyzer<edm::one::SharedResources> {
public:
  explicit TtbarSingleLeptonAnalyzer(const edm::ParameterSet&);
  ~TtbarSingleLeptonAnalyzer();

private:
  virtual void analyze(const edm::Event&, const edm::EventSetup&);

  //----------------------------------------------------------------
  bool IsTightElectron(const cat::Electron & i_electron_candidate);
  bool IsLooseElectron(const cat::Electron & i_electron_candidate);

  double transverseMass( const reco::Candidate::LorentzVector& lepton, const reco::Candidate::LorentzVector& met); 

  edm::EDGetTokenT<cat::GenTopCollection>          genTopToken_;
  edm::EDGetTokenT<reco::GenParticleCollection>  genToken_;
  edm::EDGetTokenT<cat::MuonCollection>          muonToken_;
  edm::EDGetTokenT<cat::ElectronCollection>      electronToken_;
  edm::EDGetTokenT<cat::JetCollection>           jetToken_;
  edm::EDGetTokenT<cat::METCollection>           metToken_;
  edm::EDGetTokenT<int>                            pvToken_;
  edm::EDGetTokenT<float>                          puWeight_;

  // ----------member data ---------------------------

  TTree *tree;

  //---------------------------------------------------------------------------
  //---------------------------------------------------------------------------
  // Tree Branches
  //---------------------------------------------------------------------------
  //---------------------------------------------------------------------------

  // Event info
  int EVENT; 
  int RUN; 
  int LUMI;

  // PU/Vertices
  float puweight;
  int nvertex;

  // MET
  float MET, MET_phi;

  // Leptons
  float lep_E;
  float lep_pt;
  float lep_relIso;

  int nLep;
  int nLooseLep;
  int nMuon;
  int nElec;

  //WMT

  float mt; 

  // Jets
  std::vector<float> *jet_E;
  std::vector<float> *jet_pt;

  // Flavour
  std::vector<int> *b_Jet_partonFlavour;
  std::vector<int> *b_Jet_hadronFlavour;
  // Smearing and Shifts
  std::vector<float> *b_Jet_smearedRes;
  std::vector<float> *b_Jet_smearedResDown;
  std::vector<float> *b_Jet_smearedResUp;
  std::vector<float> *b_Jet_shiftedEnUp;
  std::vector<float> *b_Jet_shiftedEnDown;
  // b-Jet discriminant
  std::vector<double> *b_Jet_CSV;

  std::vector<int> * csvid;

  //multiplicity
  int nJets30;
  int nbJets30;

  // event categorization
  int dileptonic;
  int semileptonic;
  int visible_ttjj;
  int visible_ttbb;
  int visible_ttbj;
  int visible_ttcc;
  int visible_ttLF;
  int ttjj;
  int ttbb;

};

//
// constants, enums and typedefs
//

//
// static data member definitions
//

//
// constructors and destructor
//
TtbarSingleLeptonAnalyzer::TtbarSingleLeptonAnalyzer(const edm::ParameterSet& iConfig){
  //now do what ever initialization is needed
  genTopToken_      = consumes<cat::GenTopCollection>  (iConfig.getParameter<edm::InputTag>("genTopLabel"));
  genToken_      = consumes<reco::GenParticleCollection>  (iConfig.getParameter<edm::InputTag>("genLabel"));
  muonToken_     = consumes<cat::MuonCollection>          (iConfig.getParameter<edm::InputTag>("muonLabel"));
  electronToken_ = consumes<cat::ElectronCollection>      (iConfig.getParameter<edm::InputTag>("electronLabel"));
  jetToken_      = consumes<cat::JetCollection>           (iConfig.getParameter<edm::InputTag>("jetLabel"));
  metToken_      = consumes<cat::METCollection>           (iConfig.getParameter<edm::InputTag>("metLabel"));
  pvToken_       = consumes<int>                            (iConfig.getParameter<edm::InputTag>("pvLabel"));
  puWeight_      = consumes<float>                          (iConfig.getParameter<edm::InputTag>("puWeight"));

  jet_E  = new std::vector<float>;
  jet_pt  = new std::vector<float>;

  b_Jet_partonFlavour = new std::vector<int>;
  b_Jet_hadronFlavour = new std::vector<int>;

  b_Jet_smearedRes     = new std::vector<float>;
  b_Jet_smearedResDown = new std::vector<float>;
  b_Jet_smearedResUp   = new std::vector<float>;
  b_Jet_shiftedEnUp    = new std::vector<float>;
  b_Jet_shiftedEnDown  = new std::vector<float>;

  b_Jet_CSV = new std::vector<double>;
  csvid = new std::vector<int>;

  usesResource("TFileService");
  edm::Service<TFileService> fs;
  tree = fs->make<TTree>("tree", "TopTree");

  tree->Branch("EVENT",      &EVENT,       "EVENT/I");
  tree->Branch("RUN",        &RUN,         "RUN/I");
  tree->Branch("LUMI", &LUMI, "LUMI/I");

  tree->Branch("puweight", &puweight, "puweight/F");
  tree->Branch("nvertex",   &nvertex,  "nvertex/I");

  tree->Branch("MET",     &MET,     "MET/F");
  tree->Branch("MET_phi", &MET_phi, "MET_phi/F");

  tree->Branch("lep_E" , &lep_E,  "lep_E/F" );
  tree->Branch("lep_pt" , &lep_pt,  "lep_pt/F" );
  tree->Branch("lep_relIso" , &lep_relIso,  "lep_relIso/F" );
  tree->Branch("nLep", &nLep, "nLep/i");
  tree->Branch("nLooseLep", &nLooseLep, "nLooseLep/i");
  tree->Branch("nMuon", &nMuon, "nMuon/i");
  tree->Branch("nElec", &nElec, "nElec/i");

  tree->Branch("mt", &mt, "mt/F");

  tree->Branch("jet_E" , "std::vector<float>", &jet_E );
  tree->Branch("jet_pt" , "std::vector<float>", &jet_pt );

  tree->Branch("jet_partonFlavour", "std::vector<int>", &b_Jet_partonFlavour);
  tree->Branch("jet_hadronFlavour", "std::vector<int>", &b_Jet_hadronFlavour);

  tree->Branch("jet_smearedRes",     "std::vector<float>", &b_Jet_smearedRes);
  tree->Branch("jet_smearedResDown", "std::vector<float>", &b_Jet_smearedResDown);
  tree->Branch("jet_smearedResUp",   "std::vector<float>", &b_Jet_smearedResUp);
  tree->Branch("jet_shiftedEnUp",    "std::vector<float>", &b_Jet_shiftedEnUp);
  tree->Branch("jet_shiftedEnDown",  "std::vector<float>", &b_Jet_shiftedEnDown);

  tree->Branch("jet_CSV" , "std::vector<double>", &b_Jet_CSV );
  tree->Branch("csvid","std::vector<int>", &csvid);
  tree->Branch("nJets30", &nJets30, "nJets30/i");
  tree->Branch("nbJets30", &nbJets30, "nbJets30/i");

  tree->Branch("dileptonic", &dileptonic, "dileptonic/i"); //for test
  tree->Branch("semileptonic", &semileptonic, "semileptonic/i");
  tree->Branch("visible_ttjj", &visible_ttjj, "visible_ttjj/i"); //visible for semileptonic
  tree->Branch("visible_ttbb", &visible_ttbb, "visible_ttbb/i"); //visible for semileptonic
  tree->Branch("visible_ttbj", &visible_ttbj, "visible_ttbj/i"); //visible for semileptonic
  tree->Branch("visible_ttcc", &visible_ttcc, "visible_ttcc/i"); //visible for semileptonic
  tree->Branch("visible_ttLF", &visible_ttLF, "visible_ttLF/i"); //visible for semileptonic
  tree->Branch("ttjj", &ttjj, "ttjj/i");
  tree->Branch("ttbb", &ttbb, "ttbb/i");
}


TtbarSingleLeptonAnalyzer::~TtbarSingleLeptonAnalyzer()
{
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

  delete jet_E;
  delete jet_pt;

  delete b_Jet_partonFlavour;
  delete b_Jet_hadronFlavour;

  delete b_Jet_smearedRes;
  delete b_Jet_smearedResDown;
  delete b_Jet_smearedResUp;
  delete b_Jet_shiftedEnUp;
  delete b_Jet_shiftedEnDown;

  delete b_Jet_CSV;
  delete csvid;
}

//
// member functions
//

// ------------ method called for each event  ------------
void TtbarSingleLeptonAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  using namespace edm;

  puweight = 1;
  nvertex = 0;
  MET = -1.0;
  MET_phi = -1.0;
  lep_E = -1.0; 
  lep_pt = -1.0;
  lep_relIso = -1.0;
  nLep = -1;
  nLooseLep = -1;
  nMuon = -1;
  nElec = -1;
  mt = -1.0;

  nJets30 = -1;
  nbJets30 = -1;
  dileptonic = -9;
  semileptonic = -9;
  visible_ttjj = -9;
  visible_ttbb = -9;
  visible_ttbj = -9;
  visible_ttcc = -9;
  visible_ttLF = -9;
  ttjj = -9;
  ttbb = -9;

  jet_E->clear();
  jet_pt->clear();

  b_Jet_partonFlavour->clear();
  b_Jet_hadronFlavour->clear();

  b_Jet_smearedRes->clear();
  b_Jet_smearedResDown->clear();
  b_Jet_smearedResUp->clear();
  b_Jet_shiftedEnUp->clear();
  b_Jet_shiftedEnDown->clear();

  b_Jet_CSV->clear();
  csvid->clear();

  //---------------------------------------------------------------------------
  //---------------------------------------------------------------------------
  // Event Info
  //---------------------------------------------------------------------------
  //---------------------------------------------------------------------------

  EVENT        = iEvent.id().event();
  RUN          = iEvent.id().run();
  LUMI         = iEvent.luminosityBlock();

  //---------------------------------------------------------------------------
  //---------------------------------------------------------------------------
  // PU Info
  //---------------------------------------------------------------------------
  //---------------------------------------------------------------------------
  if(!iEvent.isRealData()) {
    edm::Handle<float> PUWeight;

    iEvent.getByToken(puWeight_, PUWeight);
    puweight = *PUWeight;
  }
  else puweight = 1;

  //event categorization
  edm::Handle<cat::GenTopCollection> genTops;
  iEvent.getByToken(genTopToken_, genTops);

  // init
  dileptonic = 0;
  semileptonic = 0;
  visible_ttjj = 0;
  visible_ttbb = 0;
  visible_ttbj = 0;
  visible_ttcc = 0;
  visible_ttLF = 0;
  ttjj = 0;
  ttbb = 0;

  cat::GenTop catGenTop = genTops->at(0);

  dileptonic = catGenTop.diLeptonic(0);
  semileptonic = catGenTop.semiLeptonic(0);

  visible_ttjj = catGenTop.NbJets20() >= 2 && catGenTop.NJets20() >= 6 && (( catGenTop.lepton1().pt() > 30 && abs(catGenTop.lepton1().eta()) < 2.4) || (catGenTop.lepton2().pt() > 30 && abs(catGenTop.lepton2().eta()) < 2.4)) && catGenTop.semiLeptonic(0);
  visible_ttbb = catGenTop.NbJets20() >= 4 && visible_ttjj;
  visible_ttbj = catGenTop.NbJets20() == 3 && visible_ttjj; 
  visible_ttcc = catGenTop.NcJets20() >= 2 && visible_ttjj;
  visible_ttLF = visible_ttjj && !visible_ttbb && !visible_ttbj && !visible_ttcc;

  ttjj = catGenTop.NaddJets20() >= 2;
  ttbb = catGenTop.NaddbJets20() >= 2;

  //---------------------------------------------------------------------------
  //---------------------------------------------------------------------------
  // Primary Vertex Info
  //---------------------------------------------------------------------------
  //---------------------------------------------------------------------------

  edm::Handle<int> pvHandle;
  iEvent.getByToken( pvToken_, pvHandle );

  nvertex = *pvHandle;

  //---------------------------------------------------------------------------
  //---------------------------------------------------------------------------
  // Missing E_T
  //---------------------------------------------------------------------------
  //---------------------------------------------------------------------------

  Handle<cat::METCollection> METHandle;
  iEvent.getByToken(metToken_, METHandle);

  // MET-PF
  MET     = METHandle->begin()->pt();
  MET_phi = METHandle->begin()->phi();

  //---------------------------------------------------------------------------
  //---------------------------------------------------------------------------
  // Electrons
  //---------------------------------------------------------------------------
  //---------------------------------------------------------------------------

  cat::ElectronCollection selectedElectrons;
  cat::ElectronCollection vetoElectrons;

  Handle<cat::ElectronCollection> electrons;
  iEvent.getByToken(electronToken_, electrons);

  for (unsigned int i = 0; i < electrons->size() ; i++) {
    const cat::Electron & electron = electrons->at(i);
    if( IsTightElectron( electron ) ) selectedElectrons.push_back( electron );
    else if( IsLooseElectron( electron ) ) vetoElectrons.push_back( electron ); // does not Include selected electrons

  }

  //---------------------------------------------------------------------------
  //---------------------------------------------------------------------------
  // Muons
  //---------------------------------------------------------------------------
  //---------------------------------------------------------------------------

  cat::MuonCollection selectedMuons;
  cat::MuonCollection vetoMuons;

  Handle<cat::MuonCollection> muons;
  iEvent.getByToken(muonToken_, muons);

  for (unsigned int i = 0; i < muons->size() ; i++) {
    const cat::Muon & muon = muons->at(i);
    bool pass = muon.pt() > 30 && fabs(muon.eta()) < 2.4 && muon.isTightMuon() && muon.relIso() < 0.15;
    if( pass ) selectedMuons.push_back( muon);
    else if( muon.isLooseMuon() ) vetoMuons.push_back( muon); // does not Include selected muons
  }

  //---------------------------------------------------------------------------
  //---------------------------------------------------------------------------
  // Channel Selection
  //---------------------------------------------------------------------------
  //---------------------------------------------------------------------------

  nElec = selectedElectrons.size();
  nMuon = selectedMuons.size();
  nLep = nElec + nMuon;
  nLooseLep = vetoMuons.size() + vetoElectrons.size() ;

  TLorentzVector lepton;

  if(selectedMuons.size() == 1 && selectedElectrons.size() == 0 && nLooseLep == 0){
    lep_E = selectedMuons[0].energy();
    lep_pt = selectedMuons[0].pt();
    lep_relIso = selectedMuons[0].relIso();
    lepton.SetPxPyPzE(selectedMuons[0].px(), selectedMuons[0].py(), selectedMuons[0].pz(), selectedMuons[0].energy());
    mt = transverseMass( selectedMuons[0].p4(), METHandle->begin()->p4() ); 
  }

  if(selectedMuons.size() == 0 && selectedElectrons.size() == 1 && nLooseLep == 0){
    lep_E = selectedElectrons[0].energy();
    lep_pt = selectedElectrons[0].pt();
    lep_relIso = selectedElectrons[0].relIso();
    lepton.SetPxPyPzE(selectedElectrons[0].px(), selectedElectrons[0].py(), selectedElectrons[0].pz(), selectedElectrons[0].energy());
  }

  //---------------------------------------------------------------------------
  //---------------------------------------------------------------------------
  // Jets
  //---------------------------------------------------------------------------
  //---------------------------------------------------------------------------

  Handle<cat::JetCollection> jets;
  iEvent.getByToken(jetToken_, jets);

  std::map<int,double> mapJetBDiscriminator;

  nJets30 = 0;
  nbJets30 = 0;
  int index = 0;

  for (unsigned int i = 0; i < jets->size() ; i++) {

    const cat::Jet & jet = jets->at(i);

    bool goodJet  = false;
    bool cleanJet = false;

    // Jet Selection
    if(std::abs(jet.eta()) < 2.4 && jet.pt() > 30 && jet.LooseId()) goodJet = true;
    // Jet Cleaning
    TLorentzVector vjet(jet.px(), jet.py(), jet.pz(), jet.energy());
    double dr_LepJet = vjet.DeltaR(lepton);
    if(dr_LepJet > 0.4) cleanJet = true;

    if(goodJet && cleanJet){
      nJets30++;
      // Basic variables
      jet_E ->push_back(jet.energy());
      jet_pt ->push_back(jet.pt());

      // Parton Flavour
      b_Jet_partonFlavour->push_back(jet.partonFlavour());
      b_Jet_hadronFlavour->push_back(jet.hadronFlavour());

      // Smeared and Shifted
      b_Jet_smearedRes     ->push_back(jet.smearedRes() );
      b_Jet_smearedResDown ->push_back(jet.smearedResDown());
      b_Jet_smearedResUp   ->push_back(jet.smearedResUp());
      b_Jet_shiftedEnUp    ->push_back(jet.shiftedEnUp());
      b_Jet_shiftedEnDown  ->push_back(jet.shiftedEnDown());

      // b-tag discriminant
      double bDiscriminator = jet.bDiscriminator(BTAG_CSVv2);
      mapJetBDiscriminator[index] = bDiscriminator;
      b_Jet_CSV->push_back(bDiscriminator);
      if( bDiscriminator > WP_BTAG_CSVv2M) nbJets30++;
      // BTagSFUtil *fBTagSF;   //The BTag SF utility
      // fBTagSF = new BTagSFUtil("CSVv2", "Medium", 0);

      // if (fBTagSF->IsTagged(jet_btagDis_CSV, -999999, jet.pt(), jet.eta()) ){
      //    std::cout << "Is btag Jet....." << std::endl;
      // }
      index++;
    }
  }

  //csv order

  std::vector< std::pair<int,double> > vecJetBDisc(mapJetBDiscriminator.begin(), mapJetBDiscriminator.end());
  std::sort(vecJetBDisc.begin(), vecJetBDisc.end(), bigger_second<data_t>());
  for( std::vector< std::pair<int,double> >::iterator it = vecJetBDisc.begin() ; it != vecJetBDisc.end(); ++it)
      csvid->push_back((*it).first);


  // Fill Tree with event at 1 lepton cut level
  if( nLep > 0 ) tree->Fill();

}

//------------- Good Electron Selection -----------------------
bool TtbarSingleLeptonAnalyzer::IsTightElectron(const cat::Electron & i_electron_candidate)
{
  bool GoodElectron=true;

  GoodElectron &= (i_electron_candidate.isPF() );            // PF
  GoodElectron &= (i_electron_candidate.pt() > 30);          // pT
  GoodElectron &= (std::abs(i_electron_candidate.eta()) < 2.4);  // eta
  GoodElectron &= (std::abs(i_electron_candidate.scEta()) < 1.4442 || // eta Super-Cluster
                   std::abs(i_electron_candidate.scEta()) > 1.566);

  // From https://twiki.cern.ch/twiki/bin/viewauth/CMS/CutBasedElectronIdentificationRun2
  GoodElectron &= i_electron_candidate.electronID("cutBasedElectronID-Spring15-25ns-V1-standalone-loose") > 0.0;

  //----------------------------------------------------------------------------------------------------
  //------------- The Relative Isolation is already calculated in the CAT object -----------------------
  //----------------------------------------------------------------------------------------------------
  // relIso( R ) already includes AEff and RhoIso
  // float relIso = ( chIso + std::max(0.0, nhIso + phIso - rhoIso*AEff) )/ ecalpt;

  GoodElectron &=( i_electron_candidate.relIso( 0.3 ) < 0.12 );

  // Effective Area Parametrization can be found in:
  // Last recommendation: https://twiki.cern.ch/twiki/bin/viewauth/CMS/SWGuideMuonId2015 Slide 8
  //----------------------------------------------------------------------------------------------------
  //----------------------------------------------------------------------------------------------------

  return GoodElectron;

}
//------------- Loose Electron Selection -----------------------
bool TtbarSingleLeptonAnalyzer::IsLooseElectron(const cat::Electron & i_electron_candidate)
{
  bool GoodElectron=true;

  GoodElectron &= (i_electron_candidate.isPF() );            // PF
  GoodElectron &= (i_electron_candidate.pt() > 15);          // pT
  GoodElectron &= (std::abs(i_electron_candidate.eta()) < 2.4);  // eta
  GoodElectron &= (std::abs(i_electron_candidate.scEta()) < 1.4442 || // eta Super-Cluster
                   std::abs(i_electron_candidate.scEta()) > 1.566);

  // From https://twiki.cern.ch/twiki/bin/viewauth/CMS/CutBasedElectronIdentificationRun2
  GoodElectron &= i_electron_candidate.electronID("cutBasedElectronID-Spring15-25ns-V1-standalone-medium") > 0.0;

  //----------------------------------------------------------------------------------------------------
  //------------- The Relative Isolation is already calculated in the CAT object -----------------------
  //----------------------------------------------------------------------------------------------------
  GoodElectron &=( i_electron_candidate.relIso( 0.3 ) < 0.12 );
  //----------------------------------------------------------------------------------------------------
  //----------------------------------------------------------------------------------------------------

  return GoodElectron;

}

double TtbarSingleLeptonAnalyzer::transverseMass( const reco::Candidate::LorentzVector& lepton, 
				     const reco::Candidate::LorentzVector& met) {  
  reco::Candidate::LorentzVector leptonT(lepton.Px(),lepton.Py(),0.,lepton.E()*sin(lepton.Theta()));
  reco::Candidate::LorentzVector sumT=leptonT+met;
  return std::sqrt(sumT.M2());
}

//define this as a plug-in
DEFINE_FWK_MODULE(TtbarSingleLeptonAnalyzer);
