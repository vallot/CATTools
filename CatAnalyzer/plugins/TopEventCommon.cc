#include "CATTools/CatAnalyzer/interface/TopEventCommon.h"

#include "FWCore/Framework/interface/LuminosityBlock.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "CATTools/CatAnalyzer/interface/TopTriggerSF.h"
#include "CATTools/CatAnalyzer/interface/TopEventGlobalVar.h"
#include "CATTools/CatAnalyzer/interface/KinematicReconstructionSolution.h"

#include "DataFormats/Math/interface/deltaR.h"
#include "DataFormats/Math/interface/deltaPhi.h"
#include "CATTools/CatAnalyzer/interface/TTDiLeptonEventSelector.h"

using namespace std;
using namespace cat;
using namespace TopEventCommonGlobal;

void TopEventCommon::paramInit(const edm::ParameterSet& iConfig) {
  partonTop_channel_ = consumes<int>(iConfig.getParameter<edm::InputTag>("partonTop_channel"));
  topPtWeight_ = consumes<float>(iConfig.getParameter<edm::InputTag>("topPtWeight"));

  puweightToken_ = consumes<float>(iConfig.getParameter<edm::InputTag>("puweight"));
  puweightToken_up_ = consumes<float>(iConfig.getParameter<edm::InputTag>("puweight_up"));
  puweightToken_dn_ = consumes<float>(iConfig.getParameter<edm::InputTag>("puweight_dn"));

  partonTop_modes_   = consumes<vector<int> >(iConfig.getParameter<edm::InputTag>("partonTop_modes"));
  partonTop_genParticles_   = consumes<reco::GenParticleCollection>(iConfig.getParameter<edm::InputTag>("partonTop_genParticles"));

  pseudoTop_leptons_   = consumes<edm::View<reco::Candidate> >(edm::InputTag("pseudoTop", "leptons"));
  pseudoTop_neutrinos_   = consumes<edm::View<reco::Candidate> >(edm::InputTag("pseudoTop", "neutrinos"));
  pseudoTop_jets_ = consumes<edm::View<reco::Candidate> >(edm::InputTag("pseudoTop", "jets"));

  genWeightToken_ = consumes<float>(iConfig.getParameter<edm::InputTag>("genweight"));
  pdfweightToken_ = consumes<vector<float>>(iConfig.getParameter<edm::InputTag>("pdfweight"));
  scaleupweightsToken_ = consumes<vector<float>>(iConfig.getParameter<edm::InputTag>("scaleupweight"));
  scaledownweightsToken_ = consumes<vector<float>>(iConfig.getParameter<edm::InputTag>("scaledownweight"));

  TopEventInfo& evInfo_ = TopEventInfo::getInstance();

  usesResource("TFileService");
  edm::Service<TFileService> fs;

  h_nevents = fs->make<TH1D>("nevents","nevents",1,0,1);
  for (int sys = 0; sys < nsys_e; ++sys){
    ttree_.push_back(fs->make<TTree>(sys_name[sys].c_str(), sys_name[sys].c_str()));
    auto tr = ttree_.back();
    setBranch(tr, sys);
  }
  for (int i = 0; i < NCutflow; i++) evInfo_.cutflow_.push_back({0,0,0,0});

  evInfo_.kinematicReconstruction.reset(new KinematicReconstruction(1, true));
}

TopEventCommon::TopEventCommon(const edm::ParameterSet& iConfig ):iConfig_(iConfig){
  paramInit(iConfig);
}

TopEventCommon::~TopEventCommon(){}

void TopEventCommon::beginJob(){
}

void TopEventCommon::endJob(){
  showSummary();
}

void TopEventCommon::showSummary() {
  TopEventInfo& evInfo_ = TopEventInfo::getInstance();
  cout <<setw(10)<<"cut flow"<<setw(10)<<"no ll"<<setw(10)<<"emu"<<setw(10)<<"ee"<<setw(10)<<"mumu"<< endl;
  for ( int i=0; i<NCutflow; ++i ) {
    cout <<setw(10)<<"step "<<i<< setw(10)<<evInfo_.cutflow_[i][0] << setw(10)<<evInfo_.cutflow_[i][1] << setw(10)<<evInfo_.cutflow_[i][2] <<setw(10)<<evInfo_.cutflow_[i][3]<< endl;
  }
}

void TopEventCommon::setBranch(TTree* tr, int sys) {
  TopEventInfo& evInfo_ = TopEventInfo::getInstance();
  tr->Branch("run",     &(evInfo_.run)    , "run/I");
  tr->Branch("event",   &(evInfo_.event)  , "event/I");

  tr->Branch("nvertex", &(evInfo_.nvertex), "nvertex/I");
  tr->Branch("step",    &(evInfo_.step)   , "step/I");
  tr->Branch("channel", &(evInfo_.channel), "channel/I");
  tr->Branch("njet",    &(evInfo_.njet), "njet/I");
  tr->Branch("nbjet",   &(evInfo_.nbjet), "nbjet/I");
  tr->Branch("step1",   &(evInfo_.step1), "step1/O");
  tr->Branch("step2",   &(evInfo_.step2), "step2/O");
  tr->Branch("step3", &(evInfo_.step3), "step3/O");
  tr->Branch("step4", &(evInfo_.step4), "step4/O");
  tr->Branch("step5", &(evInfo_.step5), "step5/O");
  tr->Branch("step6", &(evInfo_.step6), "step6/O");
  tr->Branch("step7", &(evInfo_.step7), "step7/O");
  tr->Branch("step8", &(evInfo_.step8), "step8/O");
  tr->Branch("tri", &(evInfo_.tri), "tri/F");
  tr->Branch("tri_up", &(evInfo_.tri_up), "tri_up/F");
  tr->Branch("tri_dn", &(evInfo_.tri_dn), "tri_dn/F");
  tr->Branch("filtered", &(evInfo_.filtered), "filtered/O");
  tr->Branch("met", &(evInfo_.met), "met/F");
  tr->Branch("weight", &(evInfo_.weight), "weight/F");
  tr->Branch("topPtWeight", &(evInfo_.topPtWeight), "topPtWeight/F");
  tr->Branch("puweight", &(evInfo_.puweight), "puweight/F");
  tr->Branch("puweight_up", &(evInfo_.puweight_up), "puweight_up/F");
  tr->Branch("puweight_dn", &(evInfo_.puweight_dn), "puweight_dn/F");
  tr->Branch("genweight", &(evInfo_.genweight), "genweight/F");
  tr->Branch("mueffweight", &(evInfo_.mueffweight), "mueffweight/F");
  tr->Branch("mueffweight_up", &(evInfo_.mueffweight_up), "mueffweight_up/F");
  tr->Branch("mueffweight_dn", &(evInfo_.mueffweight_dn), "mueffweight_dn/F");
  tr->Branch("eleffweight", &(evInfo_.eleffweight), "eleffweight/F");
  tr->Branch("eleffweight_up", &(evInfo_.eleffweight_up), "eleffweight_up/F");
  tr->Branch("eleffweight_dn", &(evInfo_.eleffweight_dn), "eleffweight_dn/F");
  tr->Branch("btagweight", &(evInfo_.btagweight), "btagweight/F");
  tr->Branch("btagweight_up", &(evInfo_.btagweight_up), "btagweight_up/F");
  tr->Branch("btagweight_dn", &(evInfo_.btagweight_dn), "btagweight_dn/F");

  if (sys == 0){
    tr->Branch("csvweights","std::vector<float>",&(evInfo_.csvweights));
    tr->Branch("pdfWeights","std::vector<float>",&(evInfo_.pdfWeights));
    tr->Branch("scaleWeights_up","std::vector<float>",&(evInfo_.scaleWeights_up));
    tr->Branch("scaleWeights_dn","std::vector<float>",&(evInfo_.scaleWeights_dn));
  }

  tr->Branch("gen_partonChannel", &(evInfo_.gen_partonChannel), "gen_partonChannel/I");
  tr->Branch("gen_partonMode1", &(evInfo_.gen_partonMode1), "gen_partonMode1/I");
  tr->Branch("gen_partonMode2", &(evInfo_.gen_partonMode2), "gen_partonMode2/I");
  tr->Branch("gen_partonMode", &(evInfo_.gen_partonMode), "gen_partonMode/I");
  tr->Branch("gen_partonInPhase", &(evInfo_.gen_partonInPhase), "gen_partonInPhase/O");
  tr->Branch("gen_partonInPhaseLep", &(evInfo_.gen_partonInPhaseLep), "gen_partonInPhaseLep/O");
  tr->Branch("gen_partonInPhaseJet", &(evInfo_.gen_partonInPhaseJet), "gen_partonInPhaseJet/O");
  tr->Branch("gen_partonlep1_pid", &(evInfo_.gen_partonlep1_pid), "gen_partonlep1_pid/I");
  tr->Branch("gen_partonlep2_pid", &(evInfo_.gen_partonlep2_pid), "gen_partonlep2_pid/I");

  if (sys == 0){
    tr->Branch("gen_partonlep1", "TLorentzVector", &(evInfo_.gen_partonlep1));
    tr->Branch("gen_partonlep2", "TLorentzVector", &(evInfo_.gen_partonlep2));
    tr->Branch("gen_partondilep", "TLorentzVector", &(evInfo_.gen_partondilep));
    tr->Branch("gen_partonjet1", "TLorentzVector", &(evInfo_.gen_partonjet1));
    tr->Branch("gen_partonjet2", "TLorentzVector", &(evInfo_.gen_partonjet2));
    tr->Branch("gen_partontop1", "TLorentzVector", &(evInfo_.gen_partontop1));
    tr->Branch("gen_partontop2", "TLorentzVector", &(evInfo_.gen_partontop2));
    tr->Branch("gen_partonnu1", "TLorentzVector", &(evInfo_.gen_partonnu1));
    tr->Branch("gen_partonnu2", "TLorentzVector", &(evInfo_.gen_partonnu2));
    tr->Branch("gen_partonttbar", "TLorentzVector", &(evInfo_.gen_partonttbar));
    tr->Branch("gen_partonttbar_dphi", &(evInfo_.gen_partonttbar_dphi), "gen_partonttbar_dphi/F");
  }

  tr->Branch("gen_pseudoChannel", &(evInfo_.gen_pseudoChannel), "gen_pseudoChannel/I");
  tr->Branch("gen_pseudoInPhase", &(evInfo_.gen_pseudoInPhase), "gen_pseudoInPhase/O");
  tr->Branch("gen_pseudolep1_pid", &(evInfo_.gen_pseudolep1_pid), "gen_pseudolep1_pid/I");
  tr->Branch("gen_pseudolep2_pid", &(evInfo_.gen_pseudolep2_pid), "gen_pseudolep2_pid/I");

  if (sys == 0){
    tr->Branch("gen_pseudolep1", "TLorentzVector", &(evInfo_.gen_pseudolep1));
    tr->Branch("gen_pseudolep2", "TLorentzVector", &(evInfo_.gen_pseudolep2));
    tr->Branch("gen_pseudodilep", "TLorentzVector", &(evInfo_.gen_pseudodilep));
    tr->Branch("gen_pseudojet1", "TLorentzVector", &(evInfo_.gen_pseudojet1));
    tr->Branch("gen_pseudojet2", "TLorentzVector", &(evInfo_.gen_pseudojet2));
    tr->Branch("gen_pseudotop1", "TLorentzVector", &(evInfo_.gen_pseudotop1));
    tr->Branch("gen_pseudotop2", "TLorentzVector", &(evInfo_.gen_pseudotop2));
    tr->Branch("gen_pseudonu1", "TLorentzVector", &(evInfo_.gen_pseudonu1));
    tr->Branch("gen_pseudonu2", "TLorentzVector", &(evInfo_.gen_pseudonu2));
    tr->Branch("gen_pseudottbar", "TLorentzVector", &(evInfo_.gen_pseudottbar));
    tr->Branch("gen_pseudottbar_dphi", &(evInfo_.gen_pseudottbar_dphi), "gen_pseudottbar_dphi/F");
  }

  tr->Branch("lep1", "TLorentzVector", &(evInfo_.lep1));
  tr->Branch("lep1_pid", &(evInfo_.lep1_pid), "lep1_pid/I");
  tr->Branch("lep2", "TLorentzVector", &(evInfo_.lep2));
  tr->Branch("lep2_pid", &(evInfo_.lep2_pid), "lep2_pid/I");
  tr->Branch("dilep", "TLorentzVector", &(evInfo_.dilep));

  tr->Branch("pseudojet1", "TLorentzVector", &(evInfo_.pseudojet1));
  tr->Branch("pseudojet1_CSVInclV2", &(evInfo_.pseudojet1_CSVInclV2), "pseudojet1_CSVInclV2/F");
  tr->Branch("pseudojet2", "TLorentzVector", &(evInfo_.pseudojet2));
  tr->Branch("pseudojet2_CSVInclV2", &(evInfo_.pseudojet2_CSVInclV2), "pseudojet2_CSVInclV2/F");
  tr->Branch("pseudotop1", "TLorentzVector", &(evInfo_.pseudotop1));
  tr->Branch("pseudotop2", "TLorentzVector", &(evInfo_.pseudotop2));
  tr->Branch("pseudonu1", "TLorentzVector", &(evInfo_.pseudonu1));
  tr->Branch("pseudonu2", "TLorentzVector", &(evInfo_.pseudonu2));
  tr->Branch("pseudottbar", "TLorentzVector", &(evInfo_.pseudottbar));
  tr->Branch("pseudottbar_dphi", &(evInfo_.pseudottbar_dphi), "pseudottbar_dphi/F");

  tr->Branch("partonjet1", "TLorentzVector", &(evInfo_.partonjet1));
  tr->Branch("partonjet1_CSVInclV2", &(evInfo_.partonjet1_CSVInclV2), "partonjet1_CSVInclV2/F");
  tr->Branch("partonjet2", "TLorentzVector", &(evInfo_.partonjet2));
  tr->Branch("partonjet2_CSVInclV2", &(evInfo_.partonjet2_CSVInclV2), "partonjet2_CSVInclV2/F");
  tr->Branch("partontop1", "TLorentzVector", &(evInfo_.partontop1));
  tr->Branch("partontop2", "TLorentzVector", &(evInfo_.partontop2));
  tr->Branch("partonnu1", "TLorentzVector", &(evInfo_.partonnu1));
  tr->Branch("partonnu2", "TLorentzVector", &(evInfo_.partonnu2));
  tr->Branch("partonttbar", "TLorentzVector", &(evInfo_.partonttbar));
  tr->Branch("partonttbar_dphi", &(evInfo_.partonttbar_dphi), "partonttbar_dphi/F");

  /*
  tr->Branch("is3lep", &(evInfo_.is3lep), "is3lep/I");
  tr->Branch("desyjet1", "TLorentzVector", &(evInfo_.desyjet1));
  tr->Branch("desyjet1_CSVInclV2", &(evInfo_.desyjet1_CSVInclV2), "desyjet1_CSVInclV2/F");
  tr->Branch("desyjet2", "TLorentzVector", &(evInfo_.desyjet2));
  tr->Branch("desyjet2_CSVInclV2", &(evInfo_.desyjet2_CSVInclV2), "desyjet2_CSVInclV2/F");
  tr->Branch("desytop1", "TLorentzVector", &(evInfo_.desytop1));
  tr->Branch("desytop2", "TLorentzVector", &(evInfo_.desytop2));
  tr->Branch("desyttbar", "TLorentzVector", &(evInfo_.desyttbar));
  tr->Branch("desyttbar_dphi", &(evInfo_.desyttbar_dphi), "desyttbar_dphi/F");
  */



}

void TopEventCommon::beginLuminosityBlock(const edm::LuminosityBlock& lumi, const edm::EventSetup&)
{
  /*
  if ( dynamic_cast<DESYSmearedSolver*>(evInfo_.solver_.get()) != 0 ) {
    edm::Service<edm::RandomNumberGenerator> rng;
    CLHEP::HepRandomEngine& engine = rng->getEngine(lumi.index());
    dynamic_cast<DESYSmearedSolver*>(evInfo_.solver_.get())->setRandom(&engine);
  }
  if ( dynamic_cast<DESYSmearedSolver*>(evInfo_.solverPT_.get()) != 0 ) {
    edm::Service<edm::RandomNumberGenerator> rngPT;
    CLHEP::HepRandomEngine& enginePT = rngPT->getEngine(lumi.index());
    dynamic_cast<DESYSmearedSolver*>(evInfo_.solverPT_.get())->setRandom(&enginePT);
  }
  */
}

void TopEventCommon::genInfo(const edm::Event& iEvent, const edm::EventSetup& iSetup){
  TopEventInfo& evInfo_ = TopEventInfo::getInstance();
  evInfo_.keepTtbarSignal = false;
  const bool runOnMC = !iEvent.isRealData();

  if (runOnMC){
    edm::Handle<float> puweightHandle;
    iEvent.getByToken(puweightToken_, puweightHandle);
    evInfo_.puweight = *puweightHandle;

    edm::Handle<float> puweightHandle_up;
    iEvent.getByToken(puweightToken_up_, puweightHandle_up);
    evInfo_.puweight_up = *puweightHandle_up;

    edm::Handle<float> puweightHandle_dn;
    iEvent.getByToken(puweightToken_dn_, puweightHandle_dn);
    evInfo_.puweight_dn = *puweightHandle_dn;

    edm::Handle<float> genweightHandle;
    iEvent.getByToken(genWeightToken_, genweightHandle);
    evInfo_.genweight = (*genweightHandle);
    evInfo_.weight = evInfo_.genweight*evInfo_.puweight;

    // edm::Handle<cat::GenWeights> genweightHandle;
    // iEvent.getByToken(genweightsToken_, genweightHandle);
    // evInfo_.genweight = genweightHandle->genWeight();
    // evInfo_.weight = evInfo_.genweight*evInfo_.puweight;
  }

  h_nevents->Fill(0.5,evInfo_.puweight*evInfo_.genweight);


  edm::Handle<int> partonTop_channel;
  if ( iEvent.getByToken(partonTop_channel_, partonTop_channel)){

    edm::Handle<float> topPtWeightHandle;
    iEvent.getByToken(topPtWeight_, topPtWeightHandle);
    evInfo_.topPtWeight = *topPtWeightHandle;

    edm::Handle<vector<float>> pdfweightHandle;
    iEvent.getByToken(pdfweightToken_, pdfweightHandle);
    for (const float & aPdfWeight : *pdfweightHandle){
      evInfo_.pdfWeights.push_back(aPdfWeight);
    }
    edm::Handle<vector<float>> scaleupweightsHandle, scaledownweightsHandle;
    iEvent.getByToken(scaleupweightsToken_, scaleupweightsHandle);
    for (const float & aScaleWeight : *scaleupweightsHandle){
      evInfo_.scaleWeights_up.push_back(aScaleWeight);
    }
    iEvent.getByToken(scaledownweightsToken_, scaledownweightsHandle);
    for (const float & aScaleWeight : *scaledownweightsHandle){
      evInfo_.scaleWeights_dn.push_back(aScaleWeight);
    }

    edm::Handle<vector<int> > partonTop_modes;
    edm::Handle<reco::GenParticleCollection> partonTop_genParticles;
    iEvent.getByToken(partonTop_modes_, partonTop_modes);
    iEvent.getByToken(partonTop_genParticles_, partonTop_genParticles);
    if ( (*partonTop_modes).size() == 0 ) {
      evInfo_.gen_partonMode1 = 0;
      evInfo_.gen_partonMode2 = 0;
    }
    else if ( (*partonTop_modes).size() == 1 ) { evInfo_.gen_partonMode2 = 0; }
    else{
      evInfo_.gen_partonChannel = *partonTop_channel;
      evInfo_.gen_partonMode1 = (*partonTop_modes)[0];
      evInfo_.gen_partonMode2 = (*partonTop_modes)[1];
    }
    if (evInfo_.gen_partonChannel == CH_FULLLEPTON) evInfo_.keepTtbarSignal = true;
    if (evInfo_.gen_partonChannel == CH_SEMILEPTON) evInfo_.keepTtbarSignal = true;

    if(evInfo_.gen_partonMode1==1 && evInfo_.gen_partonMode2==2) evInfo_.gen_partonMode=1;
    if(evInfo_.gen_partonMode1==2 && evInfo_.gen_partonMode2==1) evInfo_.gen_partonMode=1;
    if(evInfo_.gen_partonMode1==1 && evInfo_.gen_partonMode2==1) evInfo_.gen_partonMode=3;
    if(evInfo_.gen_partonMode1==2 && evInfo_.gen_partonMode2==2) evInfo_.gen_partonMode=2;
    if(evInfo_.gen_partonMode1>3 || evInfo_.gen_partonMode2>3)   evInfo_.gen_partonMode=4;

    if ( !(partonTop_genParticles->empty()) ){

      // Get Top quark pairs
      auto parton1 = &partonTop_genParticles->at(0);
      auto parton2 = &partonTop_genParticles->at(1);
      if (parton1->charge() < 0) swap(parton1, parton2);

      evInfo_.gen_partontop1 = ToTLorentzVector(*parton1);
      evInfo_.gen_partontop2 = ToTLorentzVector(*parton2);
      evInfo_.gen_partonttbar = evInfo_.gen_partontop1 + evInfo_.gen_partontop2;
      evInfo_.gen_partonttbar_dphi = evInfo_.gen_partontop1.DeltaPhi(evInfo_.gen_partontop2);

      // Get W and b quarks
      if ( parton1 and parton2 ) {
        const auto partonW1 = parton1->daughter(0);
        const auto partonB1 = parton1->daughter(1);
        const auto partonW2 = parton2->daughter(0);
        const auto partonB2 = parton2->daughter(1);

        if ( (partonB1->pt() > 30 && std::abs(partonB1->eta()) < 2.4) &&
            (partonB2->pt() > 30 && std::abs(partonB2->eta()) < 2.4))
          evInfo_.gen_partonInPhaseJet = true;

        // Get W daughters
        if ( partonW1 and partonW2 and partonB1 and partonB2 ) {
          const auto partonW11 = partonW1->daughter(0);
          const auto partonW12 = partonW1->daughter(1);
          const auto partonW21 = partonW2->daughter(0);
          const auto partonW22 = partonW2->daughter(1);
          if ( (partonW11->pt() > 20 && std::abs(partonW11->eta()) < 2.4 && (std::abs(partonW11->pdgId()) == 11 || std::abs(partonW11->pdgId()) == 13) ) &&
              (partonW21->pt() > 20 && std::abs(partonW21->eta()) < 2.4 && (std::abs(partonW11->pdgId()) == 11 || std::abs(partonW11->pdgId()) == 13) ))
            evInfo_.gen_partonInPhaseLep = true;

          // Fill lepton informations
          evInfo_.gen_partonnu1 = ToTLorentzVector(*partonW12);
          evInfo_.gen_partonnu2 = ToTLorentzVector(*partonW22);
          evInfo_.gen_partonlep1 = ToTLorentzVector(*partonW11);
          evInfo_.gen_partonlep2 = ToTLorentzVector(*partonW21);
          evInfo_.gen_partonlep1_pid = partonW11->pdgId();
          evInfo_.gen_partonlep2_pid = partonW21->pdgId();
          evInfo_.gen_partondilep = evInfo_.gen_partonlep1 + evInfo_.gen_partonlep2;
          evInfo_.gen_partonjet1 = ToTLorentzVector(*partonB1);
          evInfo_.gen_partonjet2 = ToTLorentzVector(*partonB2);
        }
      }
      if (evInfo_.gen_partonInPhaseJet && evInfo_.gen_partonInPhaseLep) evInfo_.gen_partonInPhase = true;
    }

    // Start to build pseudo top
    evInfo_.gen_pseudoChannel = CH_NOLL;
    edm::Handle<edm::View<reco::Candidate> > pseudoTopLeptonHandle;
    edm::Handle<edm::View<reco::Candidate> > pseudoTopNeutrinoHandle;
    edm::Handle<edm::View<reco::Candidate> > pseudoTopJetHandle;
    iEvent.getByToken(pseudoTop_leptons_, pseudoTopLeptonHandle);
    iEvent.getByToken(pseudoTop_neutrinos_, pseudoTopNeutrinoHandle);
    iEvent.getByToken(pseudoTop_jets_, pseudoTopJetHandle);
    do {
      // Basic lepton, jet multiplicity
      if ( pseudoTopLeptonHandle->size() < 2 or pseudoTopJetHandle->size() < 2 or pseudoTopNeutrinoHandle->size() < 2 ) break;

      std::vector<size_t> leptonIdxs, neutrinoIdxs, bjetIdxs;
      // Lepton acceptance cuts
      for ( size_t i=0, n=pseudoTopLeptonHandle->size(); i<n; ++i ) {
        const auto& x = pseudoTopLeptonHandle->at(i);
        if ( x.pt() < 20 or std::abs(x.eta()) > 2.4 ) continue;
        if ( abs(x.pdgId()) != 11 and abs(x.pdgId()) != 13 ) continue;
        leptonIdxs.push_back(i);
      }
      if ( leptonIdxs.size() !=2 ) break;
      std::nth_element(leptonIdxs.begin(), leptonIdxs.begin()+2, leptonIdxs.end(),
          [&](size_t i, size_t j){return pseudoTopLeptonHandle->at(i).pt() > pseudoTopLeptonHandle->at(j).pt();});
      auto lepton1 = pseudoTopLeptonHandle->at(leptonIdxs[0]).p4();
      auto lepton2 = pseudoTopLeptonHandle->at(leptonIdxs[1]).p4();
      const int pseudoW1DauId = abs(pseudoTopLeptonHandle->at(leptonIdxs[0]).pdgId());
      const int pseudoW2DauId = abs(pseudoTopLeptonHandle->at(leptonIdxs[1]).pdgId());
      switch ( pseudoW1DauId+pseudoW2DauId ) {
        case 22: evInfo_.gen_pseudoChannel = CH_ELEL; break;
        case 26: evInfo_.gen_pseudoChannel = CH_MUMU; break;
        case 24: evInfo_.gen_pseudoChannel = CH_MUEL; break;
        default: evInfo_.gen_pseudoChannel = CH_NOLL;
      }
      if (evInfo_.gen_pseudoChannel > 0) evInfo_.keepTtbarSignal = true;
      //std::nth_element(neutrinoIdxs.begin(), neutrinoIdxs.begin()+2, neutrinoIdxs.end(),
      //                 [&](size_t i, size_t j){return pseudoTopLeptonHandle->at(i).pt() > pseudoTopLeptonHandle->at(j).pt();});
      auto nu1 = pseudoTopNeutrinoHandle->at(0).p4(), nu2 = pseudoTopNeutrinoHandle->at(1).p4();

      // Jet acceptance and generator level b tag
      for ( size_t i=0, n=pseudoTopJetHandle->size(); i<n; ++i ) {
        const auto& x = pseudoTopJetHandle->at(i);
        if ( x.pt() < 30 or std::abs(x.eta()) > 2.4 ) continue;
        if ( abs(x.pdgId()) != 5 ) continue;
        bjetIdxs.push_back(i);
      }
      if ( bjetIdxs.size() > 1 )evInfo_.gen_pseudoInPhase=true;
      if ( bjetIdxs.size() < 2 ){evInfo_.gen_pseudoChannel = CH_NOLL; break;}

      std::nth_element(bjetIdxs.begin(), bjetIdxs.begin()+2, bjetIdxs.end(),
          [&](size_t i, size_t j){return pseudoTopJetHandle->at(i).pt() > pseudoTopJetHandle->at(j).pt();});
      auto bjet1 = pseudoTopJetHandle->at(bjetIdxs[0]).p4(), bjet2 = pseudoTopJetHandle->at(bjetIdxs[1]).p4();

      // Do the W combinations
      auto w1 = lepton1 + nu1;
      auto w2 = lepton2 + nu2;
      if ( true ) {
        const auto w1Alt = lepton1 + nu2;
        const auto w2Alt = lepton2 + nu1;

        const double wMass = 80.4;
        const double dm = std::abs(w1.mass()-wMass)+std::abs(w2.mass()-wMass);
        const double dmAlt = std::abs(w1Alt.mass()-wMass)+std::abs(w2Alt.mass()-wMass);
        if ( dm > dmAlt ) { w1 = w1Alt; w2 = w2Alt; std::swap(nu1, nu2); }
      }
      // Do the top combinations
      auto gentop1 = w1 + bjet1;
      auto gentop2 = w2 + bjet2;

      if ( true ) {
        const auto t1Alt = w1 + bjet2;
        const auto t2Alt = w2 + bjet1;

        const double tMass = 172.5;
        const double dm = std::abs(gentop1.mass()-tMass)+std::abs(gentop2.mass()-tMass);
        const double dmAlt = std::abs(t1Alt.mass()-tMass)+std::abs(t2Alt.mass()-tMass);
        if ( dm > dmAlt ) { gentop1 = t1Alt; gentop2 = t2Alt; std::swap(bjet1, bjet2); }
      }
      //if (gentop1.Pt() < gentop2.Pt()) { swap(gentop1, gentop2); }
      if (pseudoTopLeptonHandle->at(leptonIdxs[0]).charge() < 0){
    	swap(nu1, nu2);
    	swap(bjet1, bjet2);
    	swap(lepton1, lepton2);
    	swap(gentop1, gentop2);
      }	
      evInfo_.gen_pseudonu1 = ToTLorentzVector(nu1);
      evInfo_.gen_pseudonu2 = ToTLorentzVector(nu2);
      evInfo_.gen_pseudotop1 = ToTLorentzVector(gentop1);
      evInfo_.gen_pseudotop2 = ToTLorentzVector(gentop2);
      evInfo_.gen_pseudottbar = evInfo_.gen_pseudotop1 + evInfo_.gen_pseudotop2;
      evInfo_.gen_pseudottbar_dphi = evInfo_.gen_pseudotop1.DeltaPhi(evInfo_.gen_pseudotop2);

      evInfo_.gen_pseudolep1 = ToTLorentzVector(lepton1);
      evInfo_.gen_pseudolep2 = ToTLorentzVector(lepton2);
      evInfo_.gen_pseudodilep = evInfo_.gen_pseudolep1 + evInfo_.gen_pseudolep2;
      //evInfo_.gen_pseudolep1_pid = lepton1.pdgId();
      //evInfo_.gen_pseudolep2_pid = lepton2.pdgId();

      //if (bjet1.Pt() < bjet2.Pt()) { swap(bjet1, bjet2); }
      evInfo_.gen_pseudojet1 = ToTLorentzVector(bjet1);
      evInfo_.gen_pseudojet2 = ToTLorentzVector(bjet2);

    } while ( false );
  }
  return;
}
void TopEventCommon::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  TopEventInfo& evInfo_ = TopEventInfo::getInstance();
  evInfo_.run = iEvent.id().run();
  evInfo_.event = iEvent.id().event();

  const bool runOnMC = !iEvent.isRealData();

  for (int sys = 0; sys < nsys_e; ++sys){
    if (sys > 0 && !runOnMC) break;
    resetBr();
    if( sys == 0 ) genInfo(iEvent, iSetup);
    int terminate = runEventSelection(iEvent, iSetup, ttree_[sys] , sys);
    if ( terminate == -1 ) continue;
    else if ( terminate == -2 ) return;

    analyzeCustom(iEvent, iSetup, sys);
    ttree_[sys]->Fill();
  }
}
void TopEventCommon::analyzeCustom(const edm::Event& iEvent, const edm::EventSetup& iSetup, int sys) {

}

//define this as a plug-in
#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(TopEventCommon);
