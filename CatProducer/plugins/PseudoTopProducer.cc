#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "DataFormats/Math/interface/LorentzVector.h"
#include "DataFormats/Candidate/interface/Candidate.h"
#include "DataFormats/Candidate/interface/CompositeCandidate.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/Math/interface/deltaR.h"

#include "CommonTools/Utils/interface/PtComparator.h"

#include "RecoJets/JetProducers/interface/JetSpecific.h"
#include "fastjet/JetDefinition.hh"
#include "fastjet/ClusterSequence.hh"
#include "fastjet/PseudoJet.hh"

#include "TMath.h"

#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "CommonTools/Utils/interface/TFileDirectory.h"
#include "TH1F.h"
#include "TH2F.h"

#include <memory>
#include <vector>
#include <set>

using namespace std;
using namespace edm;
using namespace reco;

class PseudoTopProducer : public edm::EDProducer
{
public:
  PseudoTopProducer(const edm::ParameterSet& pset);
  void produce(edm::Event& event, const edm::EventSetup& eventSetup) override;

private:
  bool isFromHadron(const reco::Candidate* p) const;
  bool isBHadron(const reco::Candidate* p) const;
  bool isBHadron(const unsigned int pdgId) const;
  bool isFromW(const reco::Candidate&) const;

  const reco::Candidate* getLast(const reco::Candidate* p);
  reco::GenParticleRef buildGenParticle(const reco::Candidate* p, reco::GenParticleRefProd& refHandle,
                                        std::auto_ptr<reco::GenParticleCollection>& outColl) const;

  typedef reco::Particle::LorentzVector LorentzVector;

private:
  edm::EDGetTokenT<edm::View<reco::Candidate> > finalStateToken_;
  edm::EDGetTokenT<edm::View<reco::Candidate> > genParticleToken_;
  const double leptonMinPt_, leptonMaxEta_, jetMinPt_, jetMaxEta_;
  const double wMass_, tMass_;

  typedef fastjet::JetDefinition JetDef;
  std::shared_ptr<JetDef> fjLepDef_, fjJetDef_;
  reco::Particle::Point genVertex_;

  bool doHist_;
  TH2F* hDecayMode_;
  TH1F* hSL_W1Mass_, * hSL_W2Mass_;
  TH1F* hSL_Top1Mass_, * hSL_Top2Mass_;
  TH1F* hDL_W1Mass_, * hDL_W2Mass_;
  TH1F* hDL_Top1Mass_, * hDL_Top2Mass_;
};

PseudoTopProducer::PseudoTopProducer(const edm::ParameterSet& pset):
  leptonMinPt_(pset.getParameter<double>("leptonMinPt")),
  leptonMaxEta_(pset.getParameter<double>("leptonMaxEta")),
  jetMinPt_(pset.getParameter<double>("jetMinPt")),
  jetMaxEta_(pset.getParameter<double>("jetMaxEta")),
  wMass_(pset.getParameter<double>("wMass")),
  tMass_(pset.getParameter<double>("tMass"))
{
  finalStateToken_ = consumes<edm::View<reco::Candidate> >(pset.getParameter<edm::InputTag>("finalStates"));
  genParticleToken_ = consumes<edm::View<reco::Candidate> >(pset.getParameter<edm::InputTag>("genParticles"));

  const double leptonConeSize = pset.getParameter<double>("leptonConeSize");
  const double jetConeSize = pset.getParameter<double>("jetConeSize");
  fjLepDef_ = std::shared_ptr<JetDef>(new JetDef(fastjet::antikt_algorithm, leptonConeSize));
  fjJetDef_ = std::shared_ptr<JetDef>(new JetDef(fastjet::antikt_algorithm, jetConeSize));

  genVertex_ = reco::Particle::Point(0,0,0);

  produces<reco::GenParticleCollection>("neutrinos");
  produces<reco::GenJetCollection>("leptons");
  produces<reco::GenJetCollection>("jets");

  produces<reco::GenParticleCollection>();
  produces<reco::GenParticleCollection>("partons");

  doHist_ = false;
  edm::Service<TFileService> fs;
  if ( fs.isAvailable() )
  {
    doHist_ = true;

    hSL_W1Mass_ = fs->make<TH1F>("hSL_W1Mass", "W mass;Mass (GeV/c^{2});Entries", 300, 0, 300);
    hSL_W2Mass_ = fs->make<TH1F>("hSL_W2Mass", "W mass;Mass (GeV/c^{2});Entries", 300, 0, 300);
    hSL_Top1Mass_ = fs->make<TH1F>("hSL_Top1Mass", "Top mass;Mass (GeV/c^{2});Entries", 300, 0, 300);
    hSL_Top2Mass_ = fs->make<TH1F>("hSL_Top2Mass", "Top mass;Mass (GeV/c^{2});Entries", 300, 0, 300);

    hDL_W1Mass_ = fs->make<TH1F>("hDL_W1Mass", "W mass;Mass (GeV/c^{2});Entries", 300, 0, 300);
    hDL_W2Mass_ = fs->make<TH1F>("hDL_W2Mass", "W mass;Mass (GeV/c^{2});Entries", 300, 0, 300);
    hDL_Top1Mass_ = fs->make<TH1F>("hDL_Top1Mass", "Top mass;Mass (GeV/c^{2});Entries", 300, 0, 300);
    hDL_Top2Mass_ = fs->make<TH1F>("hDL_Top2Mass", "Top mass;Mass (GeV/c^{2});Entries", 300, 0, 300);

    hDecayMode_ = fs->make<TH2F>("hDecayMode", "Decay mode;Parton nLepton;Particle nLepton", 4, 0, 4, 4, 0, 4);
  }
}

void PseudoTopProducer::produce(edm::Event& event, const edm::EventSetup& eventSetup)
{
  edm::Handle<edm::View<reco::Candidate> > finalStateHandle;
  event.getByToken(finalStateToken_, finalStateHandle);

  edm::Handle<edm::View<reco::Candidate> > genParticleHandle;
  event.getByToken(genParticleToken_, genParticleHandle);

  std::auto_ptr<reco::GenParticleCollection> neutrinos(new reco::GenParticleCollection);
  std::auto_ptr<reco::GenJetCollection> leptons(new reco::GenJetCollection);
  std::auto_ptr<reco::GenJetCollection> jets(new reco::GenJetCollection);
  auto neutrinosRefHandle = event.getRefBeforePut<reco::GenParticleCollection>("neutrinos");
  auto leptonsRefHandle = event.getRefBeforePut<reco::GenJetCollection>("leptons");
  auto jetsRefHandle = event.getRefBeforePut<reco::GenJetCollection>("jets");

  std::auto_ptr<reco::GenParticleCollection> pseudoTop(new reco::GenParticleCollection);
  auto pseudoTopRefHandle = event.getRefBeforePut<reco::GenParticleCollection>();

  std::auto_ptr<reco::GenParticleCollection> partons(new reco::GenParticleCollection);
  auto partonRefHandle = event.getRefBeforePut<reco::GenParticleCollection>("partons");

  // Collect top quarks and unstable B-hadrons
  std::vector<const reco::Candidate*> tQuarks;
  std::set<size_t> bHadronIdxs;
  for ( size_t i=0, n=genParticleHandle->size(); i<n; ++i )
  {
    const reco::Candidate& p = genParticleHandle->at(i);
    const int status = p.status();
    if ( status == 1 ) continue;

    // Collect B-hadrons, to be used in b tagging
    if ( isBHadron(&p) ) bHadronIdxs.insert(i);

    // Collect parton level objects. Ignore gluons and photons
    const int absPdgId = abs(p.pdgId());
    if ( absPdgId == 6 )
    {
      // top quark : select one 'after radiations'
      bool toKeep = true;
      if ( p.numberOfDaughters() == 0 ) toKeep = false;
      for ( size_t j=0, m=p.numberOfDaughters(); j<m; ++j )
      {
        const int dauId = p.daughter(j)->pdgId();
        if ( dauId == p.pdgId() ) { toKeep = false; break; }
      }
      if ( toKeep ) tQuarks.push_back(&p);
    }
  }
  // Build top quark decay tree in parton level
  // Also determine decay mode from parton level information
  size_t nElectron = 0, nMuon = 0, nTau = 0, nTauToLepton = 0;
  for ( int i=0, n=tQuarks.size(); i<n; ++i )
  {
    const reco::Candidate* t = tQuarks.at(i);
    const reco::Candidate* tLast = getLast(t);
    reco::GenParticleRef tRef = buildGenParticle(tLast, partonRefHandle, partons);

    const reco::Candidate* w = 0;
    const reco::Candidate* b = 0;
    for ( int j=0, m=tLast->numberOfDaughters(); j<m; ++j )
    {
      const reco::Candidate* dau = tLast->daughter(j);
      const unsigned int dauAbsId = abs(dau->pdgId());
      if ( dauAbsId == 24 and !w ) w = dau;
      else if ( dauAbsId == 5 and !b ) b = dau;
    }
    if ( !w or !b ) continue;
    reco::GenParticleRef wRef = buildGenParticle(w, partonRefHandle, partons);
    reco::GenParticleRef bRef = buildGenParticle(b, partonRefHandle, partons);
    partons->at(wRef.key()).addMother(tRef);
    partons->at(bRef.key()).addMother(tRef);
    partons->at(tRef.key()).addDaughter(wRef);
    partons->at(tRef.key()).addDaughter(bRef);

    // W decay products
    const reco::Candidate* wLast = getLast(w);
    const reco::Candidate* wDau1 = 0;
    const reco::Candidate* wDau2 = 0;
    for ( int j=0, m=wLast->numberOfDaughters(); j<m; ++j )
    {
      const reco::Candidate* dau = wLast->daughter(j);
      const unsigned int dauAbsId = abs(dau->pdgId());
      if ( dauAbsId > 16 ) continue; // Consider quarks and leptons only for W decays

      if ( !wDau1 ) wDau1 = dau;
      else if ( !wDau2 ) wDau2 = dau;
      else
      {
        cout << "--------------------------------------" << endl;
        cout << "WDECAY with more than 2 body!!! " << wLast->numberOfDaughters() << endl;
        cout << " dau1 = " << wDau1->pdgId() << " dau2 = " << wDau2->pdgId() << " skipped = " << dau->pdgId() << endl;
        cout << "--------------------------------------" << endl;
      }
    }
    if ( !wDau1 or !wDau2 ) continue;
    if ( abs(wDau1->pdgId()) > abs(wDau2->pdgId()) ) swap(wDau1, wDau2);
    reco::GenParticleRef wDauRef1 = buildGenParticle(wDau1, partonRefHandle, partons);
    reco::GenParticleRef wDauRef2 = buildGenParticle(wDau2, partonRefHandle, partons);
    partons->at(wDauRef1.key()).addMother(wRef);
    partons->at(wDauRef2.key()).addMother(wRef);
    partons->at(wRef.key()).addDaughter(wDauRef1);
    partons->at(wRef.key()).addDaughter(wDauRef2);

    // Special care for tau->lepton decays
    // Note : we do not keep neutrinos from tau decays (tau->W, nu_tau, W->l, nu_l)
    // Note : Up to 6 neutrinos from top decays if both W decays to taus and all taus go into leptonic decay chain
    const reco::Candidate* lepFromTau = 0;
    if ( abs(wDau1->pdgId()) == 15 )
    {
      const reco::Candidate* tauLast = getLast(wDau1);
      for ( int j=0, m=tauLast->numberOfDaughters(); j<m; ++j )
      {
        const reco::Candidate* dau = tauLast->daughter(j);
        const unsigned int dauAbsId = abs(dau->pdgId());
        if ( dauAbsId == 11 or dauAbsId == 13 )
        {
          if ( !lepFromTau ) lepFromTau = dau;
          else
          {
            cout << "--------------------------------------" << endl;
            cout << "TAU decay with more than 1 leptons!!!, nDau=" << tauLast->numberOfDaughters() << endl;
            cout << " dau = " << lepFromTau->pdgId() << " skipped = " << dau->pdgId() << endl;
            cout << "--------------------------------------" << endl;
          }
        }
      }
      if ( lepFromTau )
      {
        reco::GenParticleRef lepRef = buildGenParticle(lepFromTau, partonRefHandle, partons);
        partons->at(lepRef.key()).addMother(wDauRef1);
        partons->at(wDauRef1.key()).addDaughter(lepRef);
      }
    }
    switch ( abs(wDau1->pdgId()) )
    {
      case 11: ++nElectron; break;
      case 13: ++nMuon; break;
      case 15:
        ++nTau;
        if ( lepFromTau ) ++nTauToLepton;
        break;
    }
  }

  // Collect stable leptons and neutrinos
  size_t nStables = 0;
  std::vector<size_t> leptonIdxs;
  std::set<size_t> neutrinoIdxs;
  for ( size_t i=0, n=finalStateHandle->size(); i<n; ++i )
  {
    const reco::Candidate& p = finalStateHandle->at(i);
    const int absPdgId = abs(p.pdgId());
    if ( p.status() != 1 ) continue;

    ++nStables;
    if ( isFromHadron(&p) ) continue;
    if ( !isFromW(p) ) continue;
    switch ( absPdgId )
    {
      case 11: case 13: // Leptons
      case 22: // Photons
        leptonIdxs.push_back(i);
        break;
      case 12: case 14: case 16:
        neutrinoIdxs.insert(i);
        neutrinos->push_back(reco::GenParticle(p.charge(), p.p4(), p.vertex(), p.pdgId(), p.status(), true));
        break;
    }
  }

  // Sort neutrinos by pT.
  std::sort(neutrinos->begin(), neutrinos->end(), GreaterByPt<reco::Candidate>());

  // Make dressed leptons with anti-kt(0.1) algorithm
  //// Prepare input particle list
  std::vector<fastjet::PseudoJet> fjLepInputs;
  fjLepInputs.reserve(leptonIdxs.size());
  for ( auto index : leptonIdxs )
  {
    const reco::Candidate& p = finalStateHandle->at(index);
    if ( std::isnan(p.pt()) or p.pt() <= 0 ) continue;

    fjLepInputs.push_back(fastjet::PseudoJet(p.px(), p.py(), p.pz(), p.energy()));
    fjLepInputs.back().set_user_index(index);
  }

  //// Run the jet algorithm
  fastjet::ClusterSequence fjLepClusterSeq(fjLepInputs, *fjLepDef_);
  std::vector<fastjet::PseudoJet> fjLepJets = fastjet::sorted_by_pt(fjLepClusterSeq.inclusive_jets(leptonMinPt_));

  //// Build dressed lepton objects from the FJ output
  leptons->reserve(fjLepJets.size());
  std::set<size_t> lepDauIdxs; // keep lepton constituents to remove from GenJet construction
  for ( auto& fjJet : fjLepJets )
  {
    if ( abs(fjJet.eta()) > leptonMaxEta_ ) continue;

    // Get jet constituents from fastJet
    const std::vector<fastjet::PseudoJet> fjConstituents = fastjet::sorted_by_pt(fjJet.constituents());
    // Convert to CandidatePtr
    std::vector<reco::CandidatePtr> constituents;
    reco::CandidatePtr lepCand;
    for ( auto& fjConstituent : fjConstituents )
    {
      const size_t index = fjConstituent.user_index();
      reco::CandidatePtr cand = finalStateHandle->ptrAt(index);
      const int absPdgId = abs(cand->pdgId());
      if ( absPdgId == 11 or absPdgId == 13 )
      {
        if ( lepCand.isNonnull() and lepCand->pt() > cand->pt() ) continue; // Choose one with highest pt
        lepCand = cand;
      }
      constituents.push_back(cand);
    }
    if ( lepCand.isNull() ) continue;
    if ( lepCand->pt() < fjJet.pt()/2 ) continue; // Central lepton must be the major component

    const LorentzVector jetP4(fjJet.px(), fjJet.py(), fjJet.pz(), fjJet.E());
    reco::GenJet lepJet;
    reco::writeSpecific(lepJet, jetP4, genVertex_, constituents, eventSetup);

    lepJet.setPdgId(lepCand->pdgId());
    lepJet.setCharge(lepCand->charge());

    const double jetArea = fjJet.has_area() ? fjJet.area() : 0;
    lepJet.setJetArea(jetArea);

    leptons->push_back(lepJet);

    // Keep constituent indices to be used in the next step.
    for ( auto& fjConstituent : fjConstituents )
    {
      lepDauIdxs.insert(fjConstituent.user_index());
    }
  }

  // Now proceed to jets.
  // Jets: anti-kt excluding the e, mu, nu, and photons in selected leptons.
  //// Prepare input particle list. Remove particles used in lepton clusters, neutrinos
  std::vector<fastjet::PseudoJet> fjJetInputs;
  fjJetInputs.reserve(nStables);
  for ( size_t i=0, n=finalStateHandle->size(); i<n; ++i )
  {
    const reco::Candidate& p = finalStateHandle->at(i);
    if ( p.status() != 1 ) continue;
    if ( std::isnan(p.pt()) or p.pt() <= 0 ) continue;

    if ( neutrinoIdxs.find(i) != neutrinoIdxs.end() ) continue;
    if ( lepDauIdxs.find(i) != lepDauIdxs.end() ) continue;

    fjJetInputs.push_back(fastjet::PseudoJet(p.px(), p.py(), p.pz(), p.energy()));
    fjJetInputs.back().set_user_index(i);
  }
  //// Also don't forget to put B hadrons
  for ( auto index : bHadronIdxs )
  {
    const reco::Candidate& p = genParticleHandle->at(index);
    if ( std::isnan(p.pt()) or p.pt() <= 0 ) continue;

    const double scale = 1e-20/p.p();
    fjJetInputs.push_back(fastjet::PseudoJet(p.px()*scale, p.py()*scale, p.pz()*scale, p.energy()*scale));
    fjJetInputs.back().set_user_index(index);
  }

  //// Run the jet algorithm
  fastjet::ClusterSequence fjJetClusterSeq(fjJetInputs, *fjJetDef_);
  std::vector<fastjet::PseudoJet> fjJets = fastjet::sorted_by_pt(fjJetClusterSeq.inclusive_jets(jetMinPt_));

  /// Build jets
  jets->reserve(fjJets.size());
  std::vector<size_t> bjetIdxs, ljetIdxs;
  for ( auto& fjJet : fjJets )
  {
    if ( abs(fjJet.eta()) > jetMaxEta_ ) continue;

    // Get jet constituents from fastJet
    const std::vector<fastjet::PseudoJet> fjConstituents = fastjet::sorted_by_pt(fjJet.constituents());
    // Convert to CandidatePtr
    std::vector<reco::CandidatePtr> constituents;
    bool hasBHadron = false;
    for ( size_t j=0, m=fjConstituents.size(); j<m; ++j )
    { 
      const size_t index = fjConstituents[j].user_index();
      if ( bHadronIdxs.find(index) != bHadronIdxs.end() ) hasBHadron = true;
      reco::CandidatePtr cand = finalStateHandle->ptrAt(index);
      constituents.push_back(cand);
    }
    
    const LorentzVector jetP4(fjJet.px(), fjJet.py(), fjJet.pz(), fjJet.E());
    reco::GenJet genJet;
    reco::writeSpecific(genJet, jetP4, genVertex_, constituents, eventSetup);

    const double jetArea = fjJet.has_area() ? fjJet.area() : 0;
    genJet.setJetArea(jetArea);
    if ( hasBHadron )
    {
      genJet.setPdgId(5);
      bjetIdxs.push_back(jets->size());
    }
    else
    {
      ljetIdxs.push_back(jets->size());
    }

    jets->push_back(genJet);
  }

  // Every building blocks are ready. Continue to pseudo-W and pseudo-top combination
  // NOTE : A C++ trick, use do-while instead of long-nested if-statements.
  do
  {
    // should be generated with nTop=2. This is baseline.
    if ( tQuarks.size() != 2 ) break;

    // Note : we will do dilepton or semilepton channel only, no tau.
    //        channels are identified from parton level information
    //        I believe this is not MC dependent too much.
    if ( nTau > 0 ) break;
    const size_t nLepton = nElectron+nMuon;

    // Collect leptonic-decaying W's
    std::vector<std::pair<int, int> > wCandIdxs;
    for ( auto lep = leptons->begin(); lep != leptons->end(); ++lep )
    {
      const size_t iLep = lep-leptons->begin();
      for ( auto nu = neutrinos->begin(); nu != neutrinos->end(); ++nu )
      {
        if ( abs(lep->pdgId()+nu->pdgId()) != 1 ) continue; // Enforce to conserve flavour, this reduces combinatorial bkg
        const double m = (lep->p4()+nu->p4()).mass();
        if ( m > 300 ) continue; // Raw mass cut, reduce combinatorial background

        const size_t iNu = nu-neutrinos->begin();
        wCandIdxs.push_back(make_pair(iLep, iNu));
      }
    }
    if ( doHist_ ) hDecayMode_->Fill(nLepton, wCandIdxs.size());

    if      ( nLepton == 0 or wCandIdxs.empty() ) break; // Skip full hadronic channel
    else if ( nLepton == 1 and wCandIdxs.size() >= 1 ) // Semi-leptonic channel
    {
      double dm = 1e9; // Note: this will be re-used later.
      int bestLepIdx = -1, bestNuIdx = -1;
      for ( auto wCandIdx : wCandIdxs )
      {
        const int lepIdx = wCandIdx.first;
        const int nuIdx  = wCandIdx.second;
        const LorentzVector& lepLVec = leptons->at(lepIdx).p4();
        const LorentzVector& nuLVec = neutrinos->at(nuIdx).p4();
        const double mW = (lepLVec + nuLVec).mass();
        if ( mW > 300 ) continue;

        const double dmNew = abs(mW-wMass_);
        if ( dmNew < dm )
        {
          dm = dmNew;
          bestLepIdx = lepIdx;
          bestNuIdx = nuIdx;
        }
      }
      if ( bestLepIdx < 0 ) break; // Actually this should never happen
      const auto& lepton = leptons->at(bestLepIdx);
      const auto& neutrino = neutrinos->at(bestNuIdx);
      const LorentzVector w1LVec = lepton.p4()+neutrino.p4();
      if ( doHist_ ) hSL_W1Mass_->Fill(w1LVec.mass());

      // Continue to hadronic W
      // We also collect b jets to be used in the next step.
      dm = 1e9; // Reset dM to very large value.
      int bestJ1Idx = -1, bestJ2Idx = -1;
      for ( auto j1Idx = ljetIdxs.begin(); j1Idx != ljetIdxs.end(); ++j1Idx )
      {
        const auto& j1 = jets->at(*j1Idx);
        for ( auto j2Idx = j1Idx+1; j2Idx != ljetIdxs.end(); ++j2Idx )
        {
          const auto& j2 = jets->at(*j2Idx);
          const double mW = (j1.p4()+j2.p4()).mass();
          if ( mW > 300 ) continue;

          const double dmNew = abs(mW-wMass_);
          if ( dmNew < dm )
          {
            dm = dmNew;
            bestJ1Idx = *j1Idx;
            bestJ2Idx = *j2Idx;
          }
        }
      }
      if ( bestJ1Idx < 0 ) break;
      const auto& wJet1 = jets->at(bestJ1Idx);
      const auto& wJet2 = jets->at(bestJ2Idx);
      const LorentzVector w2LVec = wJet1.p4() + wJet2.p4();
      if ( doHist_ ) hSL_W2Mass_->Fill(w2LVec.mass());

      // Now we have leptonic W and hadronic W.
      // Contiue to top quarks
      dm = 1e9; // Reset once again for top combination.
      int bestB1Idx = -1, bestB2Idx = -1;
      if ( bjetIdxs.size() < 2 ) break;
      for ( auto b1Idx : bjetIdxs )
      {
        const double t1Mass = (w1LVec + jets->at(b1Idx).p4()).mass();
        if ( t1Mass > 300 ) continue;
        for ( auto b2Idx : bjetIdxs )
        {
          if ( b1Idx == b2Idx ) continue;
          const double t2Mass = (w2LVec + jets->at(b2Idx).p4()).mass();
          if ( t2Mass > 300 ) continue;

          const double dmNew = abs(t1Mass-tMass_) + abs(t2Mass-tMass_);
          if ( dmNew < dm )
          {
            dm = dmNew;
            bestB1Idx = b1Idx;
            bestB2Idx = b2Idx;
          }
        }
      }
      if ( bestB1Idx < 0 ) break;
      const auto& bJet1 = jets->at(bestB1Idx);
      const auto& bJet2 = jets->at(bestB2Idx);
      const LorentzVector t1LVec = w1LVec + bJet1.p4();
      const LorentzVector t2LVec = w2LVec + bJet2.p4();
      if ( doHist_ )
      {
        hSL_Top1Mass_->Fill(t1LVec.mass());
        hSL_Top2Mass_->Fill(t2LVec.mass());
      }

      // Put all of them into candidate collection
      if ( true ) // Trick to restrict variables' scope to avoid collision
      {
        const int lepQ = lepton.charge();

        reco::GenParticle t1(lepQ*2/3, t1LVec, genVertex_, lepQ*6, 3, false);
        reco::GenParticle w1(lepQ, w1LVec, genVertex_, lepQ*24, 3, true);
        reco::GenParticle b1(0, bJet1.p4(), genVertex_, lepQ*5, 1, true);
        reco::GenParticle l1(lepQ, lepton.p4(), genVertex_, lepton.pdgId(), 1, true);
        reco::GenParticle n1(0, neutrino.p4(), genVertex_, neutrino.pdgId(), 1, true);

        reco::GenParticle t2(-lepQ*2/3, t2LVec, genVertex_, -lepQ*6, 3, false);
        reco::GenParticle w2(-lepQ, w2LVec, genVertex_, -lepQ*24, 3, true);
        reco::GenParticle b2(0, bJet2.p4(), genVertex_, -lepQ*5, 1, true);
        reco::GenParticle qA(0, wJet1.p4(), genVertex_, 1, 1, true);
        reco::GenParticle qB(0, wJet2.p4(), genVertex_, 2, 1, true);

        pseudoTop->push_back(t1);
        pseudoTop->push_back(w1);
        pseudoTop->push_back(b1);
        pseudoTop->push_back(l1);
        pseudoTop->push_back(n1);

        pseudoTop->push_back(t2);
        pseudoTop->push_back(w2);
        pseudoTop->push_back(b2);
        pseudoTop->push_back(qA);
        pseudoTop->push_back(qB);
      }
    }
    else if ( nLepton == 2 and wCandIdxs.size() >= 2 ) // Dilepton channel
    {
      double dm = 1e9;
      int bestLep1Idx = -1, bestLep2Idx = -1, bestNu1Idx = -1, bestNu2Idx = -1;
      for ( auto i : wCandIdxs )
      {
        const auto& lepton1 = leptons->at(i.first);
        if ( lepton1.charge() < 0 ) continue;
        const auto& neutrino1 = neutrinos->at(i.second);
        const double mW1 = (lepton1.p4()+neutrino1.p4()).mass();
        if ( mW1 > 300 ) continue;
        for ( auto j : wCandIdxs )
        {
          if ( i == j ) continue;
          const auto& lepton2 = leptons->at(j.first);
          if ( lepton2.charge() > 0 ) continue;
          const auto& neutrino2 = neutrinos->at(j.second);
          const double mW2 = (lepton2.p4()+neutrino2.p4()).mass();
          if ( mW2 > 300 ) continue;

          const double dmNew = abs(mW1-wMass_) + abs(mW2-wMass_);
          if ( dmNew < dm )
          {
            dm = dmNew;
            bestLep1Idx = i.first;
            bestLep2Idx = j.first;
            bestNu1Idx = i.second;
            bestNu2Idx = j.second;
          }
        }
      }
      if ( bestLep1Idx < 0 ) break;
      const auto& lepton1 = leptons->at(bestLep1Idx);
      const auto& lepton2 = leptons->at(bestLep2Idx);
      const auto& neutrino1 = neutrinos->at(bestNu1Idx);
      const auto& neutrino2 = neutrinos->at(bestNu2Idx);
      const LorentzVector w1LVec = lepton1.p4()+neutrino1.p4();
      const LorentzVector w2LVec = lepton2.p4()+neutrino2.p4();
      if ( doHist_ )
      {
        hDL_W1Mass_->Fill(w1LVec.mass());
        hDL_W2Mass_->Fill(w2LVec.mass());
      }

      // Contiue to top quarks
      dm = 1e9; // Reset once again for top combination.
      int bestB1Idx = -1, bestB2Idx = -1;
      if ( bjetIdxs.size() < 2 ) break;
      for ( auto b1Idx : bjetIdxs )
      {
        const double t1Mass = (w1LVec + jets->at(b1Idx).p4()).mass();
        if ( t1Mass > 300 ) continue;
        for ( auto b2Idx : bjetIdxs )
        {
          if ( b1Idx == b2Idx ) continue;
          const double t2Mass = (w2LVec + jets->at(b2Idx).p4()).mass();
          if ( t2Mass > 300 ) continue;

          const double dmNew = abs(t1Mass-tMass_) + abs(t2Mass-tMass_);
          if ( dmNew < dm )
          {
            dm = dmNew;
            bestB1Idx = b1Idx;
            bestB2Idx = b2Idx;
          }
        }
      }
      if ( bestB1Idx < 0 ) break;
      const auto& bJet1 = jets->at(bestB1Idx);
      const auto& bJet2 = jets->at(bestB2Idx);
      const LorentzVector t1LVec = w1LVec + bJet1.p4();
      const LorentzVector t2LVec = w2LVec + bJet2.p4();
      if ( doHist_ )
      {
        hDL_Top1Mass_->Fill(t1LVec.mass());
        hDL_Top2Mass_->Fill(t2LVec.mass());
      }

      // Put all of them into candidate collection
      if ( true ) // Trick to restrict variables' scope to avoid collision
      {
        const int lep1Q = lepton1.charge();
        const int lep2Q = lepton2.charge();

        reco::GenParticle t1(lep1Q*2/3, t1LVec, genVertex_, lep1Q*6, 3, true);
        reco::GenParticle w1(lep1Q, w1LVec, genVertex_, lep1Q*24, 3, true);
        reco::GenParticle b1(0, bJet1.p4(), genVertex_, lep1Q*5, 1, true);
        reco::GenParticle l1(lep1Q, lepton1.p4(), genVertex_, lepton1.pdgId(), 1, true);
        reco::GenParticle n1(0, neutrino1.p4(), genVertex_, neutrino1.pdgId(), 1, true);

        reco::GenParticle t2(lep2Q*2/3, t2LVec, genVertex_, lep2Q*6, 3, true);
        reco::GenParticle w2(lep2Q, w2LVec, genVertex_, lep2Q*24, 3, true);
        reco::GenParticle b2(0, bJet2.p4(), genVertex_, lep2Q*5, 1, true);
        reco::GenParticle l2(0, lepton2.p4(), genVertex_, lepton2.pdgId(), 1, true);
        reco::GenParticle n2(0, neutrino2.p4(), genVertex_, neutrino2.pdgId(), 1, true);

        pseudoTop->push_back(t1);
        pseudoTop->push_back(w1);
        pseudoTop->push_back(b1);
        pseudoTop->push_back(l1);
        pseudoTop->push_back(n1);

        pseudoTop->push_back(t2);
        pseudoTop->push_back(w2);
        pseudoTop->push_back(b2);
        pseudoTop->push_back(l2);
        pseudoTop->push_back(n2);
      }
    }

    if ( pseudoTop->size() == 10 ) // If pseudtop decay tree is completed
    {
      // t->W+b
      pseudoTop->at(0).addDaughter(reco::GenParticleRef(pseudoTopRefHandle, 1)); // t->W
      pseudoTop->at(0).addDaughter(reco::GenParticleRef(pseudoTopRefHandle, 2)); // t->b
      pseudoTop->at(1).addMother(reco::GenParticleRef(pseudoTopRefHandle, 0)); // t->W
      pseudoTop->at(2).addMother(reco::GenParticleRef(pseudoTopRefHandle, 0)); // t->b

      // W->lv or W->jj
      pseudoTop->at(1).addDaughter(reco::GenParticleRef(pseudoTopRefHandle, 3));
      pseudoTop->at(1).addDaughter(reco::GenParticleRef(pseudoTopRefHandle, 4));
      pseudoTop->at(3).addMother(reco::GenParticleRef(pseudoTopRefHandle, 1));
      pseudoTop->at(4).addMother(reco::GenParticleRef(pseudoTopRefHandle, 1));

      // tbar->W-b
      pseudoTop->at(5).addDaughter(reco::GenParticleRef(pseudoTopRefHandle, 6));
      pseudoTop->at(5).addDaughter(reco::GenParticleRef(pseudoTopRefHandle, 7));
      pseudoTop->at(6).addMother(reco::GenParticleRef(pseudoTopRefHandle, 5));
      pseudoTop->at(7).addMother(reco::GenParticleRef(pseudoTopRefHandle, 5));

      // W->jj
      pseudoTop->at(6).addDaughter(reco::GenParticleRef(pseudoTopRefHandle, 8));
      pseudoTop->at(6).addDaughter(reco::GenParticleRef(pseudoTopRefHandle, 9));
      pseudoTop->at(8).addMother(reco::GenParticleRef(pseudoTopRefHandle, 6));
      pseudoTop->at(9).addMother(reco::GenParticleRef(pseudoTopRefHandle, 6));
    }
  } while (false);

  event.put(neutrinos, "neutrinos");
  event.put(leptons, "leptons");
  event.put(jets, "jets");

  event.put(pseudoTop);
  event.put(partons, "partons");
}

bool PseudoTopProducer::isFromW(const reco::Candidate& p) const
{
  const reco::Candidate* mother = &p;
  while ( (mother = mother->mother()) != 0 )
  {
    if ( abs(mother->pdgId()) != 24 ) continue;
    if ( mother->status() == 3 ) return true; // Pythia6 hard interaction
    if ( mother->status()>20 && mother->status()<30 ) return true; // Pythia8 hard interaction
  }
  return false;
}

const reco::Candidate* PseudoTopProducer::getLast(const reco::Candidate* p)
{
  for ( size_t i=0, n=p->numberOfDaughters(); i<n; ++i )
  {
    const reco::Candidate* dau = p->daughter(i);
    if ( p->pdgId() == dau->pdgId() ) return getLast(dau);
  }
  return p;
}

bool PseudoTopProducer::isFromHadron(const reco::Candidate* p) const
{
  for ( size_t i=0, n=p->numberOfMothers(); i<n; ++i )
  {
    const reco::Candidate* mother = p->mother(i);
    if ( mother->numberOfMothers() == 0 ) continue; // Skip incident beam
    const int pdgId = abs(mother->pdgId());

    if ( pdgId > 100 ) return true;
    else if ( isFromHadron(mother) ) return true;
  }
  return false;
}

bool PseudoTopProducer::isBHadron(const reco::Candidate* p) const
{
  const unsigned int absPdgId = abs(p->pdgId());
  if ( !isBHadron(absPdgId) ) return false;

  // Do not consider this particle if it has B hadron daughter
  // For example, B* -> B0 + photon; then we drop B* and take B0 only
  for ( int i=0, n=p->numberOfDaughters(); i<n; ++i )
  {
    const reco::Candidate* dau = p->daughter(i);
    if ( isBHadron(abs(dau->pdgId())) ) return false;
  }

  return true;
}

bool PseudoTopProducer::isBHadron(const unsigned int absPdgId) const
{
  if ( absPdgId <= 100 ) return false; // Fundamental particles and MC internals
  if ( absPdgId >= 1000000000 ) return false; // Nuclei, +-10LZZZAAAI

  // General form of PDG ID is 7 digit form
  // +- n nr nL nq1 nq2 nq3 nJ
  //const int nJ = absPdgId % 10; // Spin
  const int nq3 = (absPdgId / 10) % 10;
  const int nq2 = (absPdgId / 100) % 10;
  const int nq1 = (absPdgId / 1000) % 10;

  if ( nq3 == 0 ) return false; // Diquarks
  if ( nq1 == 0 and nq2 == 5 ) return true; // B mesons
  if ( nq1 == 5 ) return true; // B baryons

  return false;
}

reco::GenParticleRef PseudoTopProducer::buildGenParticle(const reco::Candidate* p, reco::GenParticleRefProd& refHandle,
                                                               std::auto_ptr<reco::GenParticleCollection>& outColl) const
{
  reco::GenParticle pOut(*dynamic_cast<const reco::GenParticle*>(p));
  pOut.clearMothers();
  pOut.clearDaughters();
  pOut.resetMothers(refHandle.id());
  pOut.resetDaughters(refHandle.id());

  outColl->push_back(pOut);

  return reco::GenParticleRef(refHandle, outColl->size()-1);
}

#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(PseudoTopProducer);

