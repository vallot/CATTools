#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/stream/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "fastjet/JetDefinition.hh"
#include "fastjet/ClusterSequence.hh"
#include "RecoJets/JetProducers/interface/JetSpecific.h"

#include "CATTools/CommonTools/interface/TTbarModeDefs.h"

using namespace std;

class PartonTopProducer : public edm::stream::EDProducer<>
{
public:
  PartonTopProducer(const edm::ParameterSet& pset);
  void produce(edm::Event& event, const edm::EventSetup& eventSetup) override;

private:
  const reco::Candidate* getLast(const reco::Candidate* p) const;
  reco::GenParticleRef buildGenParticle(const reco::Candidate* p, reco::GenParticleRefProd& refHandle,
                                        std::unique_ptr<reco::GenParticleCollection>& outColl) const;

  typedef reco::Particle::LorentzVector LorentzVector;

  const double jetMinPt_, jetMaxEta_;
  typedef fastjet::JetDefinition JetDef;
  std::shared_ptr<JetDef> fjDef_;
  const reco::Particle::Point genVertex_;

private:
  edm::EDGetTokenT<edm::View<reco::Candidate> > genParticleToken_;
};

PartonTopProducer::PartonTopProducer(const edm::ParameterSet& pset):
  jetMinPt_(pset.getParameter<double>("jetMinPt")),
  jetMaxEta_(pset.getParameter<double>("jetMaxEta")),
  genVertex_(0,0,0)
{
  genParticleToken_ = consumes<edm::View<reco::Candidate> >(pset.getParameter<edm::InputTag>("genParticles"));
  const double jetConeSize = pset.getParameter<double>("jetConeSize");
  fjDef_ = std::shared_ptr<JetDef>(new JetDef(fastjet::antikt_algorithm, jetConeSize));

  produces<reco::GenParticleCollection>();
  produces<int>("channel");
  produces<std::vector<int> >("modes");

  produces<reco::GenJetCollection>("qcdJets");
}

void PartonTopProducer::produce(edm::Event& event, const edm::EventSetup& eventSetup)
{
  std::unique_ptr<reco::GenParticleCollection> partons(new reco::GenParticleCollection);
  auto partonRefHandle = event.getRefBeforePut<reco::GenParticleCollection>();

  std::unique_ptr<int> channel(new int(cat::CH_NOTT));
  std::unique_ptr<std::vector<int> > modes(new std::vector<int>());

  std::unique_ptr<reco::GenJetCollection> qcdJets(new reco::GenJetCollection);

  edm::Handle<edm::View<reco::Candidate> > genParticleHandle;
  if ( event.isRealData() or !event.getByToken(genParticleToken_, genParticleHandle) ) {
    event.put(std::move(partons));
    event.put(std::move(channel), "channel");
    event.put(std::move(modes), "modes");
    event.put(std::move(qcdJets), "qcdJets");
    return;
  }
  
  // Collect top quarks and unstable B-hadrons
  std::vector<const reco::Candidate*> tQuarks;
  std::vector<int> qcdParticleIdxs;
  for ( size_t i=0, n=genParticleHandle->size(); i<n; ++i ) {
    const reco::Candidate& p = genParticleHandle->at(i);
    const int status = p.status();
    if ( status == 1 ) continue;

    // Collect parton level objects.
    const int absPdgId = abs(p.pdgId());
    if ( absPdgId == 6 ) {
      // top quark : select one 'after radiations'
      bool toKeep = true;
      if ( p.numberOfDaughters() == 0 ) toKeep = false;
      for ( size_t j=0, m=p.numberOfDaughters(); j<m; ++j ) {
        const int dauId = p.daughter(j)->pdgId();
        if ( dauId == p.pdgId() ) { toKeep = false; break; }
      }
      if ( toKeep ) tQuarks.push_back(&p);
    }
    else if ( absPdgId < 6 or absPdgId == 21 ) {
      // QCD particles : select one after parton shower, before hadronization
      bool toKeep = true;
      for ( size_t j=0, m=p.numberOfDaughters(); j<m; ++j ) {
        const int absDauId = abs(p.daughter(j)->pdgId());
        if ( absDauId < 6 or absPdgId == 21 ) { toKeep = false; break; }
      }
      if ( toKeep ) qcdParticleIdxs.push_back(i);
    }
  }
  // Build top quark decay tree in parton level
  // Also determine decay mode from parton level information
  size_t nElectron = 0, nMuon = 0, nTau = 0, nTauToLepton = 0;
  for ( int i=0, n=tQuarks.size(); i<n; ++i ) {
    const reco::Candidate* t = tQuarks.at(i);
    const reco::Candidate* tLast = getLast(t);
    reco::GenParticleRef tRef = buildGenParticle(tLast, partonRefHandle, partons);
  }

  for ( int i=0, n=tQuarks.size(); i<n; ++i ) {
    const reco::Candidate* tLast = getLast(tQuarks.at(i));
    reco::GenParticleRef tRef(partonRefHandle, i);

    const reco::Candidate* w = 0;
    const reco::Candidate* b = 0;
    for ( int j=0, m=tLast->numberOfDaughters(); j<m; ++j ) {
      const reco::Candidate* dau = tLast->daughter(j);
      const unsigned int dauAbsId = abs(dau->pdgId());
      if ( (dauAbsId == 24 or dauAbsId == 25 or dauAbsId == 23) and !w ) w = dau; // Include top-FCNC
      else if ( dauAbsId < 6 and !b ) b = dau;
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
    for ( int j=0, m=wLast->numberOfDaughters(); j<m; ++j ) {
      const reco::Candidate* dau = wLast->daughter(j);
      const unsigned int dauAbsId = abs(dau->pdgId());
      if ( abs(wLast->pdgId()) != 25 and dauAbsId > 16 ) continue; // Consider quarks and leptons only for W/Z decays
      // With the line above, we allow H->ff/GG/ZZ/WW.
      // wLast should be W/Z/H, nothing else.

      if ( !wDau1 ) wDau1 = dau;
      else if ( !wDau2 ) wDau2 = dau;
      else {
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
    std::vector<const reco::Candidate*> lepsFromTau;
    if ( abs(wDau1->pdgId()) == 15 ) {
      const reco::Candidate* tauLast = getLast(wDau1);
      for ( int j=0, m=tauLast->numberOfDaughters(); j<m; ++j ) {
        const reco::Candidate* dau = tauLast->daughter(j);
        const unsigned int dauAbsId = abs(dau->pdgId());
        if ( dauAbsId == 11 or dauAbsId == 13 ) lepsFromTau.push_back(dau);
      }
      // Cleanup the daughter lepton if net charge is zero. This happens if a conversion photon is radiated (tau->gamma+X, gamma->e+e-)
      // This happens in sub-per-mil level (observed 6 among 10000 events)
      const int sumQ = std::accumulate(lepsFromTau.begin(), lepsFromTau.end(), 0, [](int b, const reco::Candidate* a){return a->charge()+b;});
      if ( sumQ == 0 ) lepsFromTau.clear();
      // Sort daughter leptons, largest pT with consistent electric charge to the original tau at the front.
      std::sort(lepsFromTau.begin(), lepsFromTau.end(), [](const reco::Candidate* a, const reco::Candidate* b){return a->pt() > b->pt();});
      std::stable_sort(lepsFromTau.begin(), lepsFromTau.end(), [&](const reco::Candidate* a, const reco::Candidate* b){return a->charge() == tauLast->charge();});
      for ( auto lepFromTau : lepsFromTau ) {
        reco::GenParticleRef lepRef = buildGenParticle(lepFromTau, partonRefHandle, partons);
        partons->at(lepRef.key()).addMother(wDauRef1);
        partons->at(wDauRef1.key()).addDaughter(lepRef);
      }
    }
    int mode = cat::CH_HADRON;
    switch ( abs(wDau1->pdgId()) ) {
      case 11: ++nElectron; mode = cat::CH_ELECTRON; break;
      case 13: ++nMuon; mode = cat::CH_MUON; break;
      case 15:
        ++nTau; mode = cat::CH_TAU_HADRON;
        if ( !lepsFromTau.empty() ) {

          const reco::Candidate* lepFromTau = lepsFromTau.front();
          ++nTauToLepton;
          if ( abs(lepFromTau->pdgId()) == 13 ) {
            mode += 1;
            ++nMuon;
          }
          else {
            mode += 2;
            ++nElectron;
          }
        }
        break;
    }
    modes->push_back(mode);
  }

  if ( modes->size() == 2 ) {
    const int nLepton = nElectron + nMuon;
    if      ( nLepton == 0 ) *channel = cat::CH_FULLHADRON;
    else if ( nLepton == 1 ) *channel = cat::CH_SEMILEPTON;
    else if ( nLepton == 2 ) *channel = cat::CH_FULLLEPTON;
  }

  // Make genJets using particles after PS, but before hadronization
  std::vector<fastjet::PseudoJet> fjInputs;
  fjInputs.reserve(qcdParticleIdxs.size());
  for ( int i : qcdParticleIdxs ) {
    const auto& p = genParticleHandle->at(i);
    fjInputs.push_back(fastjet::PseudoJet(p.px(), p.py(), p.pz(), p.energy()));
    fjInputs.back().set_user_index(i);
  }
  fastjet::ClusterSequence fjClusterSeq(fjInputs, *fjDef_);
  std::vector<fastjet::PseudoJet> fjJets = fastjet::sorted_by_pt(fjClusterSeq.inclusive_jets(jetMinPt_));
  qcdJets->reserve(fjJets.size());
  for ( auto& fjJet : fjJets ) {
    if ( abs(fjJet.eta()) > jetMaxEta_ ) continue;
    const auto& fjCons = fjJet.constituents();
    std::vector<reco::CandidatePtr> cons;
    cons.reserve(fjCons.size());
    for ( auto con : fjCons ) {
      const size_t index = con.user_index();
      cons.push_back(genParticleHandle->ptrAt(index));
    }

    const LorentzVector jetP4(fjJet.px(), fjJet.py(), fjJet.pz(), fjJet.E());
    reco::GenJet qcdJet;
    reco::writeSpecific(qcdJet, jetP4, genVertex_, cons, eventSetup);

    qcdJets->push_back(qcdJet);
  }

  event.put(std::move(partons));
  event.put(std::move(channel), "channel");
  event.put(std::move(modes), "modes");
  event.put(std::move(qcdJets), "qcdJets");
}

const reco::Candidate* PartonTopProducer::getLast(const reco::Candidate* p) const
{
  int nDecay = 0;
  std::vector<const reco::Candidate*> sameCopies;
  for ( size_t i=0, n=p->numberOfDaughters(); i<n; ++i ) {
    const reco::Candidate* dau = p->daughter(i);
    const int dauId = dau->pdgId();
    if ( dauId == 22 or dauId == 21 ) continue;
    if ( p->pdgId() == dau->pdgId() ) sameCopies.push_back(dau);
    else ++nDecay;
  }
  if ( nDecay == 0 ) {
    for ( const auto dau : sameCopies ) {
      if ( p->pdgId() == dau->pdgId() ) return getLast(dau);
    }
  }
  return p;
}

reco::GenParticleRef PartonTopProducer::buildGenParticle(const reco::Candidate* p, reco::GenParticleRefProd& refHandle,
                                                         std::unique_ptr<reco::GenParticleCollection>& outColl) const
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
DEFINE_FWK_MODULE(PartonTopProducer);

