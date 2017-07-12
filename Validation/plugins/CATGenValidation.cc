#include "FWCore/Framework/interface/one/EDAnalyzer.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/Math/interface/deltaR.h"
#include "DataFormats/Math/interface/deltaPhi.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "CATTools/DataFormats/interface/GenWeights.h"

#include "TH1F.h"
#include "TH2F.h"

#include <iostream>

using namespace std;

class CATGenValidation : public edm::one::EDAnalyzer<edm::one::SharedResources>
{
public:
  CATGenValidation(const edm::ParameterSet& pset);
  void analyze(const edm::Event& event, const edm::EventSetup& eventSetup) override;

private:
  bool isFirst(const reco::Candidate& p) const {
    const int nMo = p.numberOfMothers();
    if ( nMo == 0 ) return true;

    for ( int i=0; i<nMo; ++i ) {
      const reco::Candidate* mo = p.mother(i);
      if ( !mo ) continue;
      if ( mo->pdgId() == p.pdgId() ) return false;
    }
    return true;
  };
  bool isLast(const reco::Candidate& p) const {
    const int nDa = p.numberOfDaughters();
    if ( nDa == 0 ) return true;

    for ( int i=0; i<nDa; ++i ) {
      const reco::Candidate* da = p.daughter(i);
      if ( !da ) continue;
      if ( da->pdgId() == p.pdgId() ) return false;
    }
    return true;
  };

private:
  typedef std::vector<float> vfloat;

  edm::EDGetTokenT<cat::GenWeights> genWeightsToken_;
  edm::EDGetTokenT<reco::GenParticleCollection> genParticlesToken_;
  edm::EDGetTokenT<vfloat> scaleupWeightsToken_, scaledownWeightsToken_;
  edm::EDGetTokenT<vfloat> pdfWeightsToken_, otherWeightsToken_;

  TH1D* hQscale_;
  TH2D* hPartonId_;

  TH1D* hWeight_, * hWeight_LHE_, * hWeight_Norm_;
  TH1D* hWeights_scaleup_, * hWeights_scaledown_, * hWeights_pdf_, * hWeights_others_;

  // GenParticle multiplicity
  TH1D* hN_All_; // Everything in the event
  TH1D* hN_Final_; // Final state particles
  TH1D* hN_FirstCopy_, * hN_LastCopy_; // First/last copy of particles
  TH1D* hN_UnstableFinal_; // Final state particles with wrong status code

  TH1D* hN_FirstT_, * hN_FirstW_, * hN_FirstZ_, * hN_FirstH_;
  TH1D* hN_LastT_, * hN_LastW_, * hN_LastZ_, * hN_LastH_;
  TH1D* hN_FirstB_, * hN_FirstG_, * hN_FirstL_;
  TH1D* hN_LastB_, * hN_LastG_, * hN_LastL_;

  TH1D* hM_FirstT_, * hM_FirstW_, * hM_FirstZ_, * hM_FirstH_;
  TH1D* hM_LastT_, * hM_LastW_, * hM_LastZ_, * hM_LastH_;

};

CATGenValidation::CATGenValidation(const edm::ParameterSet& pset)
{
  genWeightsToken_ = consumes<cat::GenWeights>(pset.getParameter<edm::InputTag>("weight"));
  genParticlesToken_ = consumes<reco::GenParticleCollection>(pset.getParameter<edm::InputTag>("genParticles"));

  scaleupWeightsToken_ = consumes<vfloat>(pset.getParameter<edm::InputTag>("scaleupWeights"));
  scaledownWeightsToken_ = consumes<vfloat>(pset.getParameter<edm::InputTag>("scaledownWeights"));
  pdfWeightsToken_ = consumes<vfloat>(pset.getParameter<edm::InputTag>("pdfWeights"));
  otherWeightsToken_ = consumes<vfloat>(pset.getParameter<edm::InputTag>("otherWeights"));

  usesResource("TFileService");
  edm::Service<TFileService> fs;

  hQscale_ = fs->make<TH1D>("hQscale", "Qscale;Q scale;Events", 100, 0, 1000);
  hPartonId_ = fs->make<TH2D>("hPartonId", "PartonIds;parton1;parton2;Events", 15, -7.5, 7.5, 15, -7.5, 7.5);

  hWeight_ = fs->make<TH1D>("hWeight", "GenWeight;GenWeight;Events", 100, -10, 10);
  hWeight_LHE_ = fs->make<TH1D>("hWeight_LHE", "From LHE;GenWeight;Events", 100, -10, 10);
  hWeight_Norm_ = fs->make<TH1D>("hWeight_Norm", "Normalized by 1;GenWeight;Events", 100, -10, 10);

  hWeights_scaleup_ = fs->make<TH1D>("hWeights_scaleup", "Scaleup;GenWeight;Events", 100, -10, 10);
  hWeights_scaledown_ = fs->make<TH1D>("hWeights_scaledown", "Scaledown;GenWeight;Events", 100, -10, 10);
  hWeights_pdf_ = fs->make<TH1D>("hWeights_pdf", "PDF;GenWeight;Events", 100, -10, 10);
  hWeights_others_ = fs->make<TH1D>("hWeights_other", "Others;GenWeight;Events", 100, -10, 10);

  hN_All_ = fs->make<TH1D>("hN_All", "All particles;GenParticle multiplicity;Events", 500, 0, 500);
  hN_Final_ = fs->make<TH1D>("hN_Final", "Final state particles;GenParticle multiplicity;Events", 500, 0, 500);
  hN_UnstableFinal_ = fs->make<TH1D>("hN_UnstableFinal", "Unstable final state particles;GenParticle multiplicity;Events", 500, 0, 500);
  hN_FirstCopy_ = fs->make<TH1D>("hN_FirstCopy", "First copy of particles;GenParticle multiplicity;Events", 500, 0, 500);
  hN_LastCopy_ = fs->make<TH1D>("hN_LastCopy", "Last copy of particles;GenParticle multiplicity;Events", 500, 0, 500);

  hN_FirstT_ = fs->make<TH1D>("hN_FirstT", "Top quarks;Particle multiplicity;Events", 10, 0, 10);
  hN_FirstW_ = fs->make<TH1D>("hN_FirstW", "W bosons;Particle multiplicity;Events", 10, 0, 10);
  hN_FirstZ_ = fs->make<TH1D>("hN_FirstZ", "Z bosons;Particle multiplicity;Events", 10, 0, 10);
  hN_FirstH_ = fs->make<TH1D>("hN_FirstH", "Higgs;Particle multiplicity;Events", 10, 0, 10);
  hN_FirstB_ = fs->make<TH1D>("hN_FirstB", "B quarks;Particle multiplicity;Events", 10, 0, 10);
  hN_FirstG_ = fs->make<TH1D>("hN_FirstG", "prompt photons;Particle multiplicity;Events", 10, 0, 10);
  hN_FirstL_ = fs->make<TH1D>("hN_FirstL", "prompt leptons;Particle multiplicity;Events", 10, 0, 10);

  hM_FirstT_ = fs->make<TH1D>("hM_FirstT", "Top quarks;Mass (GeV);Entries per 1GeV", 200, 0, 200);
  hM_FirstW_ = fs->make<TH1D>("hM_FirstW", "W bosons;Mass (GeV);Entries per 1GeV", 200, 0, 200);
  hM_FirstZ_ = fs->make<TH1D>("hM_FirstZ", "Z bosons;Mass (GeV);Entries per 1GeV", 200, 0, 200);
  hM_FirstH_ = fs->make<TH1D>("hM_FirstH", "Higgs;Mass (GeV);Entries per 1GeV", 200, 0, 200);

  hN_LastT_ = fs->make<TH1D>("hN_LastT", "Top quarks;Particle multiplicity;Events", 10, 0, 10);
  hN_LastW_ = fs->make<TH1D>("hN_LastW", "W bosons;Particle multiplicity;Events", 10, 0, 10);
  hN_LastZ_ = fs->make<TH1D>("hN_LastZ", "Z bosons;Particle multiplicity;Events", 10, 0, 10);
  hN_LastH_ = fs->make<TH1D>("hN_LastH", "Higgs;Particle multiplicity;Events", 10, 0, 10);
  hN_LastB_ = fs->make<TH1D>("hN_LastB", "B quarks;Particle multiplicity;Events", 10, 0, 10);
  hN_LastG_ = fs->make<TH1D>("hN_LastG", "prompt photons;Particle multiplicity;Events", 10, 0, 10);
  hN_LastL_ = fs->make<TH1D>("hN_LastL", "prompt leptons;Particle multiplicity;Events", 10, 0, 10);

  hM_LastT_ = fs->make<TH1D>("hM_LastT", "Top quarks;Mass (GeV);Entries per 1GeV", 200, 0, 200);
  hM_LastW_ = fs->make<TH1D>("hM_LastW", "W bosons;Mass (GeV);Entries per 1GeV", 200, 0, 200);
  hM_LastZ_ = fs->make<TH1D>("hM_LastZ", "Z bosons;Mass (GeV);Entries per 1GeV", 200, 0, 200);
  hM_LastH_ = fs->make<TH1D>("hM_LastH", "Higgs;Mass (GeV);Entries per 1GeV", 200, 0, 200);

}

void CATGenValidation::analyze(const edm::Event& event, const edm::EventSetup& eventSetup)
{
  edm::Handle<cat::GenWeights> genWeightsHandle;
  event.getByToken(genWeightsToken_, genWeightsHandle);
  const float weight = genWeightsHandle->genWeight();

  edm::Handle<vfloat> scaleupWeightsHandle, scaledownWeightsHandle;
  edm::Handle<vfloat> pdfWeightsHandle, otherWeightsHandle;
  event.getByToken(scaleupWeightsToken_, scaleupWeightsHandle);
  event.getByToken(scaledownWeightsToken_, scaledownWeightsHandle);
  event.getByToken(pdfWeightsToken_, pdfWeightsHandle);
  event.getByToken(otherWeightsToken_, otherWeightsHandle);

  edm::Handle<reco::GenParticleCollection> genParticlesHandle;
  event.getByToken(genParticlesToken_, genParticlesHandle);

  hQscale_->Fill(genWeightsHandle->qScale(), weight);
  auto gluonToZero = [](const int id){ return id == 21 ? 0 : id; };
  hPartonId_->Fill(gluonToZero(genWeightsHandle->id1()), gluonToZero(genWeightsHandle->id2()), weight);

  hWeight_->Fill(weight);
  hWeight_LHE_->Fill(genWeightsHandle->lheWeight());
  hWeight_Norm_->Fill(weight == 0 ? 0 : weight/std::abs(weight));

  if ( scaleupWeightsHandle.isValid() ) { for( auto w : *scaleupWeightsHandle ) hWeights_scaleup_->Fill(w); }
  if ( scaledownWeightsHandle.isValid() ) { for( auto w : *scaledownWeightsHandle ) hWeights_scaledown_->Fill(w); }
  if ( pdfWeightsHandle.isValid() ) { for( auto w : *pdfWeightsHandle ) hWeights_pdf_->Fill(w); }
  if ( otherWeightsHandle.isValid() ) { for( auto w : *otherWeightsHandle ) hWeights_others_->Fill(w); }

  hN_All_->Fill(genParticlesHandle->size(), weight);
  std::vector<const reco::GenParticle*> fsParticles, firstCopies, lastCopies;
  for ( const auto& x : *genParticlesHandle ) {
    if ( x.pt() < 1e-3 and std::abs(x.pz()) > 100 ) continue; // Rough selection to skip incident beams and partons
    if ( x.numberOfDaughters() == 0 ) fsParticles.push_back(&x);
    if ( isFirst(x) ) firstCopies.push_back(&x);
    if ( isLast(x) ) lastCopies.push_back(&x);
  }
  int nUnstableFinal = 0; // This should be zero
  for ( const auto x : fsParticles ) {
    if ( x->status() != 1 ) ++nUnstableFinal;
  }
  hN_Final_->Fill(fsParticles.size(), weight);
  hN_UnstableFinal_->Fill(nUnstableFinal, weight);
  hN_FirstCopy_->Fill(firstCopies.size(), weight);
  hN_LastCopy_->Fill(lastCopies.size(), weight);

  int nT = 0, nW = 0, nZ = 0, nG = 0, nL = 0, nH = 0, nB = 0;
  for ( const auto x : firstCopies ) {
    const int pid = x->pdgId();
    const int aid = std::abs(pid);
    const double mass = x->mass();

    if      ( aid ==  6 ) { hM_FirstT_->Fill(mass, weight); ++nT; }
    else if ( aid == 24 ) { hM_FirstW_->Fill(mass, weight); ++nW; }
    else if ( aid == 23 ) { hM_FirstZ_->Fill(mass, weight); ++nZ; }
    else if ( aid == 25 ) { hM_FirstH_->Fill(mass, weight); ++nH; }
    else if ( aid ==  5 ) ++nB;
    else if ( aid == 22 ) ++nG;
    else if ( aid == 11 or aid == 13 ) ++nL;
  }
  hN_FirstT_->Fill(nT, weight);
  hN_FirstW_->Fill(nW, weight);
  hN_FirstZ_->Fill(nZ, weight);
  hN_FirstH_->Fill(nH, weight);
  hN_FirstB_->Fill(nB, weight);
  hN_FirstG_->Fill(nG, weight);
  hN_FirstL_->Fill(nL, weight);

  nT = nW = nZ = nG = nL = nH = nB = 0;
  for ( const auto x : lastCopies ) {
    const int pid = x->pdgId();
    const int aid = std::abs(pid);
    const double mass = x->mass();

    if      ( aid ==  6 ) { hM_LastT_->Fill(mass, weight); ++nT; }
    else if ( aid == 24 ) { hM_LastW_->Fill(mass, weight); ++nW; }
    else if ( aid == 23 ) { hM_LastZ_->Fill(mass, weight); ++nZ; }
    else if ( aid == 25 ) { hM_LastH_->Fill(mass, weight); ++nH; }
    else if ( aid ==  5 ) ++nB;
    else if ( aid == 22 ) ++nG;
    else if ( aid == 11 or aid == 13 ) ++nL;
  }
  hN_LastT_->Fill(nT, weight);
  hN_LastW_->Fill(nW, weight);
  hN_LastZ_->Fill(nZ, weight);
  hN_LastH_->Fill(nH, weight);
  hN_LastB_->Fill(nB, weight);
  hN_LastG_->Fill(nG, weight);
  hN_LastL_->Fill(nL, weight);

}

DEFINE_FWK_MODULE(CATGenValidation);

