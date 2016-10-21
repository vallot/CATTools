#ifndef CATTools_FatJet_H
#define CATTools_FatJet_H

#include "DataFormats/Candidate/interface/CandidateFwd.h"
#include "DataFormats/Candidate/interface/Candidate.h"
#include "CATTools/DataFormats/interface/Particle.h"
#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/PatCandidates/interface/JetCorrFactors.h"
#include "CATTools/DataFormats/interface/Jet.h"

#include <string>
#include <boost/array.hpp>

// Define typedefs for convenience
namespace cat {
    class FatJet;
    typedef std::vector<FatJet>              FatJetCollection;
    typedef edm::Ref<FatJetCollection>       FatJetRef;
    typedef edm::RefVector<FatJetCollection> FatJetRefVector;
}

namespace cat {

  class FatJet : public Jet {
  public:
    FatJet();
    FatJet(const reco::LeafCandidate & aJet);
    virtual ~FatJet();

	// accessing variables for jet substructure
	float tau1() const { return tau1_; }
	float tau2() const { return tau2_; }
	float tau3() const { return tau3_; }
	float prunedmass() const { return prunedMass_; }
	float softdropmass() const { return softdropMass_; }
	float puppi_pt() const { return puppi_pt_; }
	float puppi_eta() const { return puppi_eta_; }
	float puppi_phi() const { return puppi_phi_; }
	float puppi_mass() const { return puppi_mass_; }
	float puppi_tau1() const { return puppi_tau1_; }
	float puppi_tau2() const { return puppi_tau2_; }
	float puppi_tau3() const { return puppi_tau3_; }

	void set_taus(float f1, float f2, float f3) { tau1_ = f1; tau2_ = f2; tau3_ = f3; }
	void set_prunedmass(float f) { prunedMass_ = f; }
	void set_softdropmass(float f) { softdropMass_ = f;}
	void set_puppijet(float fpt, float feta, float fphi, float fm) { puppi_pt_ = fpt; puppi_eta_ = feta; puppi_phi_ = fphi, puppi_mass_ = fm; }
	void set_puppijettaus(float f1, float f2, float f3) { puppi_tau1_ = f1; puppi_tau2_ = f2; puppi_tau3_ = f3; }

  private:

	// Additional variables for fat jets
	float tau1_, tau2_, tau3_;
	float prunedMass_;
	float softdropMass_;

	// associated puppi jet info
	float puppi_pt_, puppi_eta_, puppi_phi_, puppi_mass_;
	float puppi_tau1_, puppi_tau2_, puppi_tau3_;
  };
}

#endif
