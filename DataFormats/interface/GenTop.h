#ifndef CATTools_GenTop_H
#define CATTools_GenTop_H

#include "DataFormats/Candidate/interface/CandidateFwd.h"
#include "DataFormats/Candidate/interface/Candidate.h"
#include "DataFormats/PatCandidates/interface/Particle.h"
#include "DataFormats/JetReco/interface/GenJetCollection.h"
#include "DataFormats/JetReco/interface/GenJet.h"
#include "CATTools/DataFormats/interface/Particle.h"
#include "CATTools/DataFormats/interface/MCParticle.h"
#include "CommonTools/Utils/interface/PtComparator.h"
#include "DataFormats/Math/interface/deltaR.h"
#include <map>
#include <vector>

using namespace edm;
using namespace std;
using namespace reco;

// Define typedefs for convenience
namespace cat {
  class GenTop;
  typedef std::vector<GenTop>              GenTopCollection;
  typedef edm::Ref<GenTopCollection>       GenTopRef;
  typedef edm::RefVector<GenTopCollection> GenTopRefVector;
}

namespace cat {

  class GenTop : public reco::LeafCandidate{
  public:
    GenTop();
    /* GenTop(const reco::GenParticle & aGenTop);  */
    GenTop(const reco::Candidate & aGenTop);
    virtual ~GenTop();

    typedef std::vector<math::XYZTLorentzVector> LorentzVectors;

    // status 3
    const math::XYZTLorentzVector topquark1() const { return tops_[0]; }
    const math::XYZTLorentzVector topquark2() const { return tops_[1]; }

    const math::XYZTLorentzVector bquarks1() const { return bquarks_[0]; }
    const math::XYZTLorentzVector bquarks2() const { return bquarks_[1]; }
    const math::XYZTLorentzVector bquarks3() const { return bquarks_[2]; }
    const math::XYZTLorentzVector bquarks4() const { return bquarks_[3]; }

    const math::XYZTLorentzVector lepton1() const { return leptons_[0]; }
    const math::XYZTLorentzVector lepton2() const { return leptons_[1]; }
    const math::XYZTLorentzVector nu1() const { return nus_[0]; }
    const math::XYZTLorentzVector nu2() const { return nus_[1]; }
    const math::XYZTLorentzVector taunu1() const { return taunus_[0]; }
    const math::XYZTLorentzVector taunu2() const { return taunus_[1]; }

    const math::XYZTLorentzVector cJets1() const { return cJets_[0]; }
    const math::XYZTLorentzVector cJets2() const { return cJets_[1]; }

    const math::XYZTLorentzVector bJets1() const { return bJets_[0]; }
    const math::XYZTLorentzVector bJets2() const { return bJets_[1]; }
    const math::XYZTLorentzVector bJets3() const { return bJets_[2]; }
    const math::XYZTLorentzVector bJets4() const { return bJets_[3]; }

    const math::XYZTLorentzVector addbJets1(int i = 0) const {
      if( i == 0) return addbJetsHad_[0];
      else return addbJets_[0];
    }

    const math::XYZTLorentzVector addbJets2(int i = 0) const {
      if( i == 0) return addbJetsHad_[1];
      else return addbJets_[1];
    }

    const math::XYZTLorentzVector addJets1() const { return addJets_[0]; }
    const math::XYZTLorentzVector addJets2() const { return addJets_[1]; }

    //void building( const std::vector<reco::GenJet>* genJets, const std::vector<reco::GenParticle>* genParticles );
    void building( Handle<reco::GenJetCollection> genJets, Handle<reco::GenParticleCollection> genParticles, Handle<std::vector<int> > genBHadFlavour, Handle<std::vector<int> > genBHadJetIndex, Handle<std::vector<int> > genCHadFlavour, Handle<std::vector<int> > genCHadJetIndex  );

    float ttbarmass() const { return ttbarmass_; }

    float dRaddJets() const { return dRaddJets_; }

    float dRaddbJets(int i=0) const {
      if( i == 0){ return dRaddbJetsHad_; }
      else return dRaddbJets_;
    }
    float dRaddcJets(int i=0) const {
      if( i == 0){ return dRaddcJetsHad_; }
      else return dRaddcJets_;
    }
 
    float dRcJets(int i=0) const {
      if( i == 0){ return dRcJetsHad_; }
      else return dRcJets_;
    }

    bool taunic(int i = -1) const {
      bool hasTau = false;
      if( i == -1) hasTau = taunic1_ || taunic2_;
      if( i == 0 ) hasTau = taunic1_;
      if( i == 1 ) hasTau = taunic2_;
      return hasTau;
    }

    bool allHadronic() const { return allHadronic_; }

    bool semiLeptonic(int i = -1) const {
      bool decay = false;
      if( i == -1) decay = semiLeptonic_;
      if( i == 0) decay = semiLeptonicMuo_ || semiLeptonicEle_;
      if( i == 1) decay = ( semiLeptonicMuo_ || semiLeptonicEle_ ) && !semiLeptonicTau_;
      return decay;
    }

    bool semiLeptonicMuo() const { return semiLeptonicMuo_; }
    bool semiLeptonicEle() const { return semiLeptonicEle_; }
    bool semiLeptonicTau() const { return semiLeptonicTau_; }

    bool diLeptonic(int i = -1) const {
      bool decay = false;
      if( i == -1) decay = diLeptonic_;
      if( i == 0) decay = diLeptonicMuoMuo_ || diLeptonicMuoEle_ || diLeptonicEleEle_;
      if( i == 1) decay = ( diLeptonicMuoMuo_ || diLeptonicMuoEle_ || diLeptonicEleEle_) && !( diLeptonicTauMuo_ || diLeptonicTauEle_ || diLeptonicTauTau_);
      return decay;
    }

    bool diLeptonicMuoMuo() const { return diLeptonicMuoMuo_; }
    bool diLeptonicMuoEle() const { return diLeptonicMuoEle_; }
    bool diLeptonicEleEle() const { return diLeptonicEleEle_; }
    bool diLeptonicTauMuo() const { return diLeptonicTauMuo_; }
    bool diLeptonicTauEle() const { return diLeptonicTauEle_; }
    bool diLeptonicTauTau() const { return diLeptonicTauTau_; }

    int NbQuarksTop() const { return NbQuarksTop_ ; }
    int NbQuarksNoTop() const { return NbQuarksNoTop_ ; }
    int NbQuarks() const { return NbQuarks_ ; }
    int NbQuarks20() const { return NbQuarks20_ ; }
    int NbQuarks40() const { return NbQuarks40_ ; }
    int NaddbQuarks20() const { return NaddbQuarks20_ ; }
    int NaddbQuarks40() const { return NaddbQuarks40_ ; }
    int NcQuarks() const { return NcQuarks_; }

    int NbJets(int i=0) const {
      if( i == 0 ) return NbJetsBHad_;
      else return NbJets_ ;
    }
    int NbJets10(int i=0) const {
      if( i == 0 ) return NbJets10BHad_;
      else return NbJets10_ ;
    }
    int NbJets15(int i=0) const {
      if( i == 0 ) return NbJets15BHad_;
      else return NbJets15_ ;
    }
    int NbJets20(int i=0) const {
      if( i == 0 ) return NbJets20BHad_;
      else return NbJets20_ ;
    }
    int NbJets25(int i=0) const {
      if( i == 0 ) return NbJets25BHad_;
      else return NbJets25_ ;
    }
    int NbJets30(int i=0) const {
      if( i == 0 ) return NbJets30BHad_;
      else return NbJets30_ ;
    }
    int NbJets40(int i=0) const {
      if( i == 0 ) return NbJets40BHad_;
      else return NbJets40_ ;
    }

    int NaddbJets(int i=0) const {
      if( i == 0 ) return NaddbJetsBHad_;
      else return NbJetsNoTop_;
    }
    int NaddbJets20(int i=0) const {
      if( i == 0 ) return NaddbJets20BHad_;
      else return NbJets20NoTop_;
    }
    int NaddbJets40(int i=0) const {
      if( i == 0 ) return NaddbJets40BHad_;
      else return NbJets40NoTop_;
    }

    int NaddcJets(int i=0) const {
    //int NaddcJets() const {
      if( i == 0 ) return NaddcJetsCHad_;
      else return NaddcJets_;
    }
    int NaddcJets20(int i=0) const {
    //int NaddcJets20() const {
      if( i == 0 ) return NaddcJets20CHad_;
      else return NaddcJets20_;
    }
    int NaddcJets40(int i=0) const {
    //int NaddcJets40() const {
      if( i == 0 ) return NaddcJets40CHad_;
      else return NaddcJets40_;
    }


    int NcJets(int i=0) const {
      if( i == 0 ) return NcJetsCHad_;
      else return NcJets_ ;
    }
    int NcJets10(int i=0) const {
      if( i == 0 ) return NcJets10CHad_;
      else return NcJets10_ ;
    }
    int NcJets15(int i=0) const {
      if( i == 0 ) return NcJets15CHad_;
      else return NcJets15_ ;
    }
    int NcJets20(int i=0) const {
      if( i == 0 ) return NcJets20CHad_;
      else return NcJets20_ ;
    }
    int NcJets25(int i=0) const {
      if( i == 0 ) return NcJets25CHad_;
      else return NcJets25_ ;
    }
    int NcJets30(int i=0) const {
      if( i == 0 ) return NcJets30CHad_;
      else return NcJets30_ ;
    }
    int NcJets40(int i=0) const {
      if( i == 0 ) return NcJets40CHad_;
      else return NcJets40_ ;
    }

    int NbJetsNoTop() const { return NbJetsNoTop_ ; }
    int NbJets15NoTop() const { return NbJets15NoTop_ ; }
    int NbJets20NoTop() const { return NbJets20NoTop_ ; }
    int NbJets25NoTop() const { return NbJets25NoTop_ ; }
    int NbJets30NoTop() const { return NbJets30NoTop_ ; }

    int NJets() const { return NJets_ ;}
    int NJets10() const { return NJets10_ ;}
    int NJets15() const { return NJets15_ ;}
    int NJets20() const { return NJets20_ ;}
    int NJets25() const { return NJets25_ ;}
    int NJets30() const { return NJets30_ ;}
    int NJets40() const { return NJets40_ ;}

    int NaddJets20() const { return NaddJets20_ ;}
    int NaddJets40() const { return NaddJets40_ ;}

    int is2tops() const { return is2tops_; }

  private:

    std::vector<const reco::Candidate *> getAncestors(const reco::Candidate &c);
    bool hasBottom(const reco::Candidate &c);
    bool hasCharm(const reco::Candidate &c);
    bool decayFromBHadron(const reco::Candidate &c);
    bool decayFromCHadron(const reco::Candidate &c);
    const reco::Candidate* lastBHadron(const reco::Candidate &c);
    const reco::Candidate* lastCHadron(const reco::Candidate &c);
    bool isLastbottom(const reco::GenParticle&);
    bool isLastcharm(const reco::GenParticle&);
    bool isLastParton(const reco::GenParticle&);
    bool isFromtop(const reco::GenParticle&);
    bool isFromW(const reco::GenParticle&);
    const reco::Candidate* getLast( const reco::Candidate& p );
    float deltaR( const reco::Candidate &pasObj, const reco::Candidate &proObj );

    bool is2tops_;

    LorentzVectors tops_;
    LorentzVectors bquarks_;
    LorentzVectors leptons_;
    LorentzVectors nus_;
    LorentzVectors taunus_;
    LorentzVectors cJets_;
    LorentzVectors bJets_;
    LorentzVectors addbJets_;
    LorentzVectors addcJets_;
    LorentzVectors addbJetsHad_;
    LorentzVectors addcJetsHad_;
    LorentzVectors addJets_;
    LorentzVectors quarksfromW_;

    float ttbarmass_;

    bool allHadronic_;
    bool semiLeptonic_;
    bool semiLeptonicMuo_;
    bool semiLeptonicEle_;
    bool semiLeptonicTau_;
    bool diLeptonic_;
    bool diLeptonicMuoMuo_;
    bool diLeptonicMuoEle_;
    bool diLeptonicEleEle_;
    bool diLeptonicTauMuo_;
    bool diLeptonicTauEle_;
    bool diLeptonicTauTau_;

    bool ttbbDecay_;

    bool taunic1_;
    bool taunic2_;

    int NbJets_;
    int NbJets10_;
    int NbJets15_;
    int NbJets20_;
    int NbJets25_;
    int NbJets30_;
    int NbJets40_;

    int NbJetsBHad_;
    int NbJets10BHad_;
    int NbJets15BHad_;
    int NbJets20BHad_;
    int NbJets25BHad_;
    int NbJets30BHad_;
    int NbJets40BHad_;

    int NaddbJetsBHad_;
    int NaddbJets20BHad_;
    int NaddbJets40BHad_;

    int NaddcJetsCHad_;
    int NaddcJets20CHad_;
    int NaddcJets40CHad_;

    int NbJetsNoTop_;
    int NbJets10NoTop_;
    int NbJets15NoTop_;
    int NbJets20NoTop_;
    int NbJets25NoTop_;
    int NbJets30NoTop_;
    int NbJets40NoTop_;


    int NcJets_;
    int NcJets10_;
    int NcJets15_;
    int NcJets20_;
    int NcJets25_;
    int NcJets30_;
    int NcJets40_;

    int NaddcJets_;
    int NaddcJets10_;
    int NaddcJets20_;
    int NaddcJets30_;
    int NaddcJets40_;

    int NcJetsCHad_;
    int NcJets10CHad_;
    int NcJets15CHad_;
    int NcJets20CHad_;
    int NcJets25CHad_;
    int NcJets30CHad_;
    int NcJets40CHad_;

    int NbQuarks_;
    int NbQuarksNoTop_;
    int NbQuarksTop_;
    int NbQuarks20_;
    int NbQuarks40_;
    int NaddbQuarks20_;
    int NaddbQuarks40_;

    int NcQuarks_;
    int NaddcQuarks_;
    //int NaddcQuarks20_;
    //int NaddcQuarks40_;

    int NJets_;
    int NJets10_;
    int NJets15_;
    int NJets20_;
    int NJets25_;
    int NJets30_;
    int NJets40_;

    int NaddJets20_;
    int NaddJets40_;


    float dRaddJets_;

    float dRaddbJets_;
    float dRaddcJets_;
    float dRcJets_;

    float dRaddbJetsHad_;
    float dRaddcJetsHad_;
    float dRcJetsHad_;


  };
}

#endif
