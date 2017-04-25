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

    struct CHKeys {
      enum {
        allHadronic = 0,
        semiLeptonic, semiLeptonicMuo, semiLeptonicEle, semiLeptonicTau,
        diLeptonic, diLeptonicMuoMuo, diLeptonicMuoEle, diLeptonicEleEle,
        diLeptonicTauMuo, diLeptonicTauEle, diLeptonicTauTau, 
        ttbbDecay,
        taunic1, taunic2,
        SIZE
      };
    };

    // status 3
    const math::XYZTLorentzVector topquark1() const { return tops_[0]; }
    const math::XYZTLorentzVector topquark2() const { return tops_[1]; }

    const math::XYZTLorentzVector Wquark1() const { return quarksfromW_[0]; }
    const math::XYZTLorentzVector Wquark2() const { return quarksfromW_[1]; }
    const math::XYZTLorentzVector Wquark3() const { return quarksfromW_[2]; }
    const math::XYZTLorentzVector Wquark4() const { return quarksfromW_[3]; }

    const int Wquarkflav1() const { return qflavourfromW_[0]; }
    const int Wquarkflav2() const { return qflavourfromW_[1]; }
    const int Wquarkflav3() const { return qflavourfromW_[2]; }
    const int Wquarkflav4() const { return qflavourfromW_[3]; }
 
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

    const math::XYZTLorentzVector bJetsFromTop1() const { return bJetsFromTop_[0]; }
    const math::XYZTLorentzVector bJetsFromTop2() const { return bJetsFromTop_[1]; }

    const math::XYZTLorentzVector JetsFromW1() const { return JetsFromW_[0]; }
    const math::XYZTLorentzVector JetsFromW2() const { return JetsFromW_[1]; }
    const math::XYZTLorentzVector JetsFromW3() const { return JetsFromW_[2]; }
    const math::XYZTLorentzVector JetsFromW4() const { return JetsFromW_[3]; }

    const int JetsFlavourFromW1() const { return JetsFlavourFromW_[0]; }
    const int JetsFlavourFromW2() const { return JetsFlavourFromW_[1]; }
    const int JetsFlavourFromW3() const { return JetsFlavourFromW_[2]; }
    const int JetsFlavourFromW4() const { return JetsFlavourFromW_[3]; }

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
      if( i == -1) hasTau = chBit_[CHKeys::taunic1] || chBit_[CHKeys::taunic2];
      if( i == 0 ) hasTau = chBit_[CHKeys::taunic1];
      if( i == 1 ) hasTau = chBit_[CHKeys::taunic2];
      return hasTau;
    }

    bool allHadronic() const { return chBit_[CHKeys::allHadronic]; }

    bool semiLeptonic(int i = -1) const {
      bool decay = false;
      if( i == -1) decay = chBit_[CHKeys::semiLeptonic];
      if( i == 0) decay = chBit_[CHKeys::semiLeptonicMuo] || chBit_[CHKeys::semiLeptonicEle];
      if( i == 1) decay = ( chBit_[CHKeys::semiLeptonicMuo] || chBit_[CHKeys::semiLeptonicEle] ) && !chBit_[CHKeys::semiLeptonicTau];
      return decay;
    }

    bool semiLeptonicMuo() const { return chBit_[CHKeys::semiLeptonicMuo]; }
    bool semiLeptonicEle() const { return chBit_[CHKeys::semiLeptonicEle]; }
    bool semiLeptonicTau() const { return chBit_[CHKeys::semiLeptonicTau]; }

    bool diLeptonic(int i = -1) const {
      bool decay = false;
      if( i == -1) decay = chBit_[CHKeys::diLeptonic];
      if( i == 0) decay = chBit_[CHKeys::diLeptonicMuoMuo] || chBit_[CHKeys::diLeptonicMuoEle] || chBit_[CHKeys::diLeptonicEleEle];
      if( i == 1) decay = ( chBit_[CHKeys::diLeptonicMuoMuo] || chBit_[CHKeys::diLeptonicMuoEle] || chBit_[CHKeys::diLeptonicEleEle]) && !( chBit_[CHKeys::diLeptonicTauMuo] || chBit_[CHKeys::diLeptonicTauEle] || chBit_[CHKeys::diLeptonicTauTau]);
      return decay;
    }

    bool diLeptonicMuoMuo() const { return chBit_[CHKeys::diLeptonicMuoMuo]; }
    bool diLeptonicMuoEle() const { return chBit_[CHKeys::diLeptonicMuoEle]; }
    bool diLeptonicEleEle() const { return chBit_[CHKeys::diLeptonicEleEle]; }
    bool diLeptonicTauMuo() const { return chBit_[CHKeys::diLeptonicTauMuo]; }
    bool diLeptonicTauEle() const { return chBit_[CHKeys::diLeptonicTauEle]; }
    bool diLeptonicTauTau() const { return chBit_[CHKeys::diLeptonicTauTau]; }

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

    int NWJets() const {return NWJets_;} 
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
    LorentzVectors bJetsFromTop_;
    LorentzVectors JetsFromW_;
    std::vector<short> JetsFlavourFromW_;
    LorentzVectors addbJets_;
    LorentzVectors addcJets_;
    LorentzVectors addbJetsHad_;
    LorentzVectors addcJetsHad_;
    LorentzVectors addJets_;
    LorentzVectors quarksfromW_;
    std::vector<short> qflavourfromW_;

    float ttbarmass_;

    std::bitset<CHKeys::SIZE> chBit_;

    typedef unsigned char uchar;
    uchar NbJets_;
    uchar NbJets10_;
    uchar NbJets15_;
    uchar NbJets20_;
    uchar NbJets25_;
    uchar NbJets30_;
    uchar NbJets40_;

    uchar NbJetsBHad_;
    uchar NbJets10BHad_;
    uchar NbJets15BHad_;
    uchar NbJets20BHad_;
    uchar NbJets25BHad_;
    uchar NbJets30BHad_;
    uchar NbJets40BHad_;

    uchar NaddbJetsBHad_;
    uchar NaddbJets20BHad_;
    uchar NaddbJets40BHad_;

    uchar NaddcJetsCHad_;
    uchar NaddcJets20CHad_;
    uchar NaddcJets40CHad_;

    uchar NbJetsNoTop_;
    uchar NbJets10NoTop_;
    uchar NbJets15NoTop_;
    uchar NbJets20NoTop_;
    uchar NbJets25NoTop_;
    uchar NbJets30NoTop_;
    uchar NbJets40NoTop_;


    uchar NcJets_;
    uchar NcJets10_;
    uchar NcJets15_;
    uchar NcJets20_;
    uchar NcJets25_;
    uchar NcJets30_;
    uchar NcJets40_;

    uchar NaddcJets_;
    uchar NaddcJets10_;
    uchar NaddcJets20_;
    uchar NaddcJets30_;
    uchar NaddcJets40_;

    uchar NcJetsCHad_;
    uchar NcJets10CHad_;
    uchar NcJets15CHad_;
    uchar NcJets20CHad_;
    uchar NcJets25CHad_;
    uchar NcJets30CHad_;
    uchar NcJets40CHad_;

    uchar NbQuarks_;
    uchar NbQuarksNoTop_;
    uchar NbQuarksTop_;
    uchar NbQuarks20_;
    uchar NbQuarks40_;
    uchar NaddbQuarks20_;
    uchar NaddbQuarks40_;

    uchar NcQuarks_;
    uchar NaddcQuarks_;
    //uchar NaddcQuarks20_;
    //uchar NaddcQuarks40_;

    uchar NJets_;
    uchar NJets10_;
    uchar NJets15_;
    uchar NJets20_;
    uchar NJets25_;
    uchar NJets30_;
    uchar NJets40_;

    uchar NaddJets20_;
    uchar NaddJets40_;

    uchar NWJets_;

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
