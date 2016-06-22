#include"dileptonCommon.h"
#include "CATTools/DataFormats/interface/SecVertex.h"
#include<TClonesArray.h>

using namespace std;
using namespace cat;
using namespace dileptonCommonGlobal;

class CATDstarAnalyzer : public dileptonCommon {
  public:
    explicit CATDstarAnalyzer(const edm::ParameterSet&);
    ~CATDstarAnalyzer();

  protected:
    shared_ptr<TLorentzVector> mcMatching( vector<TLorentzVector>& , TLorentzVector& ) ; 
    virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
    virtual void analyzeCustom(const edm::Event& iEvent, const edm::EventSetup& iSetup, int sys) final;   // "final" keyword to prevent inherited
    virtual void setBranchCustom(TTree* tr, int sys) final;
    virtual void resetBrCustom() final;
    //void beginLuminosityBlock(const edm::LuminosityBlock& lumi, const edm::EventSetup&) override;
    //void endLuminosityBlock(const edm::LuminosityBlock&, const edm::EventSetup&) override {};

    std::vector<bool> b_d0_true, b_d0_fit;
    //std::vector<float> b_d0_pt, b_d0_eta, b_d0_phi, b_d0_m;
    std::vector<float> b_d0_LXY, b_d0_L3D, b_d0_dRTrue, b_d0_relPtTrue, b_d0_dca;
    std::vector<float> b_d0_dau1_q;
    std::vector<float> b_d0_dau2_q;

    std::vector<bool> b_dstar_true, b_dstar_fit;
    //std::vector<float> b_dstar_pt, b_dstar_eta, b_dstar_phi, b_dstar_m;
    std::vector<float> b_dstar_q, b_dstar_LXY, b_dstar_L3D, b_dstar_dRTrue, b_dstar_relPtTrue, b_dstar_dca, b_dstar_dca2, b_dstar_dca3;
    std::vector<float> b_dstar_dau1_q;
    std::vector<float> b_dstar_dau2_q;
    std::vector<float> b_dstar_dau3_q;

    TClonesArray *b_d0,    *b_d0_dau1,    *b_d0_dau2; 
    TClonesArray *b_dstar, *b_dstar_dau1, *b_dstar_dau2, *b_dstar_dau3; 

    edm::EDGetTokenT<cat::SecVertexCollection>      d0Token_;
    edm::EDGetTokenT<cat::SecVertexCollection>      dstarToken_;
    edm::EDGetTokenT<edm::View<reco::GenParticle> >            mcSrc_;

    double matchingDeltaR_;

};
