#include"dileptonCommon.h"
#include "CATTools/DataFormats/interface/SecVertex.h"


using namespace std;
using namespace cat;
using namespace dileptonCommonGlobal;

class CATDstarAnalyzer : public dileptonCommon {
public:
  explicit CATDstarAnalyzer(const edm::ParameterSet&);
  ~CATDstarAnalyzer();

private:
  virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
  virtual void analyzeCustom(const edm::Event& iEvent, const edm::EventSetup& iSetup, int sys) override;
  virtual void setBranchCustom(TTree* tr, int sys) override;
  virtual void resetBrCustom() override;
  //void beginLuminosityBlock(const edm::LuminosityBlock& lumi, const edm::EventSetup&) override;
  //void endLuminosityBlock(const edm::LuminosityBlock&, const edm::EventSetup&) override {};

  std::vector<bool> b_d0_true, b_d0_fit;
  std::vector<float> b_d0_pt, b_d0_eta, b_d0_phi, b_d0_m, b_d0_LXY, b_d0_L3D, b_d0_dRTrue, b_d0_relPtTrue, b_d0_dca;
  std::vector<float> b_d0_dau1_pt, b_d0_dau1_eta, b_d0_dau1_phi, b_d0_dau1_m, b_d0_dau1_q;
  std::vector<float> b_d0_dau2_pt, b_d0_dau2_eta, b_d0_dau2_phi, b_d0_dau2_m, b_d0_dau2_q;

  std::vector<bool> b_dstar_true, b_dstar_fit;
  std::vector<float> b_dstar_pt, b_dstar_eta, b_dstar_phi, b_dstar_m, b_dstar_LXY, b_dstar_L3D, b_dstar_dRTrue, b_dstar_relPtTrue, b_dstar_dca, b_dstar_dca2, b_dstar_dca3;
  std::vector<float> b_dstar_dau1_pt, b_dstar_dau1_eta, b_dstar_dau1_phi, b_dstar_dau1_m, b_dstar_dau1_q;
  std::vector<float> b_dstar_dau2_pt, b_dstar_dau2_eta, b_dstar_dau2_phi, b_dstar_dau2_m, b_dstar_dau2_q;
  std::vector<float> b_dstar_dau3_pt, b_dstar_dau3_eta, b_dstar_dau3_phi, b_dstar_dau3_m, b_dstar_dau3_q;

  edm::EDGetTokenT<cat::SecVertexCollection>      d0Token_;
  edm::EDGetTokenT<cat::SecVertexCollection>      dstarToken_;


};
//
// constructors and destructor
void CATDstarAnalyzer::setBranchCustom(TTree* tr, int sys) {
    tr->Branch("d0_pt" ,"std::vector<float>",&b_d0_pt);
    tr->Branch("d0_eta","std::vector<float>",&b_d0_eta);
    tr->Branch("d0_phi","std::vector<float>",&b_d0_phi);
    tr->Branch("d0_m","std::vector<float>",&b_d0_m);
    tr->Branch("d0_true","std::vector<bool>",&b_d0_true);
    tr->Branch("d0_fit","std::vector<bool>",&b_d0_fit);
    tr->Branch("d0_L3D","std::vector<float>",&b_d0_L3D);
    tr->Branch("d0_LXY","std::vector<float>",&b_d0_LXY);
    tr->Branch("d0_dRTrue","std::vector<float>",&b_d0_dRTrue);
    tr->Branch("d0_relPtTrue","std::vector<float>",&b_d0_relPtTrue);
    tr->Branch("d0_dca","std::vector<float>",&b_d0_dca);

    tr->Branch("d0_dau1_pt" ,"std::vector<float>",&b_d0_dau1_pt );
    tr->Branch("d0_dau1_eta","std::vector<float>",&b_d0_dau1_eta);
    tr->Branch("d0_dau1_phi","std::vector<float>",&b_d0_dau1_phi);
    tr->Branch("d0_dau1_m","std::vector<float>",&b_d0_dau1_m);
    tr->Branch("d0_dau1_q","std::vector<float>",&b_d0_dau1_q);

    tr->Branch("d0_dau2_pt" ,"std::vector<float>",&b_d0_dau2_pt);
    tr->Branch("d0_dau2_eta","std::vector<float>",&b_d0_dau2_eta);
    tr->Branch("d0_dau2_phi","std::vector<float>",&b_d0_dau2_phi);
    tr->Branch("d0_dau2_m","std::vector<float>",&b_d0_dau2_m);
    tr->Branch("d0_dau2_q","std::vector<float>",&b_d0_dau2_q);
  

    tr->Branch("dstar_pt" ,"std::vector<float>",&b_dstar_pt);
    tr->Branch("dstar_eta","std::vector<float>",&b_dstar_eta);
    tr->Branch("dstar_phi","std::vector<float>",&b_dstar_phi);
    tr->Branch("dstar_m","std::vector<float>",&b_dstar_m);
    tr->Branch("dstar_true","std::vector<bool>",&b_dstar_true);
    tr->Branch("dstar_fit","std::vector<bool>",&b_dstar_fit);
    tr->Branch("dstar_L3D","std::vector<float>",&b_dstar_L3D);
    tr->Branch("dstar_LXY","std::vector<float>",&b_dstar_LXY);
    tr->Branch("dstar_dRTrue","std::vector<float>",&b_dstar_dRTrue);
    tr->Branch("dstar_relPtTrue","std::vector<float>",&b_dstar_relPtTrue);
    tr->Branch("dstar_dca","std::vector<float>",&b_dstar_dca);
    tr->Branch("dstar_dca2","std::vector<float>",&b_dstar_dca2);
    tr->Branch("dstar_dca3","std::vector<float>",&b_dstar_dca3);

    tr->Branch("dstar_dau1_pt" ,"std::vector<float>",&b_dstar_dau1_pt );
    tr->Branch("dstar_dau1_eta","std::vector<float>",&b_dstar_dau1_eta);
    tr->Branch("dstar_dau1_phi","std::vector<float>",&b_dstar_dau1_phi);
    tr->Branch("dstar_dau1_m","std::vector<float>",&b_dstar_dau1_m);
    tr->Branch("dstar_dau1_q","std::vector<float>",&b_dstar_dau1_q);

    tr->Branch("dstar_dau2_pt" ,"std::vector<float>",&b_dstar_dau2_pt );
    tr->Branch("dstar_dau2_eta","std::vector<float>",&b_dstar_dau2_eta);
    tr->Branch("dstar_dau2_phi","std::vector<float>",&b_dstar_dau2_phi);
    tr->Branch("dstar_dau2_m","std::vector<float>",&b_dstar_dau2_m);
    tr->Branch("dstar_dau2_q","std::vector<float>",&b_dstar_dau2_q);

    tr->Branch("dstar_dau3_pt" ,"std::vector<float>",&b_dstar_dau3_pt );
    tr->Branch("dstar_dau3_eta","std::vector<float>",&b_dstar_dau3_eta);
    tr->Branch("dstar_dau3_phi","std::vector<float>",&b_dstar_dau3_phi);
    tr->Branch("dstar_dau3_m","std::vector<float>",&b_dstar_dau3_m);
    tr->Branch("dstar_dau3_q","std::vector<float>",&b_dstar_dau3_q);

}

CATDstarAnalyzer::CATDstarAnalyzer(const edm::ParameterSet& iConfig) : dileptonCommon(iConfig)
{
  parameterInit(iConfig);

  d0Token_  = consumes<cat::SecVertexCollection>(iConfig.getParameter<edm::InputTag>("d0s"));
  dstarToken_  = consumes<cat::SecVertexCollection>(iConfig.getParameter<edm::InputTag>("dstars"));

  for (int sys = 0; sys < nsys_e; ++sys){
    auto tr = ttree_[sys];
    setBranchCustom(tr, sys);
  }
}

CATDstarAnalyzer::~CATDstarAnalyzer()
{

}
void CATDstarAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  b_run = iEvent.id().run();
  b_event = iEvent.id().event();

  const bool runOnMC = !iEvent.isRealData();
  cutflow_[0][0]++;

  for (int sys = 0; sys < nsys_e; ++sys){
    if (sys > 0 && !runOnMC) break;
    resetBr();
    int terminate = eventSelection(iEvent, iSetup, sys);
    if ( terminate == -1 ) continue;
    else if ( terminate == -2 ) return;

    analyzeCustom(iEvent, iSetup, sys);
    ttree_[sys]->Fill();
  }
}


void CATDstarAnalyzer::analyzeCustom(const edm::Event& iEvent, const edm::EventSetup& iSetup, int sys) {
    edm::Handle<cat::SecVertexCollection> d0s;       iEvent.getByToken(d0Token_,d0s);
    edm::Handle<cat::SecVertexCollection> dstars;    iEvent.getByToken(dstarToken_,dstars);

    for( const auto& x : *d0s) {
      b_d0_pt.push_back( x.pt());
      b_d0_eta.push_back( x.eta());
      b_d0_phi.push_back( x.phi());
      b_d0_m.push_back( x.mass());
      b_d0_dca.push_back( x.dca());
      if( x.numberOfDaughters()==3) { 
        b_d0_true.push_back(true);
        auto mcTrue = x.daughter(2);
        float dR = reco::deltaR( *mcTrue, x);
        float relPt = abs( mcTrue->pt()- x.pt())/mcTrue->pt();
        b_d0_dRTrue.push_back(dR);
        b_d0_relPtTrue.push_back(relPt);
      }

      else {
        b_d0_true.push_back(false);
      }

      double d0_vProb = x.vProb();
      if ( abs( d0_vProb ) > 1e-5 ) {
        b_d0_fit.push_back(true);
        b_d0_L3D.push_back( x.l3D());
        b_d0_LXY.push_back( x.lxy());
      }
      else b_d0_fit.push_back(false);

      b_d0_dau1_pt.push_back ( x.daughter(0)->pt());
      b_d0_dau1_eta.push_back( x.daughter(0)->eta());
      b_d0_dau1_phi.push_back( x.daughter(0)->phi());
      b_d0_dau1_m.push_back  ( x.daughter(0)->mass());
      b_d0_dau1_q.push_back  ( x.daughter(0)->charge());

      b_d0_dau2_pt.push_back ( x.daughter(1)->pt());
      b_d0_dau2_eta.push_back( x.daughter(1)->eta());
      b_d0_dau2_phi.push_back( x.daughter(1)->phi());
      b_d0_dau2_m.push_back  ( x.daughter(1)->mass());
      b_d0_dau2_q.push_back  ( x.daughter(1)->charge());
    }
    for( const auto& x : *dstars) {
      //auto& dstar_vertex = x.vertex();
      b_dstar_pt.push_back( x.pt());
      b_dstar_eta.push_back( x.eta());
      b_dstar_phi.push_back( x.phi());
      b_dstar_m.push_back( x.mass());
      b_dstar_dca.push_back( x.dca());
      b_dstar_dca2.push_back( x.dca(1));
      b_dstar_dca3.push_back( x.dca(2));
      if( x.numberOfDaughters() ==4 ) {
        b_dstar_true.push_back(true);
        auto mcTrue = x.daughter(3);
        float dR = reco::deltaR( *mcTrue, x );
        float relPt = abs( mcTrue->pt()- x.pt())/mcTrue->pt();
        b_dstar_dRTrue.push_back(dR);
        b_dstar_relPtTrue.push_back(relPt);
      }
      else {
        b_dstar_true.push_back(false);
        b_dstar_dRTrue.push_back(-9);
        b_dstar_relPtTrue.push_back(-9);
      }

      double dstar_vProb = x.vProb();
      if ( abs( dstar_vProb) > 1e-5) {
        b_dstar_fit.push_back(true);
        b_dstar_L3D.push_back( x.l3D());
        b_dstar_LXY.push_back( x.lxy());
      }
      else b_dstar_fit.push_back(false);


      b_dstar_dau1_pt.push_back ( x.daughter(0)->pt());
      b_dstar_dau1_eta.push_back( x.daughter(0)->eta());
      b_dstar_dau1_phi.push_back( x.daughter(0)->phi());
      b_dstar_dau1_m.push_back  ( x.daughter(0)->mass());
      b_dstar_dau1_q.push_back  ( x.daughter(0)->charge());

      b_dstar_dau2_pt.push_back ( x.daughter(1)->pt());
      b_dstar_dau2_eta.push_back( x.daughter(1)->eta());
      b_dstar_dau2_phi.push_back( x.daughter(1)->phi());
      b_dstar_dau2_m.push_back  ( x.daughter(1)->mass());
      b_dstar_dau2_q.push_back  ( x.daughter(1)->charge());

      b_dstar_dau3_pt.push_back ( x.daughter(2)->pt());
      b_dstar_dau3_eta.push_back( x.daughter(2)->eta());
      b_dstar_dau3_phi.push_back( x.daughter(2)->phi());
      b_dstar_dau3_m.push_back  ( x.daughter(2)->mass());
      b_dstar_dau3_q.push_back  ( x.daughter(2)->charge());
    }

}

void CATDstarAnalyzer::resetBrCustom()
{
  b_d0_pt.clear();  b_d0_eta.clear();  b_d0_phi.clear(); b_d0_m.clear(); b_d0_true.clear() ; 
  b_d0_LXY.clear(); b_d0_L3D.clear(); b_d0_fit.clear(); b_d0_dRTrue.clear(); b_d0_relPtTrue.clear(); b_d0_dca.clear();
  b_dstar_pt.clear(); b_dstar_eta.clear() ; b_dstar_phi.clear(); b_dstar_m.clear(); b_dstar_true.clear(); 
  b_dstar_LXY.clear(); b_dstar_L3D.clear(); b_dstar_fit.clear(); b_dstar_dRTrue.clear(); b_dstar_relPtTrue.clear(); b_dstar_dca.clear(); b_dstar_dca2.clear(); b_dstar_dca3.clear();

  b_d0_dau1_pt.clear(); b_d0_dau1_eta.clear(); b_d0_dau1_phi.clear(); b_d0_dau1_m.clear(); b_d0_dau1_q.clear();
  b_d0_dau2_pt.clear(); b_d0_dau2_eta.clear(); b_d0_dau2_phi.clear(); b_d0_dau2_m.clear(); b_d0_dau2_q.clear();

  b_dstar_dau1_pt.clear() ; b_dstar_dau1_eta.clear(); b_dstar_dau1_phi.clear(); b_dstar_dau1_m.clear(); b_dstar_dau1_q.clear();
  b_dstar_dau2_pt.clear() ; b_dstar_dau2_eta.clear(); b_dstar_dau2_phi.clear(); b_dstar_dau2_m.clear(); b_dstar_dau2_q.clear();
  b_dstar_dau3_pt.clear() ; b_dstar_dau3_eta.clear(); b_dstar_dau3_phi.clear(); b_dstar_dau3_m.clear(); b_dstar_dau3_q.clear();


}

//define this as a plug-in
DEFINE_FWK_MODULE(CATDstarAnalyzer);
