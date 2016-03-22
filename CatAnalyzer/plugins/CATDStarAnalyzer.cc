#include "TtbarDiLeptonAnalyzer.h"
#include "CATTools/DataFormats/interface/SecVertex.h"



using namespace std;
using namespace cat;

class CATDStarAnalyzer : public TtbarDiLeptonAnalyzer {
public:
  virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
  explicit CATDStarAnalyzer(const edm::ParameterSet&);
  ~CATDStarAnalyzer();


private:
  //void endLuminosityBlock(const edm::LuminosityBlock&, const edm::EventSetup&) override {};

  enum sys_e {sys_nom,
    sys_jes_u, sys_jes_d, sys_jer_u, sys_jer_d,
    sys_mu_u, sys_mu_d, sys_el_u, sys_el_d,
	      //sys_mueff_u, sys_mueff_d, sys_eleff_u, sys_eleff_d,
	      //sys_btag_u, sys_btag_d,
    nsys_e
  };
  edm::EDGetTokenT<cat::SecVertexCollection>      d0Token_;
  edm::EDGetTokenT<cat::SecVertexCollection>      dstarToken_;
  

  void resetBr();

  std::vector<bool> b_d0_true, b_d0_fit;
  std::vector<float> b_d0_pt, b_d0_eta, b_d0_phi, b_d0_m, b_d0_LXY, b_d0_L3D;
  std::vector<float> b_d0_dau1_pt, b_d0_dau1_eta, b_d0_dau1_phi, b_d0_dau1_m, b_d0_dau1_q;
  std::vector<float> b_d0_dau2_pt, b_d0_dau2_eta, b_d0_dau2_phi, b_d0_dau2_m, b_d0_dau2_q;

  std::vector<bool> b_dstar_true, b_dstar_fit;
  std::vector<float> b_dstar_pt, b_dstar_eta, b_dstar_phi, b_dstar_m, b_dstar_LXY, b_dstar_L3D;
  std::vector<float> b_dstar_dau1_pt, b_dstar_dau1_eta, b_dstar_dau1_phi, b_dstar_dau1_m, b_dstar_dau1_q;
  std::vector<float> b_dstar_dau2_pt, b_dstar_dau2_eta, b_dstar_dau2_phi, b_dstar_dau2_m, b_dstar_dau2_q;
  std::vector<float> b_dstar_dau3_pt, b_dstar_dau3_eta, b_dstar_dau3_phi, b_dstar_dau3_m, b_dstar_dau3_q;

  int NCutflow; 
  std::vector<std::vector<int> > cutflow_;
  std::vector<TTree*> ttree_;

  

};
//
// constructors and destructor
//
CATDStarAnalyzer::CATDStarAnalyzer(const edm::ParameterSet& iConfig) : TtbarDiLeptonAnalyzer( iConfig )
{
  NCutflow = getNCutflow();
  cutflow_ = getCutFlow();
  d0Token_  = consumes<cat::SecVertexCollection>(iConfig.getParameter<edm::InputTag>("d0s"));
  dstarToken_  = consumes<cat::SecVertexCollection>(iConfig.getParameter<edm::InputTag>("dstars"));

  const std::string sys_name[nsys_e] = {
    "nom",
    "jes_u", "jes_d", "jer_u", "jer_d",
    "mu_u", "mu_d", "el_u", "el_d",
    //    "mueff_u", "mueff_d", "eleff_u", "eleff_d",
    //    "btag_u", "btag_d"
  };
  ttree_ = getTree();
  for (int sys = 0; sys < nsys_e; ++sys){
    auto tr = ttree_[sys];
    tr->Branch("d0_pt" ,"std::vector<float>",&b_d0_pt);
    tr->Branch("d0_eta","std::vector<float>",&b_d0_eta);
    tr->Branch("d0_phi","std::vector<float>",&b_d0_phi);
    tr->Branch("d0_m","std::vector<float>",&b_d0_m);
    tr->Branch("d0_true","std::vector<bool>",&b_d0_true);
    tr->Branch("d0_fit","std::vector<bool>",&b_d0_fit);
    tr->Branch("d0_L3D","std::vector<float>",&b_d0_L3D);
    tr->Branch("d0_LXY","std::vector<float>",&b_d0_LXY);

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
}

CATDStarAnalyzer::~CATDStarAnalyzer()
{
  cout <<"     cut flow   emu    ee    mumu"<< endl;
  for ( int i=0; i<NCutflow; ++i ) {
    cout <<"step "<< i << "    "<< cutflow_[i][0] <<  "   "<< cutflow_[i][1] << "   " << cutflow_[i][2] << "   " << cutflow_[i][3]<< endl;
  }
}

void CATDStarAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  const bool runOnMC = !iEvent.isRealData();
  TtbarDiLeptonAnalyzer::analyze( iEvent, iSetup);

  for (int sys = 0; sys < nsys_e; ++sys){
    if (sys > 0 && !runOnMC) break;
    resetBr();

    edm::Handle<cat::SecVertexCollection> d0s;       iEvent.getByToken(d0Token_,d0s);
    edm::Handle<cat::SecVertexCollection> dstars;    iEvent.getByToken(dstarToken_,dstars);


    for( const auto& x : *d0s) {
      b_d0_pt.push_back( x.pt());
      b_d0_eta.push_back( x.eta());
      b_d0_phi.push_back( x.phi());
      b_d0_m.push_back( x.mass());
      if( x.numberOfDaughters()==3) b_d0_true.push_back(true);
      else b_d0_true.push_back(false);

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
      if( x.numberOfDaughters() ==4 ) b_dstar_true.push_back(true);
      else b_dstar_true.push_back(false);

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


    getTree()[sys]->Fill();
  }
}


void CATDStarAnalyzer::resetBr()
{
  TtbarDiLeptonAnalyzer::resetBr();
  b_d0_pt.clear();  b_d0_eta.clear();  b_d0_phi.clear(); b_d0_m.clear(); b_d0_true.clear() ; b_d0_fit.clear();
  b_dstar_pt.clear(); b_dstar_eta.clear() ; b_dstar_phi.clear(); b_dstar_m.clear(); b_dstar_true.clear(), b_dstar_fit.clear(); 

  b_d0_dau1_pt.clear(); b_d0_dau1_eta.clear(); b_d0_dau1_phi.clear(); b_d0_dau1_m.clear(); b_d0_dau1_q.clear();
  b_d0_dau2_pt.clear(); b_d0_dau2_eta.clear(); b_d0_dau2_phi.clear(); b_d0_dau2_m.clear(); b_d0_dau2_q.clear();

  b_dstar_dau1_pt.clear() ; b_dstar_dau1_eta.clear(); b_dstar_dau1_phi.clear(); b_dstar_dau1_m.clear(); b_dstar_dau1_q.clear();
  b_dstar_dau2_pt.clear() ; b_dstar_dau2_eta.clear(); b_dstar_dau2_phi.clear(); b_dstar_dau2_m.clear(); b_dstar_dau2_q.clear();
  b_dstar_dau3_pt.clear() ; b_dstar_dau3_eta.clear(); b_dstar_dau3_phi.clear(); b_dstar_dau3_m.clear(); b_dstar_dau3_q.clear();

}

//define this as a plug-in
DEFINE_FWK_MODULE(CATDStarAnalyzer);
