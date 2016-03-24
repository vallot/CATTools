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

private:
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


};
//
// constructors and destructor
void CATDstarAnalyzer::setBranchCustom(TTree* tr, int sys) {

    b_d0       = new TClonesArray("TLorentzVector",100);
    b_d0_dau1   = new TClonesArray("TLorentzVector",100);
    b_d0_dau2   = new TClonesArray("TLorentzVector",100);
    b_dstar = new TClonesArray("TLorentzVector",100); 
    b_dstar_dau1   = new TClonesArray("TLorentzVector",100);
    b_dstar_dau2   = new TClonesArray("TLorentzVector",100);
    b_dstar_dau3   = new TClonesArray("TLorentzVector",100);

    tr->Branch("d0","TClonesArray",&b_d0,32000,0);
    tr->Branch("d0_dau1","TClonesArray",&b_d0_dau1,32000,0);
    tr->Branch("d0_dau2","TClonesArray",&b_d0_dau2,32000,0);

    tr->Branch("d0_true","std::vector<bool>",&b_d0_true);
    tr->Branch("d0_fit","std::vector<bool>",&b_d0_fit);

    tr->Branch("d0_L3D","std::vector<float>",&b_d0_L3D);
    tr->Branch("d0_LXY","std::vector<float>",&b_d0_LXY);
    tr->Branch("d0_dRTrue","std::vector<float>",&b_d0_dRTrue);
    tr->Branch("d0_relPtTrue","std::vector<float>",&b_d0_relPtTrue);
    tr->Branch("d0_dca","std::vector<float>",&b_d0_dca);

    tr->Branch("d0_dau1_q","std::vector<float>",&b_d0_dau1_q);

    tr->Branch("d0_dau2_q","std::vector<float>",&b_d0_dau2_q);
  

    tr->Branch("dstar",    "TClonesArray",&b_dstar    ,32000,0);
    tr->Branch("dstar_dau1","TClonesArray",&b_dstar_dau1,32000,0);
    tr->Branch("dstar_dau2","TClonesArray",&b_dstar_dau2,32000,0);
    tr->Branch("dstar_dau3","TClonesArray",&b_dstar_dau3,32000,0);
    /*
    tr->Branch("dstar_pt" ,"std::vector<float>",&b_dstar_pt);
    tr->Branch("dstar_eta","std::vector<float>",&b_dstar_eta);
    tr->Branch("dstar_phi","std::vector<float>",&b_dstar_phi);
    tr->Branch("dstar_m","std::vector<float>",&b_dstar_m);
    */
    tr->Branch("dstar_true","std::vector<bool>",&b_dstar_true);
    tr->Branch("dstar_fit","std::vector<bool>",&b_dstar_fit);
    tr->Branch("dstar_L3D","std::vector<float>",&b_dstar_L3D);
    tr->Branch("dstar_LXY","std::vector<float>",&b_dstar_LXY);
    tr->Branch("dstar_dRTrue","std::vector<float>",&b_dstar_dRTrue);
    tr->Branch("dstar_relPtTrue","std::vector<float>",&b_dstar_relPtTrue);
    tr->Branch("dstar_dca","std::vector<float>",&b_dstar_dca);
    tr->Branch("dstar_dca2","std::vector<float>",&b_dstar_dca2);
    tr->Branch("dstar_dca3","std::vector<float>",&b_dstar_dca3);

    tr->Branch("dstar_dau1_q","std::vector<float>",&b_dstar_dau1_q);

    tr->Branch("dstar_dau2_q","std::vector<float>",&b_dstar_dau2_q);

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

    

    int d0_count=-1;
    int dstar_count=-1;

    TClonesArray& br_d0 = *b_d0;
    TClonesArray& br_d0_dau1 = *b_d0_dau1;
    TClonesArray& br_d0_dau2 = *b_d0_dau2;

    TClonesArray& br_dstar = *b_dstar;
    TClonesArray& br_dstar_dau1 = *b_dstar_dau1;
    TClonesArray& br_dstar_dau2 = *b_dstar_dau2;
    TClonesArray& br_dstar_dau3 = *b_dstar_dau3;

    for( const auto& x : *d0s) {
      d0_count++; 


      new( br_d0[d0_count]) TLorentzVector(ToTLorentzVector(x));
      new( br_d0_dau1[d0_count]) TLorentzVector(ToTLorentzVector(*(x.daughter(0)))); 
      new( br_d0_dau2[d0_count]) TLorentzVector(ToTLorentzVector(*(x.daughter(1))));
      b_d0_dca.push_back( x.dca());

      double d0_vProb = x.vProb();
      if ( abs( d0_vProb ) > 1e-5 ) {
        b_d0_fit.push_back(true);
        b_d0_L3D.push_back( x.l3D());
        b_d0_LXY.push_back( x.lxy());
      }
      else {
        b_d0_fit.push_back(false);
        b_d0_L3D.push_back( -9 );
        b_d0_LXY.push_back( -9 );
      }        
      
 
      b_d0_dau1_q.push_back  ( x.daughter(0)->charge());
      b_d0_dau2_q.push_back  ( x.daughter(1)->charge());
    }
    for( const auto& x : *dstars) {
      dstar_count++;
      new( br_dstar[dstar_count]) TLorentzVector(ToTLorentzVector(x));
      new( br_dstar_dau1[dstar_count]) TLorentzVector(ToTLorentzVector(*(x.daughter(0)))); 
      new( br_dstar_dau2[dstar_count]) TLorentzVector(ToTLorentzVector(*(x.daughter(1))));
      new( br_dstar_dau3[dstar_count]) TLorentzVector(ToTLorentzVector(*(x.daughter(2))));
 
      b_dstar_dca.push_back( x.dca());
      b_dstar_dca2.push_back( x.dca(1));
      b_dstar_dca3.push_back( x.dca(2));

      double dstar_vProb = x.vProb();
      if ( abs( dstar_vProb) > 1e-5) {
        b_dstar_fit.push_back(true);
        b_dstar_L3D.push_back( x.l3D());
        b_dstar_LXY.push_back( x.lxy());
      }
      else {
        b_dstar_fit.push_back(false);
        b_dstar_L3D.push_back( -9 );
        b_dstar_LXY.push_back( -9 );
      }


      b_dstar_dau1_q.push_back  ( x.daughter(0)->charge());

      b_dstar_dau2_q.push_back  ( x.daughter(1)->charge());

      b_dstar_dau3_q.push_back  ( x.daughter(2)->charge());
    }

}

void CATDstarAnalyzer::resetBrCustom()
{
  b_d0->Clear();    b_d0_dau1->Clear();    b_d0_dau2->Clear();
  b_dstar->Clear(); b_dstar_dau1->Clear(); b_dstar_dau2->Clear(); b_dstar_dau3->Clear();


  b_d0_true.clear() ; 
  b_d0_LXY.clear(); b_d0_L3D.clear(); b_d0_fit.clear(); b_d0_dRTrue.clear(); b_d0_relPtTrue.clear(); b_d0_dca.clear();
  b_dstar_true.clear(); 
  b_dstar_LXY.clear(); b_dstar_L3D.clear(); b_dstar_fit.clear(); b_dstar_dRTrue.clear(); b_dstar_relPtTrue.clear(); b_dstar_dca.clear(); b_dstar_dca2.clear(); b_dstar_dca3.clear();

  b_d0_dau1_q.clear();
  b_d0_dau2_q.clear();

  b_dstar_dau1_q.clear();
  b_dstar_dau2_q.clear();
  b_dstar_dau3_q.clear();


}

//define this as a plug-in
DEFINE_FWK_MODULE(CATDstarAnalyzer);
