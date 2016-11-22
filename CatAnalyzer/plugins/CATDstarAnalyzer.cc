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
    std::vector<float> b_d0_vProb;

    std::vector<float> b_d0_lepSV_lowM;
    std::vector<float> b_d0_lepSV_dRM;
    std::vector<float> b_d0_lepSV_correctM; // for test

    std::vector<bool> b_dstar_true, b_dstar_fit;
    //std::vector<float> b_dstar_pt, b_dstar_eta, b_dstar_phi, b_dstar_m;
    std::vector<float> b_dstar_q, b_dstar_LXY, b_dstar_L3D, b_dstar_dRTrue, b_dstar_relPtTrue, b_dstar_dca, b_dstar_dca2, b_dstar_dca3;
    std::vector<float> b_dstar_dau1_q;
    std::vector<float> b_dstar_dau2_q;
    std::vector<float> b_dstar_dau3_q;
    std::vector<float> b_dstar_vProb;
    std::vector<float> b_dstar_diffMass;

    std::vector<float> b_dstar_lepSV_lowM;
    std::vector<float> b_dstar_lepSV_dRM;
    std::vector<float> b_dstar_opCharge_M;
    std::vector<float> b_dstar_lepSV_correctM;

    std::vector<bool> b_Jpsi_true, b_Jpsi_fit;
    //std::vector<float> b_Jpsi_pt, b_Jpsi_eta, b_Jpsi_phi, b_Jpsi_m;
    std::vector<float> b_Jpsi_LXY, b_Jpsi_L3D, b_Jpsi_dRTrue, b_Jpsi_relPtTrue, b_Jpsi_dca;
    std::vector<float> b_Jpsi_dau1_q;
    std::vector<float> b_Jpsi_dau2_q;
    std::vector<int>   b_Jpsi_dau_pid;
    std::vector<float> b_Jpsi_vProb;

    std::vector<float> b_Jpsi_lepSV_lowM;
    std::vector<float> b_Jpsi_lepSV_dRM;
    std::vector<float> b_Jpsi_lepSV_correctM; // for test

    TClonesArray *b_d0,    *b_d0_dau1,    *b_d0_dau2; 
    TClonesArray *b_dstar, *b_dstar_dau1, *b_dstar_dau2, *b_dstar_dau3; 
    TClonesArray *b_Jpsi,    *b_Jpsi_dau1,    *b_Jpsi_dau2; 

    edm::EDGetTokenT<cat::SecVertexCollection>      d0Token_;
    edm::EDGetTokenT<cat::SecVertexCollection>      dstarToken_;
    edm::EDGetTokenT<cat::SecVertexCollection>      JpsiToken_;
    edm::EDGetTokenT<edm::View<reco::GenParticle> >            mcSrc_;

    double matchingDeltaR_;

};


int isFromtop( const reco::GenParticle& p){
  int nIDMother;
  
  int nTopMother = 0;

  //  string pt = Form("%f", p.pt());
  //  string pdgid = Form("%i",p.pdgId());
  const reco::GenParticle* mother = dynamic_cast<const reco::GenParticle*>(p.mother());
  while( mother != nullptr ) {
    //    string id = Form("%i", mother->pdgId());
    //    string mopt = Form("%f", mother->pt());
    nIDMother = mother->pdgId();
    
    if( abs(nIDMother) == 6 ) {
      nTopMother = nIDMother;
      break;
    }

    mother = dynamic_cast<const reco::GenParticle*>(mother->mother());
  }

  return nTopMother;
}


//
// constructors and destructor
void CATDstarAnalyzer::setBranchCustom(TTree* tr, int sys) {

  b_d0         = new TClonesArray("TLorentzVector",100);
  b_d0_dau1    = new TClonesArray("TLorentzVector",100);
  b_d0_dau2    = new TClonesArray("TLorentzVector",100);
  b_dstar      = new TClonesArray("TLorentzVector",100); 
  b_dstar_dau1 = new TClonesArray("TLorentzVector",100);
  b_dstar_dau2 = new TClonesArray("TLorentzVector",100);
  b_dstar_dau3 = new TClonesArray("TLorentzVector",100);
  b_Jpsi       = new TClonesArray("TLorentzVector",100);
  b_Jpsi_dau1  = new TClonesArray("TLorentzVector",100);
  b_Jpsi_dau2  = new TClonesArray("TLorentzVector",100);

  // D0
  tr->Branch("d0","TClonesArray",&b_d0,32000,0);
  tr->Branch("d0_dau1","TClonesArray",&b_d0_dau1,32000,0);
  tr->Branch("d0_dau2","TClonesArray",&b_d0_dau2,32000,0);
  tr->Branch("d0_vProb","std::vector<float>",&b_d0_vProb);

  tr->Branch("d0_true","std::vector<bool>",&b_d0_true);
  tr->Branch("d0_fit","std::vector<bool>",&b_d0_fit);

  tr->Branch("d0_L3D","std::vector<float>",&b_d0_L3D);
  tr->Branch("d0_LXY","std::vector<float>",&b_d0_LXY);
  tr->Branch("d0_dRTrue","std::vector<float>",&b_d0_dRTrue);
  tr->Branch("d0_relPtTrue","std::vector<float>",&b_d0_relPtTrue);
  tr->Branch("d0_dca","std::vector<float>",&b_d0_dca);

  tr->Branch("d0_dau1_q","std::vector<float>",&b_d0_dau1_q);
  tr->Branch("d0_dau2_q","std::vector<float>",&b_d0_dau2_q);
  
  tr->Branch("d0_lepSV_lowM","std::vector<float>",&b_d0_lepSV_lowM);
  tr->Branch("d0_lepSV_dRM","std::vector<float>",&b_d0_lepSV_dRM);
  tr->Branch("d0_lepSV_correctM","std::vector<float>",&b_d0_lepSV_correctM); // for test

  // Dstar
  tr->Branch("dstar",    "TClonesArray",&b_dstar    ,32000,0);
  tr->Branch("dstar_dau1","TClonesArray",&b_dstar_dau1,32000,0);
  tr->Branch("dstar_dau2","TClonesArray",&b_dstar_dau2,32000,0);
  tr->Branch("dstar_dau3","TClonesArray",&b_dstar_dau3,32000,0);

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
  tr->Branch("dstar_vProb","std::vector<float>",&b_dstar_vProb);
  tr->Branch("dstar_diffMass","std::vector<float>",&b_dstar_diffMass);
  
  tr->Branch("dstar_lepSV_lowM","std::vector<float>",&b_dstar_lepSV_lowM);
  tr->Branch("dstar_lepSV_dRM","std::vector<float>",&b_dstar_lepSV_dRM);
  tr->Branch("dstar_opCharge_M","std::vector<float>",&b_dstar_opCharge_M);
  tr->Branch("dstar_lepSV_correctM","std::vector<float>",&b_dstar_lepSV_correctM);

  // Jpsi
  tr->Branch("Jpsi","TClonesArray",&b_Jpsi,32000,0);
  tr->Branch("Jpsi_dau1","TClonesArray",&b_Jpsi_dau1,32000,0);
  tr->Branch("Jpsi_dau2","TClonesArray",&b_Jpsi_dau2,32000,0);
  tr->Branch("Jpsi_vProb","std::vector<float>",&b_Jpsi_vProb);

  tr->Branch("Jpsi_true","std::vector<bool>",&b_Jpsi_true);
  tr->Branch("Jpsi_fit","std::vector<bool>",&b_Jpsi_fit);

  tr->Branch("Jpsi_L3D","std::vector<float>",&b_Jpsi_L3D);
  tr->Branch("Jpsi_LXY","std::vector<float>",&b_Jpsi_LXY);
  tr->Branch("Jpsi_dRTrue","std::vector<float>",&b_Jpsi_dRTrue);
  tr->Branch("Jpsi_relPtTrue","std::vector<float>",&b_Jpsi_relPtTrue);
  tr->Branch("Jpsi_dca","std::vector<float>",&b_Jpsi_dca);

  tr->Branch("Jpsi_dau1_q","std::vector<float>",&b_Jpsi_dau1_q);
  tr->Branch("Jpsi_dau2_q","std::vector<float>",&b_Jpsi_dau2_q);
  tr->Branch("Jpsi_dau_pid","std::vector<int>",&b_Jpsi_dau_pid);
  
  tr->Branch("Jpsi_lepSV_lowM","std::vector<float>",&b_Jpsi_lepSV_lowM);
  tr->Branch("Jpsi_lepSV_dRM","std::vector<float>",&b_Jpsi_lepSV_dRM);
  tr->Branch("Jpsi_lepSV_correctM","std::vector<float>",&b_Jpsi_lepSV_correctM); // for test

}

CATDstarAnalyzer::CATDstarAnalyzer(const edm::ParameterSet& iConfig) : dileptonCommon(iConfig)
{
  parameterInit(iConfig);

  d0Token_  = consumes<cat::SecVertexCollection>(iConfig.getParameter<edm::InputTag>("d0s"));
  dstarToken_  = consumes<cat::SecVertexCollection>(iConfig.getParameter<edm::InputTag>("dstars"));
  JpsiToken_  = consumes<cat::SecVertexCollection>(iConfig.getParameter<edm::InputTag>("Jpsis"));
  mcSrc_ = consumes<edm::View<reco::GenParticle> >(iConfig.getParameter<edm::InputTag>("mcLabel"));
  matchingDeltaR_  = iConfig.getParameter<double>("matchingDeltaR");


  for (int sys = 0; sys < nsys_e; ++sys){
    auto tr = ttree_[sys];
    setBranchCustom(tr, sys);
  }
}

CATDstarAnalyzer::~CATDstarAnalyzer()
{

}

shared_ptr<TLorentzVector> CATDstarAnalyzer::mcMatching( vector<TLorentzVector>& aGens, TLorentzVector& aReco) {
  float minDR= 999.;
  //float minRelPt = 1.0;
  shared_ptr<TLorentzVector> matchedGen;
  for( auto& aGen : aGens ) {
    float deltaR = aGen.DeltaR(aReco);
    if ( deltaR < minDR ) { matchedGen=make_shared<TLorentzVector>(aGen); minDR = deltaR; }
  }
  if ( minDR < matchingDeltaR_ ) {
    return matchedGen;
  }
  return nullptr;
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
  const bool runOnMC = !iEvent.isRealData();
  edm::Handle<cat::SecVertexCollection> d0s;       iEvent.getByToken(d0Token_,d0s);
  edm::Handle<cat::SecVertexCollection> dstars;    iEvent.getByToken(dstarToken_,dstars);
  edm::Handle<cat::SecVertexCollection> Jpsis;     iEvent.getByToken(JpsiToken_,Jpsis);
  
  edm::Handle<edm::View<reco::GenParticle> > mcHandle;
  
  if ( runOnMC ) {
    iEvent.getByToken(mcSrc_, mcHandle);
  }

  vector<TLorentzVector> gen_d0s;
  vector<TLorentzVector> gen_dstars;
  vector<TLorentzVector> gen_Jpsis;
  
  int nIDMother;
  
  TLorentzVector vecSumMom;

  if ( runOnMC ) {
    for( const auto& aGenParticle : *mcHandle) {
      //If genParticle is D0,
      if ( std::abs(aGenParticle.pdgId()) == 421 ) {
        gen_d0s.push_back( ToTLorentzVector(aGenParticle));
        
        nIDMother = isFromtop(aGenParticle);
        
        if ( /*0 == 1 &&*/ abs(nIDMother) != 6 ) {
            continue;
        }
        
        //printf("Sign : %i\n", nIDMother * aGenParticle.pdgId() / abs(aGenParticle.pdgId()));
        //printf("%s\n", ( aGenParticle.charge() * aGenParticle.pdgId() >= 0 ? "S" : "O" ));
        
        vecSumMom = ToTLorentzVector(aGenParticle) + 
          ( nIDMother * b_lep1_pid < 0 ? b_lep1 : b_lep2 );
        b_d0_lepSV_correctM.push_back(vecSumMom.M());
      } else if ( std::abs(aGenParticle.pdgId()) ==  413 ) {
        gen_dstars.push_back( ToTLorentzVector(aGenParticle));
        
        nIDMother = isFromtop(aGenParticle);
        
        if ( /*0 == 1 &&*/ abs(nIDMother) != 6 ) {
            continue;
        }
        
        //printf("Sign : %i\n", nIDMother * aGenParticle.pdgId() / abs(aGenParticle.pdgId()));
        //printf("%s\n", ( aGenParticle.charge() * aGenParticle.pdgId() >= 0 ? "S" : "O" ));
        
        vecSumMom = ToTLorentzVector(aGenParticle) + 
          ( nIDMother * b_lep1_pid < 0 ? b_lep1 : b_lep2 );
        b_dstar_lepSV_correctM.push_back(vecSumMom.M());
      } else if ( std::abs(aGenParticle.pdgId()) ==  443 ) {
        gen_Jpsis.push_back( ToTLorentzVector(aGenParticle));
        
        nIDMother = isFromtop(aGenParticle);
        
        if ( /*0 == 1 &&*/ abs(nIDMother) != 6 ) {
            continue;
        }
        
        //printf("Sign : %i\n", nIDMother * aGenParticle.pdgId() / abs(aGenParticle.pdgId()));
        //printf("%s\n", ( aGenParticle.charge() * aGenParticle.pdgId() >= 0 ? "S" : "O" ));
        
        vecSumMom = ToTLorentzVector(aGenParticle) + 
          ( nIDMother * b_lep1_pid < 0 ? b_lep1 : b_lep2 );
        b_Jpsi_lepSV_correctM.push_back(vecSumMom.M());
      }
    }
  } 

  int d0_count=-1;
  int dstar_count=-1;
  int Jpsi_count=-1;

  TClonesArray& br_d0 = *b_d0;
  TClonesArray& br_d0_dau1 = *b_d0_dau1;
  TClonesArray& br_d0_dau2 = *b_d0_dau2;

  TClonesArray& br_dstar = *b_dstar;
  TClonesArray& br_dstar_dau1 = *b_dstar_dau1;
  TClonesArray& br_dstar_dau2 = *b_dstar_dau2;
  TClonesArray& br_dstar_dau3 = *b_dstar_dau3;

  TClonesArray& br_Jpsi = *b_Jpsi;
  TClonesArray& br_Jpsi_dau1 = *b_Jpsi_dau1;
  TClonesArray& br_Jpsi_dau2 = *b_Jpsi_dau2;
  
  TLorentzVector vecDMMom, vecDau12;
  float fQDau, fQDM;
  
  TLorentzVector vecSumDMLep1, vecSumDMLep2;
  float fMDMLep1, fMDMLep2;
  
  float fDeltaEta, fDeltaPhi;
  float fSqrtdRMLep1, fSqrtdRMLep2;

  for( auto& x : *d0s) {
    d0_count++; 

    auto d0_tlv = ToTLorentzVector(x);
    new( br_d0[d0_count]) TLorentzVector(d0_tlv);
    new( br_d0_dau1[d0_count]) TLorentzVector(ToTLorentzVector(*(x.daughter(0)))); 
    new( br_d0_dau2[d0_count]) TLorentzVector(ToTLorentzVector(*(x.daughter(1))));
    
    b_d0_dca.push_back( x.dca());

    double d0_vProb = x.vProb();
    b_d0_vProb.push_back(d0_vProb);

    if ( abs( d0_vProb ) > 1e-5 ) {
      b_d0_fit.push_back(true);
      b_d0_L3D.push_back( x.l3D());
      b_d0_LXY.push_back( x.lxy());
    } else {
      b_d0_fit.push_back(false);
      b_d0_L3D.push_back( -9 );
      b_d0_LXY.push_back( -9 );
    }        

    if ( runOnMC ) {
        shared_ptr<TLorentzVector> genMatched = mcMatching( gen_d0s, d0_tlv ); 
        if ( genMatched != nullptr) {
          b_d0_true.push_back( true );
          b_d0_dRTrue.push_back( genMatched->DeltaR( d0_tlv ));
          b_d0_relPtTrue.push_back( (genMatched->Pt()- d0_tlv.Pt())/genMatched->Pt());
        }
        else {
          b_d0_true.push_back( false );
          b_d0_dRTrue.push_back( -9);
          b_d0_relPtTrue.push_back(-9);
        }
    }

    fQDM = 0;

    fQDau = x.daughter(0)->charge();
    b_dstar_dau1_q.push_back(fQDau);
    fQDM += fQDau;

    fQDau = x.daughter(1)->charge();
    b_dstar_dau2_q.push_back(fQDau);
    fQDM += fQDau;
    
    vecDMMom = ToTLorentzVector(*(x.daughter(0))) + ToTLorentzVector(*(x.daughter(1)));
    
    vecSumDMLep1 = b_lep1 + vecDMMom;
    vecSumDMLep2 = b_lep2 + vecDMMom;
    
    fMDMLep1 = vecSumDMLep1.M();
    fMDMLep2 = vecSumDMLep2.M();
    
    fDeltaEta = b_lep1.Eta() - vecDMMom.Eta();
    fDeltaPhi = b_lep1.Phi() - vecDMMom.Phi();
    fSqrtdRMLep1 = fDeltaEta * fDeltaEta + fDeltaPhi * fDeltaPhi;
    
    fDeltaEta = b_lep2.Eta() - vecDMMom.Eta();
    fDeltaPhi = b_lep2.Phi() - vecDMMom.Phi();
    fSqrtdRMLep2 = fDeltaEta * fDeltaEta + fDeltaPhi * fDeltaPhi;
    
    b_d0_lepSV_lowM.push_back(( fMDMLep1 >= fMDMLep2 ? fMDMLep1 : fMDMLep2 ));
    b_d0_lepSV_dRM.push_back(( fSqrtdRMLep1 >= fSqrtdRMLep2 ? fMDMLep1 : fMDMLep2 ));
  }
  for( auto& x : *dstars) {
    dstar_count++;

    auto dstar_tlv = ToTLorentzVector(x);
    new( br_dstar[dstar_count]) TLorentzVector( dstar_tlv );
    new( br_dstar_dau1[dstar_count]) TLorentzVector(ToTLorentzVector(*(x.daughter(0)))); 
    new( br_dstar_dau2[dstar_count]) TLorentzVector(ToTLorentzVector(*(x.daughter(1))));
    new( br_dstar_dau3[dstar_count]) TLorentzVector(ToTLorentzVector(*(x.daughter(2))));

    b_dstar_dca.push_back( x.dca());
    b_dstar_dca2.push_back( x.dca(1));
    b_dstar_dca3.push_back( x.dca(2));

    double dstar_vProb = x.vProb();
    b_dstar_vProb.push_back(dstar_vProb);
	
    if ( abs( dstar_vProb) > 1e-5) {
      b_dstar_fit.push_back(true);
      b_dstar_L3D.push_back( x.l3D());
      b_dstar_LXY.push_back( x.lxy());
    } else {
      b_dstar_fit.push_back(false);
      b_dstar_L3D.push_back( -9 );
      b_dstar_LXY.push_back( -9 );
    }

    if ( runOnMC ) {
        shared_ptr<TLorentzVector> genMatched = mcMatching( gen_dstars, dstar_tlv ); 
        if ( genMatched != nullptr) {
          b_dstar_true.push_back( true );
          b_dstar_dRTrue.push_back( genMatched->DeltaR( dstar_tlv));
          b_dstar_relPtTrue.push_back( (genMatched->Pt()- dstar_tlv.Pt())/genMatched->Pt());
        }
        else {
          b_dstar_true.push_back( false );
          b_dstar_dRTrue.push_back( -9);
          b_dstar_relPtTrue.push_back(-9);
        }
    }

    fQDM = 0;

    fQDau = x.daughter(0)->charge();
    b_dstar_dau1_q.push_back(fQDau);
    fQDM += fQDau;

    fQDau = x.daughter(1)->charge();
    b_dstar_dau2_q.push_back(fQDau);
    fQDM += fQDau;

    fQDau = x.daughter(2)->charge();
    b_dstar_dau3_q.push_back(fQDau);
    fQDM += fQDau;

    vecDMMom = ToTLorentzVector(*(x.daughter(0))) + 
        ToTLorentzVector(*(x.daughter(1))) + 
        ToTLorentzVector(*(x.daughter(2)));

    vecDau12 = ToTLorentzVector(*(x.daughter(0))) + 
        ToTLorentzVector(*(x.daughter(1)));
	
    b_dstar_diffMass.push_back(vecDMMom.M() - vecDau12.M());
    
    vecSumDMLep1 = b_lep1 + vecDMMom;
    vecSumDMLep2 = b_lep2 + vecDMMom;
    
    fMDMLep1 = vecSumDMLep1.M();
    fMDMLep2 = vecSumDMLep2.M();
    
    fDeltaEta = b_lep1.Eta() - vecDMMom.Eta();
    fDeltaPhi = b_lep1.Phi() - vecDMMom.Phi();
    fSqrtdRMLep1 = fDeltaEta * fDeltaEta + fDeltaPhi * fDeltaPhi;
    
    fDeltaEta = b_lep2.Eta() - vecDMMom.Eta();
    fDeltaPhi = b_lep2.Phi() - vecDMMom.Phi();
    fSqrtdRMLep2 = fDeltaEta * fDeltaEta + fDeltaPhi * fDeltaPhi;
    
    b_dstar_lepSV_lowM.push_back(( fMDMLep1 >= fMDMLep2 ? fMDMLep1 : fMDMLep2 ));
    b_dstar_lepSV_dRM.push_back(( fSqrtdRMLep1 >= fSqrtdRMLep2 ? fMDMLep1 : fMDMLep2 ));
    b_dstar_opCharge_M.push_back(( fQDM * b_lep1_pid <= 0.0 ? fMDMLep1 : fMDMLep2 ));
  }
  for( auto& x : *Jpsis) {
    Jpsi_count++;

    auto Jpsi_tlv = ToTLorentzVector(x);
    new( br_Jpsi[Jpsi_count]) TLorentzVector( Jpsi_tlv );
    new( br_Jpsi_dau1[Jpsi_count]) TLorentzVector(ToTLorentzVector(*(x.daughter(0)))); 
    new( br_Jpsi_dau2[Jpsi_count]) TLorentzVector(ToTLorentzVector(*(x.daughter(1))));

    b_Jpsi_dca.push_back( x.dca());

    double Jpsi_vProb = x.vProb();
    b_Jpsi_vProb.push_back(Jpsi_vProb);
	
    if ( abs( Jpsi_vProb) > 1e-5) {
      b_Jpsi_fit.push_back(true);
      b_Jpsi_L3D.push_back( x.l3D());
      b_Jpsi_LXY.push_back( x.lxy());
    } else {
      b_Jpsi_fit.push_back(false);
      b_Jpsi_L3D.push_back( -9 );
      b_Jpsi_LXY.push_back( -9 );
    }

    if ( runOnMC ) {
        shared_ptr<TLorentzVector> genMatched = mcMatching( gen_Jpsis, Jpsi_tlv ); 
        if ( genMatched != nullptr) {
          b_Jpsi_true.push_back( true );
          b_Jpsi_dRTrue.push_back( genMatched->DeltaR( Jpsi_tlv));
          b_Jpsi_relPtTrue.push_back( (genMatched->Pt()- Jpsi_tlv.Pt())/genMatched->Pt());
        }
        else {
          b_Jpsi_true.push_back( false );
          b_Jpsi_dRTrue.push_back( -9);
          b_Jpsi_relPtTrue.push_back(-9);
        }
    }

    fQDM = 0;

    fQDau = x.daughter(0)->charge();
    b_Jpsi_dau1_q.push_back(fQDau);
    fQDM += fQDau;

    fQDau = x.daughter(1)->charge();
    b_Jpsi_dau2_q.push_back(fQDau);
    fQDM += fQDau;
    
    b_Jpsi_dau_pid.push_back(x.daughter(0)->pdgId());
    
    vecDMMom = ToTLorentzVector(*(x.daughter(0))) + ToTLorentzVector(*(x.daughter(1)));
    
    vecSumDMLep1 = b_lep1 + vecDMMom;
    vecSumDMLep2 = b_lep2 + vecDMMom;
    
    fMDMLep1 = vecSumDMLep1.M();
    fMDMLep2 = vecSumDMLep2.M();
    
    fDeltaEta = b_lep1.Eta() - vecDMMom.Eta();
    fDeltaPhi = b_lep1.Phi() - vecDMMom.Phi();
    fSqrtdRMLep1 = fDeltaEta * fDeltaEta + fDeltaPhi * fDeltaPhi;
    
    fDeltaEta = b_lep2.Eta() - vecDMMom.Eta();
    fDeltaPhi = b_lep2.Phi() - vecDMMom.Phi();
    fSqrtdRMLep2 = fDeltaEta * fDeltaEta + fDeltaPhi * fDeltaPhi;
    
    b_Jpsi_lepSV_lowM.push_back(( fMDMLep1 >= fMDMLep2 ? fMDMLep1 : fMDMLep2 ));
    b_Jpsi_lepSV_dRM.push_back(( fSqrtdRMLep1 >= fSqrtdRMLep2 ? fMDMLep1 : fMDMLep2 ));
  }
  
  

}

void CATDstarAnalyzer::resetBrCustom()
{
  b_d0->Clear();    b_d0_dau1->Clear();    b_d0_dau2->Clear();
  b_dstar->Clear(); b_dstar_dau1->Clear(); b_dstar_dau2->Clear(); b_dstar_dau3->Clear();
  b_Jpsi->Clear();    b_Jpsi_dau1->Clear();    b_Jpsi_dau2->Clear();


  // D0
  b_d0_true.clear() ; 
  b_d0_LXY.clear(); b_d0_L3D.clear(); b_d0_fit.clear(); b_d0_dRTrue.clear(); b_d0_relPtTrue.clear(); b_d0_dca.clear();

  b_d0_dau1_q.clear();
  b_d0_dau2_q.clear();
  b_d0_vProb.clear();

  b_d0_lepSV_lowM.clear();
  b_d0_lepSV_dRM.clear();
  b_d0_lepSV_correctM.clear();
  
  // Dstar
  b_dstar_true.clear(); 
  b_dstar_LXY.clear(); b_dstar_L3D.clear(); b_dstar_fit.clear(); b_dstar_dRTrue.clear(); b_dstar_relPtTrue.clear(); b_dstar_dca.clear(); b_dstar_dca2.clear(); b_dstar_dca3.clear();

  b_dstar_dau1_q.clear();
  b_dstar_dau2_q.clear();
  b_dstar_dau3_q.clear();
  b_dstar_vProb.clear();
  b_dstar_diffMass.clear();

  b_dstar_lepSV_lowM.clear();
  b_dstar_lepSV_dRM.clear();
  b_dstar_opCharge_M.clear();
  b_dstar_lepSV_correctM.clear();
  
  // Jpsi
  b_Jpsi_true.clear() ; 
  b_Jpsi_LXY.clear(); b_Jpsi_L3D.clear(); b_Jpsi_fit.clear(); b_Jpsi_dRTrue.clear(); b_Jpsi_relPtTrue.clear(); b_Jpsi_dca.clear();

  b_Jpsi_dau1_q.clear();
  b_Jpsi_dau2_q.clear();
  b_Jpsi_dau_pid.clear();
  b_Jpsi_vProb.clear();

  b_Jpsi_lepSV_lowM.clear();
  b_Jpsi_lepSV_dRM.clear();
  b_Jpsi_lepSV_correctM.clear();

}

//define this as a plug-in
DEFINE_FWK_MODULE(CATDstarAnalyzer);


