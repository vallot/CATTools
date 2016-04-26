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
    reco::GenParticle* mcMatching( vector<reco::GenParticle>& ,  const reco::Candidate* ) ; 
    virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
    virtual void analyzeCustom(const edm::Event& iEvent, const edm::EventSetup& iSetup, int sys) final;   // "final" keyword to prevent inherited
    virtual void setBranchCustom(TTree* tr, int sys) final;
    virtual void resetBrCustom() final;
    virtual int isFromB(const reco::Candidate* p);
    virtual int isFromTop(const reco::Candidate*  p);
    //void beginLuminosityBlock(const edm::LuminosityBlock& lumi, const edm::EventSetup&) override;
    //void endLuminosityBlock(const edm::LuminosityBlock&, const edm::EventSetup&) override {};

    std::vector<bool> b_d0_true, b_d0_fit;
    std::vector<int> b_d0_trackQuality, b_d0_isFromB, b_d0_isFromTop;
    std::vector<float> b_d0_LXY, b_d0_L3D, b_d0_vProb, b_d0_dRTrue, b_d0_relPtTrue, b_d0_dca;
    std::vector<float> b_d0_dau1_q;
    std::vector<float> b_d0_dau2_q;

    std::vector<float> b_d0_lepSV_lowM1, b_d0_lepSV_dRM1, b_d0_lepSV_correctM;

    std::vector<bool> b_dstar_true, b_dstar_fit;
    std::vector<int> b_dstar_trackQuality, b_dstar_isFromB, b_dstar_isFromTop;
    std::vector<float> b_dstar_q, b_dstar_LXY, b_dstar_L3D, b_dstar_vProb, b_dstar_dRTrue, b_dstar_relPtTrue, b_dstar_dca, b_dstar_dca2, b_dstar_dca3, b_dstar_diffMass;
    std::vector<float> b_dstar_dau1_q;
    std::vector<float> b_dstar_dau2_q;
    std::vector<float> b_dstar_dau3_q;

    std::vector<float> b_dstar_lepSV_lowM1, b_dstar_lepSV_lowM2, b_dstar_lepSV_dRM1, b_dstar_lepSV_dRM2, b_dstar_lepSV_correctM, b_dstar_lepSV_wrongM, b_dstar_opCharge_M;

    TClonesArray *b_d0,    *b_d0_dau1,    *b_d0_dau2; 
    TClonesArray *b_dstar, *b_dstar_dau1, *b_dstar_dau2, *b_dstar_dau3; 

    edm::EDGetTokenT<cat::SecVertexCollection>      d0Token_;
    edm::EDGetTokenT<cat::SecVertexCollection>      dstarToken_;
    edm::EDGetTokenT<edm::View<reco::GenParticle> >            mcSrc_;

    double matchingDeltaR_;

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
  tr->Branch("d0_vProb","std::vector<float>",&b_d0_vProb);
  tr->Branch("d0_trackQuality","std::vector<int>",&b_d0_trackQuality);

  tr->Branch("d0_dRTrue","std::vector<float>",&b_d0_dRTrue);
  tr->Branch("d0_relPtTrue","std::vector<float>",&b_d0_relPtTrue);
  tr->Branch("d0_dca","std::vector<float>",&b_d0_dca);

  tr->Branch("d0_dau1_q","std::vector<float>",&b_d0_dau1_q);

  tr->Branch("d0_dau2_q","std::vector<float>",&b_d0_dau2_q);
  tr->Branch("d0_isFromB","std::vector<int>",&b_d0_isFromB);
  tr->Branch("d0_isFromTop","std::vector<int>",&b_d0_isFromTop);
  
  tr->Branch("d0_lepSV_lowM1","std::vector<float>",&b_d0_lepSV_lowM1);
  tr->Branch("d0_lepSV_dRM1","std::vector<float>",&b_d0_lepSV_dRM1);
  tr->Branch("d0_lepSV_correctM","std::vector<float>",&b_d0_lepSV_correctM);


  tr->Branch("dstar",    "TClonesArray",&b_dstar    ,32000,0);
  tr->Branch("dstar_dau1","TClonesArray",&b_dstar_dau1,32000,0);
  tr->Branch("dstar_dau2","TClonesArray",&b_dstar_dau2,32000,0);
  tr->Branch("dstar_dau3","TClonesArray",&b_dstar_dau3,32000,0);

  tr->Branch("dstar_true","std::vector<bool>",&b_dstar_true);
  tr->Branch("dstar_fit","std::vector<bool>",&b_dstar_fit);
  tr->Branch("dstar_L3D","std::vector<float>",&b_dstar_L3D);
  tr->Branch("dstar_LXY","std::vector<float>",&b_dstar_LXY);
  tr->Branch("dstar_vProb","std::vector<float>",&b_dstar_vProb);
  tr->Branch("dstar_trackQuality","std::vector<int>",&b_dstar_trackQuality);

  tr->Branch("dstar_diffMass","std::vector<float>",&b_dstar_diffMass);

  tr->Branch("dstar_dRTrue","std::vector<float>",&b_dstar_dRTrue);
  tr->Branch("dstar_relPtTrue","std::vector<float>",&b_dstar_relPtTrue);

  tr->Branch("dstar_dca","std::vector<float>",&b_dstar_dca);
  tr->Branch("dstar_dca2","std::vector<float>",&b_dstar_dca2);
  tr->Branch("dstar_dca3","std::vector<float>",&b_dstar_dca3);

  tr->Branch("dstar_dau1_q","std::vector<float>",&b_dstar_dau1_q);
  tr->Branch("dstar_dau2_q","std::vector<float>",&b_dstar_dau2_q);
  tr->Branch("dstar_dau3_q","std::vector<float>",&b_dstar_dau3_q);

  tr->Branch("dstar_lepSV_lowM1","std::vector<float>",&b_dstar_lepSV_lowM1);
  tr->Branch("dstar_lepSV_lowM2","std::vector<float>",&b_dstar_lepSV_lowM2);
  tr->Branch("dstar_lepSV_dRM1","std::vector<float>",&b_dstar_lepSV_dRM1);
  tr->Branch("dstar_lepSV_dRM2","std::vector<float>",&b_dstar_lepSV_dRM2);
  tr->Branch("dstar_opCharge_M","std::vector<float>",&b_dstar_opCharge_M);

  tr->Branch("dstar_lepSV_correctM","std::vector<float>",&b_dstar_lepSV_correctM);
  tr->Branch("dstar_lepSV_wrongM","std::vector<float>",&b_dstar_lepSV_wrongM);

  tr->Branch("dstar_isFromB","std::vector<int>",&b_dstar_isFromB);
  tr->Branch("dstar_isFromTop","std::vector<int>",&b_dstar_isFromTop);

}

CATDstarAnalyzer::CATDstarAnalyzer(const edm::ParameterSet& iConfig) : dileptonCommon(iConfig)
{
  parameterInit(iConfig);

  d0Token_  = consumes<cat::SecVertexCollection>(iConfig.getParameter<edm::InputTag>("d0s"));
  dstarToken_  = consumes<cat::SecVertexCollection>(iConfig.getParameter<edm::InputTag>("dstars"));
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

reco::GenParticle* CATDstarAnalyzer::mcMatching( vector<reco::GenParticle>& aGens, const reco::Candidate* aReco) {
  float minDR= 999.;
  //float minRelPt = 1.0;
  reco::GenParticle* matchedGen=nullptr;
  for( auto& aGen : aGens ) {
    float deltaR = reco::deltaR(aGen, *aReco );
    if ( deltaR < minDR ) { matchedGen=&aGen; minDR = deltaR; }
  }
  if ( minDR < matchingDeltaR_ ) {
    return matchedGen;
  }
  return nullptr;
}

int CATDstarAnalyzer::isFromTop( const reco::Candidate* p)
{
  if ( !p ) return 0;

  for ( size_t i=0, n=p->numberOfMothers(); i<n; ++i ) {
    const auto m = p->mother(i);
    if ( !m ) continue;
    if ( std::abs(m->pdgId()) == 6 ) return m->pdgId();
    int mm = isFromTop( m ) ;
    if ( abs(mm)>0 ) return mm;
  }

  return 0;
}

int CATDstarAnalyzer::isFromB( const reco::Candidate* p)
{
  if ( !p ) return 0;

  for ( size_t i=0, n=p->numberOfMothers(); i<n; ++i ) {
    const auto m = p->mother(i);
    if ( !m ) continue;
    if ( std::abs(m->pdgId()) == 5 ) return m->pdgId();
    int mm = isFromB( m );
    if ( abs(mm)>0 ) return mm;
  }

  return 0;
}




void CATDstarAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  //b_run = iEvent.id().run();
  //b_event = iEvent.id().event();

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


  vector<reco::GenParticle> gen_d0s;
  vector<reco::GenParticle> gen_dstars;
  edm::Handle<edm::View<reco::GenParticle> > mcHandle;
 
  if ( runOnMC ) { 
    iEvent.getByToken(mcSrc_, mcHandle);
    for( const auto& aGenParticle : *mcHandle) {
      //If genParticle is D0,
      if ( std::abs(aGenParticle.pdgId()) == 421 )       gen_d0s.push_back( aGenParticle);  
      else if ( std::abs(aGenParticle.pdgId()) ==  413 ) gen_dstars.push_back( aGenParticle);
    } 
  }

  int d0_count=-1;
  int dstar_count=-1;

  TClonesArray& br_d0 = *b_d0;
  TClonesArray& br_d0_dau1 = *b_d0_dau1;
  TClonesArray& br_d0_dau2 = *b_d0_dau2;

  TClonesArray& br_dstar = *b_dstar;
  TClonesArray& br_dstar_dau1 = *b_dstar_dau1;
  TClonesArray& br_dstar_dau2 = *b_dstar_dau2;
  TClonesArray& br_dstar_dau3 = *b_dstar_dau3;

  bool checkDuplicate= false;
  for( auto& x : *d0s) {
    d0_count++; 
    auto d0_tlv = ToTLorentzVector(x);
    bool duplicated = false;
    for ( int i=0 ; i<d0_count ; i++) {
      TLorentzVector* past_d0 = dynamic_cast<TLorentzVector*>(br_d0.At(i));
      if ( past_d0->DeltaR( d0_tlv) < 0.05 && abs(past_d0->Pt() - d0_tlv.Pt())/d0_tlv.Pt()<0.05  ) duplicated = true;
    }
    if ( checkDuplicate && duplicated) { d0_count--; continue; }

    auto lep1_d0_tlv = b_lep1+d0_tlv;
    auto lep2_d0_tlv = b_lep2+d0_tlv;

    if ( b_lep1.DeltaR( d0_tlv)< b_lep2.DeltaR( d0_tlv)) {
      b_d0_lepSV_dRM1.push_back( lep1_d0_tlv.M());
      //b_d0_lepSV_dRM2.push_back( lep2_d0_tlv.M());
    }
    else { 
      b_d0_lepSV_dRM1.push_back( lep2_d0_tlv.M());
      //b_d0_lepSV_dRM2.push_back( lep1_d0_tlv.M());
    }
  
    if ( lep1_d0_tlv.M()< lep2_d0_tlv.M() ) { 
      b_d0_lepSV_lowM1.push_back( lep1_d0_tlv.M());
      //b_d0_lepSV_lowM2.push_back( lep2_d0_tlv.M());
    }
    else { 
      b_d0_lepSV_lowM1.push_back( lep2_d0_tlv.M());
      //b_d0_lepSV_lowM2.push_back( lep1_d0_tlv.M());
    }
    b_d0_lepSV_correctM.push_back(-9);

    new( br_d0[d0_count]) TLorentzVector(d0_tlv);
    new( br_d0_dau1[d0_count]) TLorentzVector(ToTLorentzVector(*(x.daughter(0)))); 
    new( br_d0_dau2[d0_count]) TLorentzVector(ToTLorentzVector(*(x.daughter(1))));
    b_d0_dca.push_back( x.dca());

    double d0_vProb = x.vProb();
    if ( abs( d0_vProb ) > 1e-5 && abs( d0_vProb)<1e5 ) {
      b_d0_fit.push_back(true);
      b_d0_L3D.push_back( x.l3D());
      b_d0_LXY.push_back( x.lxy());
      b_d0_vProb.push_back(d0_vProb);
    }
    else {
      b_d0_fit.push_back(false);
      b_d0_L3D.push_back( -9 );
      b_d0_LXY.push_back( -9 );
      b_d0_vProb.push_back(-9);
    }        

    reco::GenParticle* genMatched = mcMatching( gen_d0s, &x ); 
    if ( genMatched != nullptr) {
      b_d0_true.push_back( true );
      b_d0_dRTrue.push_back( reco::deltaR( *genMatched, x ) ); 
      b_d0_relPtTrue.push_back( (genMatched->pt()- d0_tlv.Pt())/genMatched->pt());
      b_d0_isFromB.push_back( isFromB(genMatched));
      b_d0_isFromTop.push_back( isFromTop(genMatched));
    }
    else {
      b_d0_true.push_back( false );
      b_d0_dRTrue.push_back( -9);
      b_d0_relPtTrue.push_back(-9);
      b_d0_isFromB.push_back( -9 );
      b_d0_isFromTop.push_back( -9 );
    }
    
    b_d0_trackQuality.push_back( x.trackQuality1() + x.trackQuality2());
    b_d0_dau1_q.push_back  ( x.daughter(0)->charge());
    b_d0_dau2_q.push_back  ( x.daughter(1)->charge());
  }
  for( auto& x : *dstars) {
    dstar_count++;
    auto dstar_tlv = ToTLorentzVector(x);
    bool duplicated = false;
    for ( int i=0 ; i<dstar_count ; i++) {
      TLorentzVector* past_dstar = dynamic_cast<TLorentzVector*>(br_dstar.At(i));
      if ( past_dstar->DeltaR( dstar_tlv) < 0.05 && abs(past_dstar->Pt() - dstar_tlv.Pt())/dstar_tlv.Pt()<0.05  ) duplicated = true;
    }
    if ( checkDuplicate && duplicated) { dstar_count--; continue; }

    auto lep1_dstar_tlv = b_lep1+dstar_tlv;
    auto lep2_dstar_tlv = b_lep2+dstar_tlv;

    TLorentzVector dstar_dau1 = ToTLorentzVector(*(x.daughter(0)));
    TLorentzVector dstar_dau2 = ToTLorentzVector(*(x.daughter(1)));
    TLorentzVector dstar_dau3 = ToTLorentzVector(*(x.daughter(2)));
    new( br_dstar[dstar_count]) TLorentzVector( dstar_tlv );
    new( br_dstar_dau1[dstar_count]) TLorentzVector( dstar_dau1 ); 
    new( br_dstar_dau2[dstar_count]) TLorentzVector( dstar_dau2 );
    new( br_dstar_dau3[dstar_count]) TLorentzVector( dstar_dau3 );

    b_dstar_dca.push_back( x.dca());
    b_dstar_dca2.push_back( x.dca(1));
    b_dstar_dca3.push_back( x.dca(2));

    double dstar_vProb = x.vProb();
    if ( abs( dstar_vProb) > 1e-5 && abs( dstar_vProb)<1e5 ) {
      b_dstar_fit.push_back(true);
      b_dstar_L3D.push_back( x.l3D());
      b_dstar_LXY.push_back( x.lxy());
      b_dstar_vProb.push_back( dstar_vProb);
    }
    else {
      b_dstar_fit.push_back(false);
      b_dstar_L3D.push_back( -9 );
      b_dstar_LXY.push_back( -9 );
      b_dstar_vProb.push_back( -9);
    }
    reco::GenParticle* genMatched = mcMatching( gen_dstars, &x ); 
    if ( genMatched != nullptr) {
      b_dstar_true.push_back( true );
      b_dstar_dRTrue.push_back( reco::deltaR( *genMatched, x ));
      b_dstar_relPtTrue.push_back( (genMatched->pt()- dstar_tlv.Pt())/genMatched->pt());
      b_dstar_isFromB.push_back( isFromB(genMatched) );
      b_dstar_isFromTop.push_back( isFromTop(genMatched));
      if ( abs(isFromTop(genMatched))==6 && isFromTop(genMatched)*b_lep1_pid<0 ) {
        b_dstar_lepSV_correctM.push_back( lep1_dstar_tlv.M());
        b_dstar_lepSV_wrongM.push_back( lep2_dstar_tlv.M());
      }
      else if( abs(isFromTop(genMatched))==6 && isFromTop(genMatched)*b_lep1_pid>0){
        b_dstar_lepSV_correctM.push_back( lep2_dstar_tlv.M());
        b_dstar_lepSV_wrongM.push_back( lep1_dstar_tlv.M());
      }
      else {
        b_dstar_lepSV_correctM.push_back( -9 );
        b_dstar_lepSV_wrongM.push_back( -9 );
      }
    }
    else {
      b_dstar_true.push_back( false );
      b_dstar_dRTrue.push_back( -9);
      b_dstar_relPtTrue.push_back(-9);
      b_dstar_isFromB.push_back( -9 );
      b_dstar_isFromTop.push_back( -9 );
      b_dstar_lepSV_correctM.push_back( -9 );
      b_dstar_lepSV_wrongM.push_back( -9 );
    }


    b_dstar_trackQuality.push_back( x.trackQuality1() + x.trackQuality2());
    b_dstar_dau1_q.push_back  ( x.daughter(0)->charge());

    b_dstar_dau2_q.push_back  ( x.daughter(1)->charge());

    b_dstar_dau3_q.push_back  ( x.daughter(2)->charge());
    
    float diffMass = (dstar_dau1+dstar_dau2+dstar_dau3).M() - (dstar_dau1+dstar_dau2).M();
    
    b_dstar_diffMass.push_back( diffMass);

    // +pid = minus charge, +pid * + charge  == opposite sign
    if ( b_lep1_pid*x.daughter(2)->charge() > 0 ) {
      b_dstar_opCharge_M.push_back( lep1_dstar_tlv.M());
    }
    else {
      b_dstar_opCharge_M.push_back( lep2_dstar_tlv.M());
    }
    if ( b_lep1.DeltaR( dstar_tlv)< b_lep2.DeltaR( dstar_tlv)) {
      b_dstar_lepSV_dRM1.push_back( lep1_dstar_tlv.M());
      b_dstar_lepSV_dRM2.push_back( lep2_dstar_tlv.M());
    }
    else { 
      b_dstar_lepSV_dRM1.push_back( lep2_dstar_tlv.M());
      b_dstar_lepSV_dRM2.push_back( lep1_dstar_tlv.M());
    }
  
    if ( lep1_dstar_tlv.M()< lep2_dstar_tlv.M() ) { 
      b_dstar_lepSV_lowM1.push_back( lep1_dstar_tlv.M());
      b_dstar_lepSV_lowM2.push_back( lep2_dstar_tlv.M());
    }
    else { 
      b_dstar_lepSV_lowM1.push_back( lep2_dstar_tlv.M());
      b_dstar_lepSV_lowM2.push_back( lep1_dstar_tlv.M());
    }
      

  }
  //std::cout<<"d0 count : "<<d0_count<<"  dstar count : "<<dstar_count<<std::endl;

}

void CATDstarAnalyzer::resetBrCustom()
{
  b_d0->Clear();    b_d0_dau1->Clear();    b_d0_dau2->Clear();
  b_dstar->Clear(); b_dstar_dau1->Clear(); b_dstar_dau2->Clear(); b_dstar_dau3->Clear();


  b_d0_true.clear() ; 
  b_d0_LXY.clear(); b_d0_L3D.clear(); b_d0_vProb.clear(); 
  b_d0_fit.clear(); b_d0_dRTrue.clear(); b_d0_relPtTrue.clear(); 
  b_d0_dca.clear();
  b_d0_trackQuality.clear(); b_d0_isFromB.clear(); b_d0_isFromTop.clear();

  b_d0_lepSV_lowM1.clear();
  b_d0_lepSV_dRM1.clear();
  b_d0_lepSV_correctM.clear();

  b_dstar_true.clear(); 
  b_dstar_LXY.clear(); b_dstar_L3D.clear(); b_dstar_vProb.clear(); 
  b_dstar_fit.clear(); b_dstar_dRTrue.clear(); b_dstar_relPtTrue.clear(); 
  b_dstar_dca.clear(); b_dstar_dca2.clear(); b_dstar_dca3.clear();
  b_dstar_diffMass.clear();
  b_dstar_trackQuality.clear(); b_dstar_isFromB.clear(); b_dstar_isFromTop.clear();

  b_d0_dau1_q.clear();
  b_d0_dau2_q.clear();

  b_dstar_dau1_q.clear();
  b_dstar_dau2_q.clear();
  b_dstar_dau3_q.clear();

  b_dstar_lepSV_lowM1.clear();
  b_dstar_lepSV_lowM2.clear();
  b_dstar_lepSV_dRM1.clear();
  b_dstar_lepSV_dRM2.clear();

  b_dstar_lepSV_correctM.clear();
  b_dstar_lepSV_wrongM.clear();

  b_dstar_opCharge_M.clear();
}

//define this as a plug-in
DEFINE_FWK_MODULE(CATDstarAnalyzer);
