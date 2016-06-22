#include "CATTools/CatAnalyzer/interface/CATDstarAnalyzer.h"

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
  //parameterInit(iConfig);

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
    if( sys == 0 ) genInfo(iEvent, iSetup);
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

  edm::Handle<edm::View<reco::GenParticle> > mcHandle;
  iEvent.getByToken(mcSrc_, mcHandle);

  vector<TLorentzVector> gen_d0s;
  vector<TLorentzVector> gen_dstars;

  for( const auto& aGenParticle : *mcHandle) {
    //If genParticle is D0,
    if ( std::abs(aGenParticle.pdgId()) == 421 )       gen_d0s.push_back( ToTLorentzVector(aGenParticle));  
    else if ( std::abs(aGenParticle.pdgId()) ==  413 ) gen_dstars.push_back( ToTLorentzVector(aGenParticle));
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

  for( auto& x : *d0s) {
    d0_count++; 

    auto d0_tlv = ToTLorentzVector(x);
    new( br_d0[d0_count]) TLorentzVector(d0_tlv);
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
    

    b_d0_dau1_q.push_back  ( x.daughter(0)->charge());
    b_d0_dau2_q.push_back  ( x.daughter(1)->charge());
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
#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(CATDstarAnalyzer);
