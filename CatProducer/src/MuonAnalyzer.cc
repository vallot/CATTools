#include "../interface/MuonAnalyzer.h"

using namespace std;
using namespace cat;
using namespace reco;
using namespace edm;
using namespace isodeposit;

MuonAnalyzer::MuonAnalyzer (const edm::ParameterSet & producersNames):
verbosity_(0),
useMC_(false)
{
	muonProducer_ = producersNames.getParameter < edm::InputTag > ("muonProducer");
  primaryVertexProducer_ = producersNames.getParameter<edm::InputTag>("primaryVertexProducer");
}

MuonAnalyzer::MuonAnalyzer (const edm::ParameterSet & producersNames, const edm::ParameterSet & myConfig, int verbosity):
verbosity_ (verbosity)
{
	muonProducer_ = producersNames.getParameter < edm::InputTag > ("muonProducer");
  primaryVertexProducer_ = producersNames.getParameter<edm::InputTag>("primaryVertexProducer");
	useMC_ = myConfig.getUntrackedParameter < bool > ("doMuonMC");
}
MuonAnalyzer::MuonAnalyzer (const edm::ParameterSet & producersNames, int iter, const edm::ParameterSet & myConfig, int verbosity):
verbosity_ (verbosity)
{
	vMuonProducer = producersNames.getUntrackedParameter<std::vector<std::string> >("vmuonProducer");
	muonProducer_ =	edm::InputTag(vMuonProducer[iter]);
  primaryVertexProducer_ = producersNames.getParameter<edm::InputTag>("primaryVertexProducer");
	useMC_ = myConfig.getUntrackedParameter < bool > ("doMuonMC");
}

MuonAnalyzer::~MuonAnalyzer ()
{
}

void
MuonAnalyzer::Process (const edm::Event & iEvent, TClonesArray * rootMuons)
{
	unsigned int nMuons = 0;

	edm::Handle < std::vector < pat::Muon > >patMuons;
	iEvent.getByLabel (muonProducer_, patMuons);
	nMuons = patMuons->size ();

  edm::Handle< reco::VertexCollection > pvHandle;
  iEvent.getByLabel(primaryVertexProducer_, pvHandle);
  
//	edm::Handle<reco::BeamSpot> beamSpotHandle;
//	iEvent.getByLabel("offlineBeamSpot", beamSpotHandle);
//	const reco::TrackBase::Point & beamSpot = reco::TrackBase::Point(beamSpotHandle->x0(), beamSpotHandle->y0(), beamSpotHandle->z0());

	if (verbosity_ > 1)
		std::cout << "   Number of muons = " << nMuons << "   Label: " << muonProducer_.label () << "   Instance: " << muonProducer_.instance () << std::endl;

    for (unsigned int j = 0; j < nMuons; j++)
    {
      const pat::Muon*  patMuon = &((*patMuons)[j]);//dynamic_cast < const pat::Muon * >(&*muon);
      const reco::Muon* muon = (const reco::Muon *) patMuon;//(&((*patMuons)[j]));

      cat::CatMuon localMuon (muon->px (), muon->py (), muon->pz (), muon->energy (), muon->vx (), muon->vy (), muon->vz (), muon->pdgId (), muon->charge ());

      //set detector-based isolation
      localMuon.setIsoR03_trackIso( muon->isolationR03().sumPt);
      localMuon.setIsoR03_ecalIso( muon->isolationR03().emEt);
      localMuon.setIsoR03_hcalIso( muon->isolationR03().hadEt);

      localMuon.setAlgo (muon->type ());
      localMuon.setID (int ( muon::isGoodMuon ( *muon, muon::AllGlobalMuons)), int ( muon::isGoodMuon ( *muon, muon::AllTrackerMuons)), int ( muon::isGoodMuon ( *muon, muon::AllStandAloneMuons)), int ( muon::isGoodMuon ( *muon, muon::TrackerMuonArbitrated)), int ( muon::isGoodMuon ( *muon, muon::AllArbitrated)), int ( muon::isGoodMuon ( *muon, muon::GlobalMuonPromptTight)), int ( muon::isGoodMuon (*muon, muon::TMLastStationLoose)), int ( muon::isGoodMuon ( *muon, muon::TMLastStationTight)),int ( muon::isGoodMuon ( *muon, muon::TMLastStationAngTight)) , int ( muon::isGoodMuon ( *muon, muon::TMOneStationLoose)), int ( muon::isGoodMuon ( *muon, muon::TMOneStationTight)), int ( muon::isGoodMuon ( *muon, muon::TMLastStationOptimizedLowPtLoose)), int ( muon::isGoodMuon ( *muon, muon::TMLastStationOptimizedLowPtTight)), int ( muon::isGoodMuon ( *muon, muon::TM2DCompatibilityLoose)), int ( muon::isGoodMuon ( *muon, muon::TM2DCompatibilityTight)));

      if(muon->innerTrack().isNonnull() && muon->innerTrack().isAvailable())
      {
        if(pvHandle.isValid() && pvHandle->size() != 0)
        {
          reco::VertexRef vtx(pvHandle,0);
          localMuon.setD0 ( muon->innerTrack()->dxy(vtx->position()) );
          localMuon.setD0Error ( sqrt(pow(muon->innerTrack()->dxyError(),2)+pow(vtx->xError(),2)+ pow(vtx->yError(),2)) );
          localMuon.setDz ( muon->innerTrack()->dz(vtx->position()) );
          localMuon.setDzError ( sqrt(pow(muon->innerTrack()->dzError(),2)+pow(vtx->zError(),2)) );
        }
        else{
          localMuon.setD0 ( muon->innerTrack()->dxy() );
          localMuon.setD0Error ( muon->innerTrack()->dxyError() );
          localMuon.setDz ( muon->innerTrack()->dz() );
          localMuon.setDzError ( muon->innerTrack()->dzError() );
          
        }
        localMuon.setNofValidHits ( muon->innerTrack()->numberOfValidHits() );
        localMuon.setNofValidPixelHits( muon->innerTrack()->hitPattern().numberOfValidPixelHits() );
        localMuon.setNofTrackerLayersWithMeasurement( muon->innerTrack()->hitPattern().trackerLayersWithMeasurement() );
      }

      if(muon->isGlobalMuon ())
      {
        localMuon.setNofValidMuHits(muon->globalTrack()->hitPattern().numberOfValidMuonHits());
        localMuon.setNofMatchedStations(muon->numberOfMatchedStations());
      }

      // Some specific methods to pat::Muon
      localMuon.setIsPFMuon(patMuon->isPFMuon());

      /*if(patMuon->innerTrack().isNonnull() && patMuon->innerTrack().isAvailable())
       {
       //cout << "Step: " << temp << endl; temp++; // 5

       localMuon.setdB(patMuon->dB());
       localMuon.setdBError(patMuon->edB());
       }*/

      try {
        localMuon.setChi2 (patMuon->normChi2 ());
      }
      catch (cms::Exception &lce) {
        cerr << "MuonAnalyzer:: WARNING, unable to access muon normChi2 value!!!! (label: " << muonProducer_.label () << ")" << endl;
        localMuon.setChi2 (+99999.);
      }

      if (patMuon->ecalIsoDeposit ()) localMuon.setVetoEm  (patMuon->ecalIsoDeposit ()->candEnergy ());
      if (patMuon->hcalIsoDeposit ()) localMuon.setVetoHad (patMuon->hcalIsoDeposit ()->candEnergy ());

      if (useMC_)
      {
        // MC truth associator index
        if ((patMuon->genParticleRef ()).isNonnull ())
        {
          localMuon.setGenParticleIndex ((patMuon->genParticleRef ()).index ());
        }
        else
        {
          localMuon.setGenParticleIndex (-1);
        }
      }
      
      //set pf-isolation
      //default cone size from PAT setup : currently r=0.4
      if(patMuon->chargedHadronIso() != -1)   localMuon.setIsoR04_ChargedHadronIso(patMuon->chargedHadronIso());
      if(patMuon->puChargedHadronIso() != -1) localMuon.setIsoR04_PuChargedHadronIso( patMuon->puChargedHadronIso() );
      if(patMuon->photonIso() != -1)          localMuon.setIsoR04_PhotonIso(patMuon->photonIso());
      if(patMuon->neutralHadronIso() != -1)   localMuon.setIsoR04_NeutralHadronIso(patMuon->neutralHadronIso());

       //this is for cone size r=0.3
       //following Muon POG isolation definition
       AbsVetos veto_ch;
       AbsVetos veto_nh;
       AbsVetos veto_ph;
       Direction Dir = Direction(patMuon->eta(), patMuon->phi());
       veto_ch.push_back(new reco::isodeposit::ConeVeto( Dir, 0.0001 ));
       veto_nh.push_back(new reco::isodeposit::ConeVeto( Dir, 0.01 ));
       veto_nh.push_back(new reco::isodeposit::ThresholdVeto( 0.5 ));
       veto_ph.push_back(new reco::isodeposit::ConeVeto( Dir, 0.01 ));
       veto_ph.push_back(new reco::isodeposit::ThresholdVeto( 0.5 ));

       const double chIso03 = patMuon->isoDeposit(pat::PfChargedHadronIso)->depositAndCountWithin(0.3, veto_ch).first;
       const double nhIso03 = patMuon->isoDeposit(pat::PfNeutralHadronIso)->depositAndCountWithin(0.3, veto_nh).first;
       const double phIso03 = patMuon->isoDeposit(pat::PfGammaIso)->depositAndCountWithin(0.3, veto_ph).first;
       const double puChIso03 = patMuon->isoDeposit(pat::PfPUChargedHadronIso)->depositAndCountWithin(0.3, veto_ch).first;

       localMuon.setIsoR03_ChargedHadronIso( chIso03 );
       localMuon.setIsoR03_PuChargedHadronIso( puChIso03 );
       localMuon.setIsoR03_PhotonIso( phIso03 );
       localMuon.setIsoR03_NeutralHadronIso( nhIso03 );

      new ((*rootMuons)[j]) cat::CatMuon (localMuon);
      if (verbosity_ > 2)
        cout << "   [" << setw (3) << j << "] " << localMuon << endl;
	}
}
