#include "../interface/JetAnalyzer.h"

using namespace std;
using namespace cat;
using namespace reco;
using namespace edm;

JetAnalyzer::JetAnalyzer():verbosity_(0),useMC_(false)
{
}

JetAnalyzer::JetAnalyzer(int verbosity):verbosity_(verbosity),useMC_(false)
{
}

JetAnalyzer::JetAnalyzer(const edm::ParameterSet& myConfig, int verbosity):verbosity_(verbosity)
{
	useMC_ = myConfig.getUntrackedParameter<bool>("doJetMC");
	isData_ = myConfig.getUntrackedParameter<bool>("isData");
}

JetAnalyzer::~JetAnalyzer()
{
}

bool Rsortrule (std::pair <double,double> p1, std::pair <double,double> p2 )
{
	return p1.second<p2.second; 
}

cat::Jet JetAnalyzer::Process(const reco::Jet* jet, const edm::EventSetup& iSetup)
{
	cat::Jet localJet(
		jet->px()
		,jet->py()
		,jet->pz()
		,jet->energy()
		,jet->vx()
		,jet->vy()
		,jet->vz()
		,jet->pdgId()
		,jet->charge()
	); 

	// Some specific methods to pat::Jet
	const pat::Jet *patJet = dynamic_cast<const pat::Jet*>(&*jet);

	localJet.setNConstituents(jet->nConstituents());
	localJet.setJetArea(jet->jetArea()); // needs to be calculated with fastjet, not yet done in the standard cfg's
	localJet.setMaxDistance(jet->maxDistance());
	
	// Variables from pat::Jet (Basic)
	localJet.setBtag_jetBProbabilityBJetTags(patJet->bDiscriminator("jetBProbabilityBJetTags"));
	localJet.setBtag_jetProbabilityBJetTags(patJet->bDiscriminator("jetProbabilityBJetTags"));
	localJet.setBtag_trackCountingHighPurBJetTags(patJet->bDiscriminator("trackCountingHighPurBJetTags"));
	localJet.setBtag_trackCountingHighEffBJetTags(patJet->bDiscriminator("trackCountingHighEffBJetTags"));
	localJet.setBtag_simpleSecondaryVertexHighEffBJetTags(patJet->bDiscriminator("simpleSecondaryVertexHighEffBJetTags"));
	localJet.setBtag_simpleSecondaryVertexHighPurBJetTags(patJet->bDiscriminator("simpleSecondaryVertexHighPurBJetTags"));
	localJet.setBtag_combinedSecondaryVertexBJetTags(patJet->bDiscriminator("combinedSecondaryVertexBJetTags"));
	localJet.setBtag_combinedSecondaryVertexRetrainedBJetTags(patJet->bDiscriminator("combinedSecondaryVertexRetrainedBJetTags"));
	localJet.setBtag_combinedSecondaryVertexMVABJetTags(patJet->bDiscriminator("combinedSecondaryVertexMVABJetTags"));
	localJet.setBtag_softMuonBJetTags(patJet->bDiscriminator("softMuonBJetTags"));
	localJet.setBtag_softMuonByPtBJetTags(patJet->bDiscriminator("softMuonByPtBJetTags"));
	localJet.setBtag_softMuonByIP3dBJetTags(patJet->bDiscriminator("softMuonByIP3dBJetTags"));
	localJet.setBtag_softElectronByPtBJetTags(patJet->bDiscriminator("softElectronByPtBJetTags"));
	localJet.setBtag_softElectronByIP3dBJetTags(patJet->bDiscriminator("softElectronByIP3dBJetTags"));
	localJet.setBtag_combinedCSVJPBJetTags(patJet->bDiscriminator("combinedCSVJPBJetTags"));
	localJet.setBtag_combinedCSVJPSLBJetTags(patJet->bDiscriminator("combinedCSVJPSLBJetTags"));
	localJet.setBtag_combinedCSVSLBJetTags(patJet->bDiscriminator("combinedCSVSLBJetTags"));
	localJet.setBtag_softPFElectronRetrainedBJetsTags(patJet->bDiscriminator("softPFElectronBJetTags"));
	localJet.setBtag_softPFMuonRetrainedBJetsTags(patJet->bDiscriminator("softPFMuonBJetTags"));

	//cout << "CSV old, new: " << patJet->bDiscriminator("combinedSecondaryVertexBJetTags") << ", " << patJet->bDiscriminator("combinedSecondaryVertexRetrainedBJetTags") << endl;

// comment out the whole b-tag scale factor setup from DB below for now	
/*
  // Save b-tag scalefactor information for each tagger
  ///////////////////////////////
  ///   Begin DB setup
 ///////////////////////////////

  //// This is needed for the DB
  std::map<std::string,PerformanceResult::ResultType> measureMap;
  //measureMap["BTAGBEFF"]=PerformanceResult::BTAGBEFF; //absolute efficiencies
  //measureMap["BTAGBERR"]=PerformanceResult::BTAGBERR;
  //measureMap["BTAGLEFF"]=PerformanceResult::BTAGLEFF;
  //measureMap["BTAGLERR"]=PerformanceResult::BTAGLERR;
  measureMap["BTAGBEFFCORR"]=PerformanceResult::BTAGBEFFCORR; //scalefactors have suffix "CORR"
  measureMap["BTAGBERRCORR"]=PerformanceResult::BTAGBERRCORR;
  measureMap["BTAGLEFFCORR"]=PerformanceResult::BTAGLEFFCORR;
  measureMap["BTAGLERRCORR"]=PerformanceResult::BTAGLERRCORR;

  std::vector<std::string> measureName;
  std::vector<std::string> measureType;
  
  // Define which Btag and Mistag algorithm you want to use. These are not user defined and need to be exact
  measureName.push_back("MISTAGTCHEL");  measureName.push_back("BTAGTCHEL");  measureName.push_back("MISTAGTCHEL");  measureName.push_back("BTAGTCHEL");
  measureName.push_back("MISTAGTCHEM");  measureName.push_back("BTAGTCHEM");  measureName.push_back("MISTAGTCHEM");  measureName.push_back("BTAGTCHEM");
  measureName.push_back("MISTAGTCHPM");  measureName.push_back("BTAGTCHPM");  measureName.push_back("MISTAGTCHPM");  measureName.push_back("BTAGTCHPM");
  measureName.push_back("MISTAGTCHPT");  measureName.push_back("BTAGTCHPT");  measureName.push_back("MISTAGTCHPT");  measureName.push_back("BTAGTCHPT");
  measureName.push_back("MISTAGJPL");  measureName.push_back("BTAGJPL");  measureName.push_back("MISTAGJPL");  measureName.push_back("BTAGJPL");
  measureName.push_back("MISTAGJPM");  measureName.push_back("BTAGJPM");  measureName.push_back("MISTAGJPM");  measureName.push_back("BTAGJPM");
  measureName.push_back("MISTAGJPT");  measureName.push_back("BTAGJPT");  measureName.push_back("MISTAGJPT");  measureName.push_back("BTAGJPT");
  measureName.push_back("MISTAGCSVL");  measureName.push_back("BTAGCSVL");  measureName.push_back("MISTAGCSVL");  measureName.push_back("BTAGCSVL");
  measureName.push_back("MISTAGCSVM");  measureName.push_back("BTAGCSVM");  measureName.push_back("MISTAGCSVM");  measureName.push_back("BTAGCSVM");
  measureName.push_back("MISTAGCSVT");  measureName.push_back("BTAGCSVT");  measureName.push_back("MISTAGCSVT");  measureName.push_back("BTAGCSVT");
  measureName.push_back("MISTAGSSVHEM");  measureName.push_back("BTAGSSVHEM");  measureName.push_back("MISTAGSSVHEM");  measureName.push_back("BTAGSSVHEM");
  measureName.push_back("MISTAGSSVHPT");  measureName.push_back("BTAGSSVHPT");  measureName.push_back("MISTAGSSVHPT");  measureName.push_back("BTAGSSVHPT");


  // Tell DB you want the SF. These are not user defined and need to be exact
  // for the light eff or scalefactor, one needs the mistag algorithm in measureName at the same location
  // for the b eff or scalefactor, one needs the btag algorithm in measureName at the same location
  measureType.push_back("BTAGLEFFCORR");  measureType.push_back("BTAGBEFFCORR");  measureType.push_back("BTAGLERRCORR");	 measureType.push_back("BTAGBERRCORR");
  measureType.push_back("BTAGLEFFCORR");  measureType.push_back("BTAGBEFFCORR");  measureType.push_back("BTAGLERRCORR");	 measureType.push_back("BTAGBERRCORR");
  measureType.push_back("BTAGLEFFCORR");  measureType.push_back("BTAGBEFFCORR");  measureType.push_back("BTAGLERRCORR");	 measureType.push_back("BTAGBERRCORR");
  measureType.push_back("BTAGLEFFCORR");  measureType.push_back("BTAGBEFFCORR");  measureType.push_back("BTAGLERRCORR");	 measureType.push_back("BTAGBERRCORR");
  measureType.push_back("BTAGLEFFCORR");  measureType.push_back("BTAGBEFFCORR");  measureType.push_back("BTAGLERRCORR");	 measureType.push_back("BTAGBERRCORR");
  measureType.push_back("BTAGLEFFCORR");  measureType.push_back("BTAGBEFFCORR");  measureType.push_back("BTAGLERRCORR");	 measureType.push_back("BTAGBERRCORR");
  measureType.push_back("BTAGLEFFCORR");  measureType.push_back("BTAGBEFFCORR");  measureType.push_back("BTAGLERRCORR");	 measureType.push_back("BTAGBERRCORR");
  measureType.push_back("BTAGLEFFCORR");  measureType.push_back("BTAGBEFFCORR");  measureType.push_back("BTAGLERRCORR");	 measureType.push_back("BTAGBERRCORR");
  measureType.push_back("BTAGLEFFCORR");  measureType.push_back("BTAGBEFFCORR");  measureType.push_back("BTAGLERRCORR");	 measureType.push_back("BTAGBERRCORR");
  measureType.push_back("BTAGLEFFCORR");  measureType.push_back("BTAGBEFFCORR");  measureType.push_back("BTAGLERRCORR");	 measureType.push_back("BTAGBERRCORR");
  measureType.push_back("BTAGLEFFCORR");  measureType.push_back("BTAGBEFFCORR");  measureType.push_back("BTAGLERRCORR");	 measureType.push_back("BTAGBERRCORR");
  measureType.push_back("BTAGLEFFCORR");  measureType.push_back("BTAGBEFFCORR");  measureType.push_back("BTAGLERRCORR");	 measureType.push_back("BTAGBERRCORR");
		
  // These are user defined maps that we will use to store the SF 
  std::map<std::string,float> MISTAG_SF;
  std::map<std::string,float> BTAG_SF;
  std::map<std::string,float> MISTAG_SFerr;
  std::map<std::string,float> BTAG_SFerr;
	
  ///////////////////////////////
  ///   End DB setup
  ///////////////////////////////
  //cout<< "NEW JET" << endl;
  //cout<< "measureName size: " << measureName.size() << endl;
  //cout<< "measureType size: " << measureType.size() << endl;
  //cout << "jet ET " << patJet->et() << endl;
  //cout << "jet eta " << abs( patJet->eta() ) << endl;
  edm::ESHandle<BtagPerformance> perfH;
  for( size_t iMeasure = 0; iMeasure < measureName.size(); iMeasure++ )
  {
    //std::cout << "Testing: " << measureName[ iMeasure ] << " of type " << measureType[ iMeasure ] << std::endl;

    //Setup our measurement
    iSetup.get<BTagPerformanceRecord>().get( measureName[ iMeasure ],perfH);
    const BtagPerformance & perf = *(perfH.product());

    //Working point
    //std::cout << "Working point: " << perf.workingPoint().cut() << std::endl;
		
    //Setup the point we wish to test!
    BinningPointByMap measurePoint;
    measurePoint.reset();
    measurePoint.insert(BinningVariables::JetEt, patJet->et());
    measurePoint.insert(BinningVariables::JetAbsEta, abs(  patJet->eta() ));

    std::string suffix = "_"+measureType[iMeasure];
    //std::cout << "measureType = " << suffix << endl;
    //std::cout << measureName[ iMeasure ] + suffix << endl;
    if(measureType[iMeasure] == "BTAGLEFFCORR"){					
      //std::cout << "mistag_SF " << perf.getResult( measureMap[ measureType[ iMeasure] ], measurePoint) << endl;
      MISTAG_SF[ measureName[ iMeasure ] + suffix ] = perf.getResult( measureMap[ measureType[ iMeasure] ], measurePoint);
    }else if(measureType[iMeasure] == "BTAGBEFFCORR"){
      //std::cout << "btag_SF " << perf.getResult( measureMap[ measureType[ iMeasure] ], measurePoint) << endl;
      BTAG_SF[ measureName[ iMeasure ] + suffix ] = perf.getResult( measureMap[ measureType[ iMeasure] ], measurePoint);
    }else if(measureType[iMeasure] == "BTAGLERRCORR"){
      //std::cout << "mistag_SFerr " << perf.getResult( measureMap[ measureType[ iMeasure] ], measurePoint) << endl;
      MISTAG_SFerr[ measureName[ iMeasure ] + suffix ] = perf.getResult( measureMap[ measureType[ iMeasure] ], measurePoint);
    }else if(measureType[iMeasure] == "BTAGBERRCORR"){
      //std::cout << "btag_SFerr " << perf.getResult( measureMap[ measureType[ iMeasure] ], measurePoint) << endl;
      BTAG_SFerr[ measureName[ iMeasure ] + suffix ] = perf.getResult( measureMap[ measureType[ iMeasure] ], measurePoint);
    }
  }

  //cout << "mistag_SF map size " << MISTAG_SF.size() << endl;
  //for(std::map<std::string,float>::const_iterator it = MISTAG_SF.begin(); it != MISTAG_SF.end(); it++) {
  //	cout << "for " << it->first << " the mistag scalefactor is " << it->second << endl;
  //}
	
  localJet.setMistag_SF(MISTAG_SF);
  localJet.setBtag_SF(BTAG_SF);
  localJet.setMistag_SFerr(MISTAG_SFerr);
  localJet.setBtag_SFerr(BTAG_SFerr);
	
	
  //////////////////// end of saving b-tagging information
*/
	
	
	// jet correction factors
	std::vector< std::string > jecLevels = patJet->availableJECLevels();
	
	pat::Jet rawJet = patJet->correctedJet("Uncorrected");
	
	localJet.setJetCorrFactor(0,jecLevels[1],rawJet.jecFactor(jecLevels[1]));
	localJet.setJetCorrFactor(1,jecLevels[1]+"L2",rawJet.jecFactor("L2Relative"));
	localJet.setJetCorrFactor(2,jecLevels[1]+"L2L3",rawJet.jecFactor("L3Absolute"));
	if(jecLevels.size() > 4 && jecLevels[4] == "L2L3Residual" )
		localJet.setJetCorrFactor(3,jecLevels[1]+"L2L3L23Residual",rawJet.jecFactor("L2L3Residual"));
		
	// Matched genParticle
	if (useMC_)
	{
		// MC truth associator index
		if ((patJet->genParticleRef()).isNonnull())
			localJet.setGenParticleIndex((patJet->genParticleRef()).index());
		else
			localJet.setGenParticleIndex(-1);

		// set the parton flavour
		localJet.setPartonFlavour(patJet->partonFlavour());

		// check if jet comes from a top
		bool IsTopJet =  false;
		if(patJet->genParton())
		{
			const reco::Candidate* gen = patJet->genParton();
			while(gen->mother())
			{
				if(abs((gen->mother())->pdgId()) == 6)
				{
					IsTopJet =  true;
					break;
				}
				else
				{
					gen = (gen->mother() );
				}
			}
		}
		localJet.setIsTopJet(IsTopJet);
	}

	return localJet;
}
