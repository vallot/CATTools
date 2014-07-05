
#include "../interface/RecoMETAnalyzer.h"

using namespace std;
using namespace cat;
using namespace reco;
using namespace edm;

RecoMETAnalyzer::RecoMETAnalyzer(const edm::ParameterSet& producersNames):verbosity_(0), useMC_(false)
{
}

RecoMETAnalyzer::RecoMETAnalyzer(const edm::ParameterSet& producersNames, const edm::ParameterSet& myConfig, int verbosity):verbosity_(verbosity)
{
	useMC_ = myConfig.getUntrackedParameter<bool>("doMETMC");
}

RecoMETAnalyzer::~RecoMETAnalyzer()
{
}

cat::CatMET RecoMETAnalyzer::Process(const reco::Candidate* met)
{

  
  cat::CatMET localMET(
		    met->px()
		    ,met->py()
		    ,met->pz()
		    ,met->energy()
		    ,met->vx()
		    ,met->vy()
		    ,met->vz()
		    ); 
  
      const reco::PFMET *patMET = dynamic_cast<const reco::PFMET*>(&*met);
      
      localMET.setSumEt(patMET->sumEt());
      
     

  
  return localMET;

}
