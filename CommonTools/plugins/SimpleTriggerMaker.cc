#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "DataFormats/Common/interface/Association.h"
#include "DataFormats/Common/interface/RefToPtr.h"

#include "DataFormats/HepMCCandidate/interface/GenParticleFwd.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/JetReco/interface/GenJetCollection.h"

#include <boost/regex.hpp>

#include <memory>
#include <vector>
#include <string>

using namespace edm;
using namespace std;

namespace cat {

  class SimpleTriggerMaker : public edm::EDProducer {
  public:
    explicit SimpleTriggerMaker(const edm::ParameterSet & iConfig);
    virtual ~SimpleTriggerMaker() { }

    virtual void produce(edm::Event & iEvent, const edm::EventSetup & iSetup);
  
  private:
    typedef std::vector< std::string > strings;
    typedef std::pair<std::string, int> stringintpair;
    typedef std::vector< stringintpair > stringintpairs;

    edm::EDGetTokenT<stringintpairs > InputTriggerLabel_;

    std::vector<std::string>  *triggerNames_;
    std::vector<std::string>  *triggerNames2_;
    //std::vector<int>  *triggerIsPathWithSacle_;


  };

} // namespace

cat::SimpleTriggerMaker::SimpleTriggerMaker(const edm::ParameterSet & iConfig) :
  InputTriggerLabel_(consumes< stringintpairs >(iConfig.getParameter<edm::InputTag>("InputTriggerLabel")))
{
  triggerNames_           = new std::vector<std::string>;
  triggerNames2_          = new std::vector<std::string>;
  //triggerIsPathWithSacle_ = new std::vector<int>;

  const boost::regex matchVersion("_v[0-9\\*]+$"); // regexp from HLTrigger/HLTCore/HLTConfigProvider

  std::cout << "List of Triggers to Save" << std::endl;
  for ( auto& hltPath : iConfig.getParameter<strings>("hltPathNames") ){
    hltPath = boost::regex_replace(hltPath, matchVersion, "");
    std::string hltSavedAs = hltPath;
    hltSavedAs.erase(std::remove(hltSavedAs.begin(), hltSavedAs.end(), '_'), hltSavedAs.end());
    std::cout << " " << hltPath << std::endl;
    produces<int >( hltSavedAs );
    triggerNames2_->push_back(hltSavedAs);
  }

  produces< std::vector<std::string> >("TriggerNames");
  produces< std::vector<int>         >("TriggerIsPathWithScale");
}
/*void CATTriggerProducer::beginRun(const edm::Run& run, const edm::EventSetup& eventSetup)
{
  
}*/
void cat::SimpleTriggerMaker::produce(edm::Event & iEvent, const edm::EventSetup & iSetup) 
{

  Handle<stringintpairs> Triggers;
  iEvent.getByToken(InputTriggerLabel_, Triggers);

  unsigned int n = Triggers->size();

  for ( unsigned int ip=0; ip<n; ++ip ) 
  { 
    const stringintpair& p = (*Triggers)[ip];
    int findIndex=-1;
    std::string name_   = p.first;
    //int         isPath_ = p.second;
    for(unsigned int it=0;it<triggerNames_->size();it++)
    {
      std::string triggerName_ = triggerNames_->at(it);
      std::size_t found = triggerName_.find(name_);
      if (found!=std::string::npos && name_.length()==triggerName_.length() )
      {
         findIndex=ip;
      }
    } 
    if (findIndex==-1)
    {
      triggerNames_->push_back(p.first);
    }
  }

  std::vector<std::string> *triggerNames = new std::vector<std::string>();
  std::vector<int>  *triggerIsPathWithSacle = new std::vector<int>();
  const boost::regex matchVersion("_v[0-9\\*]+$"); // regexp from HLTrigger/HLTCore/HLTConfigProvider
  for(unsigned int it=0;it<triggerNames_->size();it++)
  {
    std::string triggerName_ = triggerNames_->at(it);
    for ( unsigned int ip=0; ip<n; ++ip ) 
    {
      const stringintpair& p = (*Triggers)[ip];
      std::string name_   = p.first;
      int         isPath_ = p.second;
      std::size_t found = triggerName_.find(name_);
      if (found!=std::string::npos && name_.length()==triggerName_.length() )
      {
         triggerNames->push_back(name_);
         triggerIsPathWithSacle->push_back(isPath_);
         //std::cout << name_ <<":"<< isPath_ << std::endl;
         name_ = boost::regex_replace(name_, matchVersion, "");
         std::string hltSavedAs = name_;
         hltSavedAs.erase(std::remove(hltSavedAs.begin(), hltSavedAs.end(), '_'), hltSavedAs.end());

        for(unsigned int it2=0;it2<triggerNames2_->size();it2++)
        {
          std::size_t found2 = triggerNames2_->at(it2).find(hltSavedAs);
          if (found2!=std::string::npos && hltSavedAs.length()==triggerNames2_->at(it2).length() )
              iEvent.put(std::auto_ptr<int>(new int (isPath_)), hltSavedAs);    
              //std::cout << hltSavedAs <<":"<< isPath_ << std::endl;
        }
      }
    }
  }

  iEvent.put( std::auto_ptr<std::vector<std::string>> (triggerNames)          ,"TriggerNames" );
  iEvent.put( std::auto_ptr<std::vector<int>>         (triggerIsPathWithSacle),"TriggerIsPathWithScale" );

///////////////
}

#include "FWCore/Framework/interface/MakerMacros.h"
using namespace cat;
DEFINE_FWK_MODULE(SimpleTriggerMaker);
