// -*- C++ -*-
//
// Package:    PhysicsTools/GenTtbarCategorizer
// Class:      GenTtbarCategorizer
// 
/**\class GenTtbarCategorizer GenTtbarCategorizer.cc PhysicsTools/JetMCAlgos/plugins/GenTtbarCategorizer.cc
//https://twiki.cern.ch/twiki/pub/CMSPublic/GenHFHadronMatcher/GenTtbarCategorizer.cc
//
 Description: Categorization of different tt+xx processes, returning unique ID for each process as e.g. tt+bb, tt+b, tt+2b, tt+cc, ...

 Implementation:
     
     The classification scheme returns an ID per event, and works as follows:
     
     All jets in the following need to be in the acceptance as given by the config parameters |eta|, pt.
     
     First, jets from top are identified, i.e. jets containing a b hadron from top. These are excluded from the search for additional jets.
     They are encoded in the ID as numberOfBjetsFromTop*100, i.e.
     0xx: no b jets from top in acceptance
     1xx: 1 b jet from top in acceptance
     2xx: both b jets from top in acceptance
     
     From the remaining jets, the ID is formed based on the additional b jets (IDs x5x) and c jets (IDs x4x) in the following order:
     x55: at least 2 additional b jets with two of them having >= 2 b hadrons
     x54: at least 2 additional b jets with one of them having >= 2 b hadrons, the other having =1 b hadron
     x53: at least 2 additional b jets with all having =1 b hadron
     x52: exactly 1 additional b jet having >=2 b hadrons
     x51: exactly 1 additional b jet having =1 b hadron
     x45: at least 2 additional c jets with two of them having >= 2 c hadrons
     x44: at least 2 additional c jets with one of them having >= 2 c hadrons, the other having =1 c hadron
     x43: at least 2 additional c jets with all having =1 c hadron
     x42: exactly 1 additional c jet having >=2 c hadrons
     x41: exactly 1 additional c jet having =1 c hadron
     x00: No additional b or c jet, i.e. only light flavour jets
*/
//
// Original Author:  Johannes Hauk
//         Created:  Sun, 14 Jun 2015 19:42:58 GMT
//
//


// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/stream/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "DataFormats/JetReco/interface/GenJetCollection.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"


//
// class declaration
//

class GenTtbarCategorizer : public edm::stream::EDProducer<> {
    public:
        explicit GenTtbarCategorizer(const edm::ParameterSet&);
        ~GenTtbarCategorizer() {};
        
        static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);
        
    private:
        void produce(edm::Event&, const edm::EventSetup&) override;
        
        //virtual void beginRun(edm::Run const&, edm::EventSetup const&) override;
        //virtual void endRun(edm::Run const&, edm::EventSetup const&) override;
        //virtual void beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) override;
        //virtual void endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) override;
        
        // ----------member data ---------------------------
        
        // Jet configuration
        const double genJetPtMin_;
        const double genJetAbsEtaMax_;
        
        // Input tags
        const edm::EDGetTokenT<reco::GenJetCollection> genJetsToken_;
        
        const edm::EDGetTokenT<std::vector<int> > genBHadJetIndexToken_;
        const edm::EDGetTokenT<std::vector<int> > genBHadFlavourToken_;
        const edm::EDGetTokenT<std::vector<int> > genBHadFromTopWeakDecayToken_;
        const edm::EDGetTokenT<std::vector<reco::GenParticle> > genBHadPlusMothersToken_;
        const edm::EDGetTokenT<std::vector<std::vector<int> > > genBHadPlusMothersIndicesToken_;
        const edm::EDGetTokenT<std::vector<int> > genBHadIndexToken_;
        const edm::EDGetTokenT<std::vector<int> > genBHadLeptonHadronIndexToken_;
        const edm::EDGetTokenT<std::vector<int> > genBHadLeptonViaTauToken_;
        
        const edm::EDGetTokenT<std::vector<int> > genCHadJetIndexToken_;
        const edm::EDGetTokenT<std::vector<int> > genCHadFlavourToken_;
        const edm::EDGetTokenT<std::vector<int> > genCHadFromTopWeakDecayToken_;
        const edm::EDGetTokenT<std::vector<int> > genCHadBHadronIdToken_;
        
        
};

//
// constructors and destructor
//
GenTtbarCategorizer::GenTtbarCategorizer(const edm::ParameterSet& iConfig):
genJetPtMin_(iConfig.getParameter<double>("genJetPtMin")),
genJetAbsEtaMax_(iConfig.getParameter<double>("genJetAbsEtaMax")),
genJetsToken_(consumes<reco::GenJetCollection>(iConfig.getParameter<edm::InputTag>("genJets"))),
genBHadJetIndexToken_(consumes<std::vector<int> >(iConfig.getParameter<edm::InputTag>("genBHadJetIndex"))),
genBHadFlavourToken_(consumes<std::vector<int> >(iConfig.getParameter<edm::InputTag>("genBHadFlavour"))),
genBHadFromTopWeakDecayToken_(consumes<std::vector<int> >(iConfig.getParameter<edm::InputTag>("genBHadFromTopWeakDecay"))),
genBHadPlusMothersToken_(consumes<std::vector<reco::GenParticle> >(iConfig.getParameter<edm::InputTag>("genBHadPlusMothers"))),
genBHadPlusMothersIndicesToken_(consumes<std::vector<std::vector<int> > >(iConfig.getParameter<edm::InputTag>("genBHadPlusMothersIndices"))),
genBHadIndexToken_(consumes<std::vector<int> >(iConfig.getParameter<edm::InputTag>("genBHadIndex"))),
genBHadLeptonHadronIndexToken_(consumes<std::vector<int> >(iConfig.getParameter<edm::InputTag>("genBHadLeptonHadronIndex"))),
genBHadLeptonViaTauToken_(consumes<std::vector<int> >(iConfig.getParameter<edm::InputTag>("genBHadLeptonViaTau"))),
genCHadJetIndexToken_(consumes<std::vector<int> >(iConfig.getParameter<edm::InputTag>("genCHadJetIndex"))),
genCHadFlavourToken_(consumes<std::vector<int> >(iConfig.getParameter<edm::InputTag>("genCHadFlavour"))),
genCHadFromTopWeakDecayToken_(consumes<std::vector<int> >(iConfig.getParameter<edm::InputTag>("genCHadFromTopWeakDecay"))),
genCHadBHadronIdToken_(consumes<std::vector<int> >(iConfig.getParameter<edm::InputTag>("genCHadBHadronId")))
{
    produces<int>("genTtbarId");
}


//
// member functions
//

// ------------ method called to produce the data  ------------
void
GenTtbarCategorizer::produce(edm::Event& iEvent, const edm::EventSetup& iSetup)
{
    // Access gen jets
    edm::Handle<reco::GenJetCollection> genJets;
    iEvent.getByToken(genJetsToken_, genJets);
    
    
    // Access B hadrons information
    edm::Handle<std::vector<int> > genBHadFlavour;
    iEvent.getByToken(genBHadFlavourToken_, genBHadFlavour);
    
    edm::Handle<std::vector<int> > genBHadJetIndex;
    iEvent.getByToken(genBHadJetIndexToken_, genBHadJetIndex);
    
    edm::Handle<std::vector<int> > genBHadFromTopWeakDecay;
    iEvent.getByToken(genBHadFromTopWeakDecayToken_, genBHadFromTopWeakDecay);
    
    edm::Handle<std::vector<reco::GenParticle> > genBHadPlusMothers;
    iEvent.getByToken(genBHadPlusMothersToken_, genBHadPlusMothers);
    
    edm::Handle<std::vector<std::vector<int> > > genBHadPlusMothersIndices;
    iEvent.getByToken(genBHadPlusMothersIndicesToken_, genBHadPlusMothersIndices);
    
    edm::Handle<std::vector<int> > genBHadIndex;
    iEvent.getByToken(genBHadIndexToken_, genBHadIndex);
    
    edm::Handle<std::vector<int> > genBHadLeptonHadronIndex;
    iEvent.getByToken(genBHadLeptonHadronIndexToken_, genBHadLeptonHadronIndex);
    
    edm::Handle<std::vector<int> > genBHadLeptonViaTau;
    iEvent.getByToken(genBHadLeptonViaTauToken_, genBHadLeptonViaTau);
    
    
    // Access C hadrons information
    edm::Handle<std::vector<int> > genCHadFlavour;
    iEvent.getByToken(genCHadFlavourToken_, genCHadFlavour);
    
    edm::Handle<std::vector<int> > genCHadJetIndex;
    iEvent.getByToken(genCHadJetIndexToken_, genCHadJetIndex);
    
    edm::Handle<std::vector<int> > genCHadFromTopWeakDecay;
    iEvent.getByToken(genCHadFromTopWeakDecayToken_, genCHadFromTopWeakDecay);
    
    edm::Handle<std::vector<int> > genCHadBHadronId;
    iEvent.getByToken(genCHadBHadronIdToken_, genCHadBHadronId);
    
    
    // Map <jet index, number of specific hadrons in jet>
    // B jets with b hadrons directly from top quark decay
    std::map<int, int> bJetFromTopIds;
    // B jets with b hadrons before top quark decay chain
    std::map<int, int> bJetIds;
    // C jets with c hadrons before top quark decay chain
    std::map<int, int> cJetIds;
    
    
    // Count number of specific b hadrons in each jet
    for(size_t hadronId = 0; hadronId < genBHadIndex->size(); ++hadronId) {
        // Flavour of the hadron's origin
        const int flavour = genBHadFlavour->at(hadronId);
        // Skip b hadrons coming for W decays (rare case due to CKM, in some generated samples not present at all)
        if(std::abs(flavour) == 24) continue;
        // Index of jet associated to the hadron
        const int jetIndex = genBHadJetIndex->at(hadronId);
        // Skip hadrons which have no associated jet
        if(jetIndex < 0) continue;
        // Skip if jet is not in acceptance
        if(genJets->at(jetIndex).pt() < genJetPtMin_) continue;
        if(std::fabs(genJets->at(jetIndex).eta()) > genJetAbsEtaMax_) continue;
        // Jet from direct top quark decay [pdgId(top)=6]
        if(std::abs(flavour) == 6) {
            if(bJetFromTopIds.count(jetIndex) < 1) bJetFromTopIds[jetIndex] = 1;
            else bJetFromTopIds[jetIndex]++;
            continue;
        }
        // Identify jets with b hadrons not from top quark decay
        if(bJetIds.count(jetIndex) < 1) bJetIds[jetIndex] = 1;
        else bJetIds[jetIndex]++;
    }
    
    // Count number of specific c hadrons in each c jet
    for(size_t hadronId = 0; hadronId < genCHadJetIndex->size(); ++hadronId) {
        // Skip c hadrons that are coming from b hadrons
        if(genCHadBHadronId->at(hadronId) >= 0) continue;
        // Skip c hadrons coming from W decays
        if(std::abs(genCHadFlavour->at(hadronId)) == 24) continue;
        // Index of jet associated to the hadron
        const int jetIndex = genCHadJetIndex->at(hadronId);
        // Skip hadrons which have no associated jet
        if(jetIndex < 0) continue;
        // Skip if jet is not in acceptance
        if(genJets->at(jetIndex).pt() < genJetPtMin_) continue;
        if(std::fabs(genJets->at(jetIndex).eta()) > genJetAbsEtaMax_) continue;
        // Identify jets with b hadrons
        if(cJetIds.count(jetIndex) < 1) cJetIds[jetIndex] = 1;
        else cJetIds[jetIndex]++;
    }
    
    // Find additional b jets
    std::vector<int> additionalBJetIds;
    for(std::map<int, int>::iterator it = bJetIds.begin(); it != bJetIds.end(); ++it) {
        const int jetId = it->first;
        // Skip jet if it contains a b hadron directly from top quark decay
        if(bJetFromTopIds.count(jetId) > 0) continue;
        additionalBJetIds.push_back(jetId);
    }
    
    // Find additional c jets
    std::vector<int> additionalCJetIds;
    for(std::map<int, int>::iterator it = cJetIds.begin(); it != cJetIds.end(); ++it) {
        const int jetId = it->first;
        // Skip jet if it contains a b hadron, thus being a b jet
        if(bJetFromTopIds.count(jetId) > 0) continue;
        additionalCJetIds.push_back(jetId);
    }
    
    
    // Categorize event based on number of additional b/c jets
    // and number of corresponding hadrons in each of them
    int additionalJetEventId = bJetFromTopIds.size()*100;
    // tt + 1 additional b jet
    if(additionalBJetIds.size() == 1){
        int nHadronsInJet = bJetIds[additionalBJetIds.at(0)];
        // tt + 1 additional b jet from 1 additional b hadron
        if(nHadronsInJet == 1) additionalJetEventId += 51;
        // tt + 1 additional b jet from >=2 additional b hadrons
        else additionalJetEventId += 52;
    }
    // tt + >=2 additional b jets
    else if(additionalBJetIds.size() > 1){
        // Check first two additional b jets (rare cases could have more)
        int nHadronsInJet1 = bJetIds[additionalBJetIds.at(0)];
        int nHadronsInJet2 = bJetIds[additionalBJetIds.at(1)];
        // tt + >=2 additional b jets each from 1 additional b hadron
        if(std::max(nHadronsInJet1, nHadronsInJet2) == 1) additionalJetEventId += 53;
        // tt + >=2 additional b jets one of which from >=2 additional b hadrons
        else if(std::min(nHadronsInJet1, nHadronsInJet2) == 1 && std::max(nHadronsInJet1, nHadronsInJet2) > 1) additionalJetEventId += 54;
        // tt + >=2 additional b jets each from >=2 additional b hadrons
        else if(std::min(nHadronsInJet1, nHadronsInJet2) > 1) additionalJetEventId += 55;
    }
    // tt + no additional b jets
    else{
        // tt + 1 additional c jet
        if(additionalCJetIds.size() == 1){
            int nHadronsInJet = cJetIds[additionalCJetIds.at(0)];
            // tt + 1 additional c jet from 1 additional c hadron
            if(nHadronsInJet == 1) additionalJetEventId += 41;
            // tt + 1 additional c jet from >=2 additional c hadrons
            else additionalJetEventId += 42;
        }
        // tt + >=2 additional c jets
        else if(additionalCJetIds.size() > 1){
            // Check first two additional c jets (rare cases could have more)
            int nHadronsInJet1 = cJetIds[additionalCJetIds.at(0)];
            int nHadronsInJet2 = cJetIds[additionalCJetIds.at(1)];
            // tt + >=2 additional c jets each from 1 additional c hadron
            if(std::max(nHadronsInJet1, nHadronsInJet2) == 1) additionalJetEventId += 43;
            // tt + >=2 additional c jets one of which from >=2 additional c hadrons
            else if(std::min(nHadronsInJet1, nHadronsInJet2) == 1 && std::max(nHadronsInJet1, nHadronsInJet2) > 1) additionalJetEventId += 44;
            // tt + >=2 additional c jets each from >=2 additional c hadrons
            else if(std::min(nHadronsInJet1, nHadronsInJet2) > 1) additionalJetEventId += 45;
        }
        // tt + no additional c jets
        else{
            // tt + light jets
            additionalJetEventId += 0;
        }
    }
    
    std::auto_ptr<int> ttbarId(new int);
    *ttbarId = additionalJetEventId;
    iEvent.put(ttbarId, "genTtbarId");
    
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
GenTtbarCategorizer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(GenTtbarCategorizer);
