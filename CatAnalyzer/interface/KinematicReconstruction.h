#ifndef KinematicReconstruction_h
#define KinematicReconstruction_h

#include <vector>
#include <map>

#include <TLorentzVector.h>

class TH1;
class TRandom3;

#include "classes.h"
#include "sampleHelpers.h"

class KinematicReconstructionSolution;
class KinematicReconstructionSolutions;
class KinematicReconstruction_MeanSol;


struct Struct_KinematicReconstruction{
    TLorentzVector lp, lm;
    TLorentzVector jetB, jetBbar;
    size_t jetB_index, jetBbar_index;
    TLorentzVector met;
    TLorentzVector Wplus, Wminus;
    TLorentzVector top, topBar;
    TLorentzVector neutrino, neutrinoBar;
    TLorentzVector ttbar;
    double recMtop;
    double weight;
    int ntags;
};

class KinematicReconstruction{
    
public:
    
    KinematicReconstruction(const int minNumberOfBtags, const bool preferBtags, const bool massLoop =false);
    ~KinematicReconstruction(){}
    
    int getNSol()const;
    Struct_KinematicReconstruction getSol()const;
    std::vector<Struct_KinematicReconstruction> getSols()const;
    
    void loadData();
    void kinReco(const LV& leptonMinus, const LV& leptonPlus, const VLV* jets, const std::vector<double>* btags, const LV* met);
    void kinRecoMassLoop(const LV& leptonMinus, const LV& leptonPlus, const VLV* jets, const std::vector<double>* btags, const LV* met);
    
    /// Retrieve all solutions valid for setup of kinematic reconstruction
    KinematicReconstructionSolutions solutions(const std::vector<int>& leptonIndices, const std::vector<int>& antiLeptonIndices,
                                               const std::vector<int>& jetIndices, const std::vector<int>& bjetIndices,
                                               const VLV& allLeptons,
                                               const VLV& allJets, const std::vector<double>& btags,
                                               const LV& met)const;
    
    
private:
    
    /// Calculate solution for specific lepton, antilepton and pair of jets
    std::vector<KinematicReconstructionSolution> solutionsPerObjectCombination(const int leptonIndex, const int antiLeptonIndex,
                                                                               const int jetIndex1, const int jetIndex2,
                                                                               const VLV& allLeptons,
                                                                               const VLV& allJets, const std::vector<double>& btags,
                                                                               const LV& met,
                                                                               const int numberOfBtags)const;
    
    /// Set seeds for random number generators
    void setRandomNumberSeeds(const LV& lepton, const LV& antiLepton, const LV& jet1, const LV& jet2)const;
    
    /// Minimum number of b-tags required for solutions (0, 1, 2)
    const int minNumberOfBtags_;
    
    /// Prefer solutions with b-tags (2 tags if existing, else 1 tag if existing, else 0 tags)
    const bool preferBtags_;
    
    /// Whether to run mass loop for top mass, instead of smearings according to uncertainties
    const bool massLoop_;
    
    
    
    // FIXME: temporary helper variables for cleanup
    void setSolutions();
    
    void inputNoJetMerging(std::vector<int>& b1_id, std::vector<int>& b2_id, std::vector<int>& nb_tag,
                           const std::vector<double>& btags)const;
    
    /// Calculate solution using smearing
    bool solutionSmearing(KinematicReconstruction_MeanSol& meanSolution,
                          const LV& lepton, const LV& antiLepton,
                          const LV& jet1, const LV& jet2,
                          const LV& met)const;
    
    
    
    void angle_rot(const double& alpha, const double& e, const TLorentzVector& inJet, TLorentzVector& jet_sm)const;
    
    TRandom3* r3_;
    
    int nSol_;
    Struct_KinematicReconstruction sol_;
    std::vector<Struct_KinematicReconstruction> sols_;
    
    // W mass
    TH1* h_wmass_;

    // jet resolution
    TH1* h_jetAngleRes_;
    TH1* h_jetEres_;
       
    //lepton resolution
    TH1* h_lepAngleRes_;
    TH1* h_lepEres_;
    
// mbl
    TH1* h_mbl_w_;
};


class KinematicReconstructionScaleFactors{
    
public:
    
    /// Constructor
    KinematicReconstructionScaleFactors(const std::vector<Channel::Channel>& channels,
                                        const Systematic::Systematic& systematic);
    
    /// Destructor
    ~KinematicReconstructionScaleFactors(){}
    
    
    
    /// Prepare scale factors per channel
    void prepareChannel(const Channel::Channel& channel);
    
    /// Get kinematic reconstruction per-event scale factor
    double getSF()const;
    
    
    
private:
    
    /// Enumeration for possible systematics
    enum SystematicInternal{nominal, vary_up, vary_down, undefined};
    
    /** Set up the scale factor for the Kinematic Reconstruction
     *
     * Currently a flat per-channel SF is used. For the systematic KIN_UP and KIN_DOWN,
     * the SF is modified by its uncertainty.
     *
     * To calculate the SF, you need to set the SF to 1 and rerun. Then determine the
     * SF with kinRecoEfficienciesAndSF
     */
    void prepareSF(const SystematicInternal& systematic);
    
    
    
    /// Map containing the flat scale factors for all channels
    std::map<Channel::Channel, double> m_scaleFactor_;
    
    /// The per-channel scale factor
    double scaleFactor_;
};







#endif



