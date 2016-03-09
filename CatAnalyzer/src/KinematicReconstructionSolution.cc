#include <iostream>
#include <iomanip>

#include "CATTools/CatAnalyzer/interface/KinematicReconstructionSolution.h"







// -------------------------------------- Methods for KinematicReconstructionSolution --------------------------------------


KinematicReconstructionSolution::KinematicReconstructionSolution():
allLeptons_(0),
allJets_(0),
leptonIndex_(-1),
antiLeptonIndex_(-1),
bjetIndex_(-1),
antiBjetIndex_(-1),
top_(LV(0.,0.,0.,0.)),
antiTop_(LV(0.,0.,0.,0.)),
neutrino_(LV(0.,0.,0.,0.)),
antiNeutrino_(LV(0.,0.,0.,0.)),
reconstructedTopMass_(-999.),
numberOfBtags_(-1),
isNoSmearSol_(false)
{}



KinematicReconstructionSolution::KinematicReconstructionSolution(const KinematicReconstructionSolution& kinematicReconstructionSolution):
allLeptons_(kinematicReconstructionSolution.allLeptons_),
allJets_(kinematicReconstructionSolution.allJets_),
leptonIndex_(kinematicReconstructionSolution.leptonIndex_),
antiLeptonIndex_(kinematicReconstructionSolution.antiLeptonIndex_),
bjetIndex_(kinematicReconstructionSolution.bjetIndex_),
antiBjetIndex_(kinematicReconstructionSolution.antiBjetIndex_),
top_(kinematicReconstructionSolution.top_),
antiTop_(kinematicReconstructionSolution.antiTop_),
neutrino_(kinematicReconstructionSolution.neutrino_),
antiNeutrino_(kinematicReconstructionSolution.antiNeutrino_),
reconstructedTopMass_(kinematicReconstructionSolution.reconstructedTopMass_),
numberOfBtags_(kinematicReconstructionSolution.numberOfBtags_),
m_weight_(kinematicReconstructionSolution.m_weight_),
isNoSmearSol_(kinematicReconstructionSolution.isNoSmearSol_)
{}



KinematicReconstructionSolution::KinematicReconstructionSolution(const VLV* const v_allLeptons, const VLV* const v_allJets,
                                                                 const int leptonIndex, const int antiLeptonIndex,
                                                                 const int bjetIndex, const int antiBjetIndex,
                                                                 const LV& top, const LV& antiTop,
                                                                 const LV& neutrino, const LV& antiNeutrino,
                                                                 const double& reconstructedTopMass,
                                                                 const int numberOfBtags,
                                                                 const std::map<WeightType, double>& m_weight,
                                                                 const bool isNoSmearSol):
allLeptons_(v_allLeptons),
allJets_(v_allJets),
leptonIndex_(leptonIndex),
antiLeptonIndex_(antiLeptonIndex),
bjetIndex_(bjetIndex),
antiBjetIndex_(antiBjetIndex),
top_(top),
antiTop_(antiTop),
neutrino_(neutrino),
antiNeutrino_(antiNeutrino),
reconstructedTopMass_(reconstructedTopMass),
numberOfBtags_(numberOfBtags),
m_weight_(m_weight),
isNoSmearSol_(isNoSmearSol)
{}



void KinematicReconstructionSolution::print()const
{
    const LV& lepton(this->lepton());
    const LV& antiLepton(this->antiLepton());
    const LV& bjet(this->bjet());
    const LV& antiBjet(this->antiBjet());
    const LV ttbar = this->ttbar();
    const LV wMinus = this->wMinus();
    const LV wPlus = this->wPlus();
    const LV top = bjet + wPlus;
    const LV antiTop = antiBjet + wMinus;
    const LV ttbarReco = top + antiTop;
    
    std::cout<<"Solution of kinematic reconstruction:\n";
    std::cout<<"\tStored quantities\n";
    std::cout<<"\t\tNumber of b-tags:                   "<<numberOfBtags_<<"\n"
             <<"\t\tIndex (lepton, antilepton):         "<<leptonIndex_<<" , "<<antiLeptonIndex_<<"\n"
             <<"\t\tIndex (b jet, anti-b jet):          "<<bjetIndex_<<" , "<<antiBjetIndex_<<"\n"
             <<std::fixed<<std::setprecision(3)
             <<"\t\tReconstructed top mass:             "<<reconstructedTopMass_<<"\n"
             <<"\t\tNeutrino (pt, eta, phi, mass):      "<<neutrino_.pt()<<" , "<<neutrino_.eta()<<" , "<<neutrino_.phi()<<" , "<<neutrino_.mass()<<"\n"
             <<"\t\tAnti-neutrino (pt, eta, phi, mass): "<<antiNeutrino_.pt()<<" , "<<antiNeutrino_.eta()<<" , "<<antiNeutrino_.phi()<<" , "<<antiNeutrino_.mass()<<"\n"
             <<"\t\tTop (pt, eta, phi, mass):           "<<top_.pt()<<" , "<<top_.eta()<<" , "<<top_.phi()<<" , "<<top_.mass()<<"\n"
             <<"\t\tAnti-top (pt, eta, phi, mass):      "<<antiTop_.pt()<<" , "<<antiTop_.eta()<<" , "<<antiTop_.phi()<<" , "<<antiTop_.mass()<<"\n"
             <<"\t\ttt system (pt, eta, phi, mass):     "<<ttbar.pt()<<" , "<<ttbar.eta()<<" , "<<ttbar.phi()<<" , "<<ttbar.mass()<<"\n";
    std::cout<<"\tInput objects\n";
    std::cout<<"\t\tLepton (pt, eta, phi, mass):        "<<lepton.pt()<<" , "<<lepton.eta()<<" , "<<lepton.phi()<<" , "<<lepton.mass()<<"\n"
             <<"\t\tAnti-lepton (pt, eta, phi, mass):   "<<antiLepton.pt()<<" , "<<antiLepton.eta()<<" , "<<antiLepton.phi()<<" , "<<antiLepton.mass()<<"\n"
             <<"\t\tb jet (pt, eta, phi, mass):         "<<bjet.pt()<<" , "<<bjet.eta()<<" , "<<bjet.phi()<<" , "<<bjet.mass()<<"\n"
             <<"\t\tAnti-b jet (pt, eta, phi, mass):    "<<antiBjet.pt()<<" , "<<antiBjet.eta()<<" , "<<antiBjet.phi()<<" , "<<antiBjet.mass()<<"\n";
    std::cout<<"\tFrom input objects reconstructed objects\n";
    std::cout<<"\t\tW+ (pt, eta, phi, mass):            "<<wPlus.pt()<<" , "<<wPlus.eta()<<" , "<<wPlus.phi()<<" , "<<wPlus.mass()<<"\n"
             <<"\t\tW- jet (pt, eta, phi, mass):        "<<wMinus.pt()<<" , "<<wMinus.eta()<<" , "<<wMinus.phi()<<" , "<<wMinus.mass()<<"\n"
             <<"\t\tTop (pt, eta, phi, mass):           "<<top.pt()<<" , "<<top.eta()<<" , "<<top.phi()<<" , "<<top.mass()<<"\n"
             <<"\t\tAnti-top (pt, eta, phi, mass):      "<<antiTop.pt()<<" , "<<antiTop.eta()<<" , "<<antiTop.phi()<<" , "<<antiTop.mass()<<"\n"
             <<"\t\ttt system (pt, eta, phi, mass):     "<<ttbarReco.pt()<<" , "<<ttbarReco.eta()<<" , "<<ttbarReco.phi()<<" , "<<ttbarReco.mass()<<"\n";
    std::cout<<"\tNumber of weights stored: "<<m_weight_.size()<<"\n";
    for(const auto& weight : m_weight_) std::cout<<"\t\tWeight (enum ID, value):            "<<weight.first<<" , "<<weight.second<<"\n";
}










// -------------------------------------- Methods for KinematicReconstructionSolutions --------------------------------------


KinematicReconstructionSolutions::KinematicReconstructionSolutions()
{}



void KinematicReconstructionSolutions::addSolutions(const std::vector<KinematicReconstructionSolution>& solutions)
{
    for(const auto& solution : solutions) this->addSolution(solution);
}



void KinematicReconstructionSolutions::addSolution(const KinematicReconstructionSolution& solution)
{
    // Initialise all maps of weight-ordered solution indices, in first solution only
    if(!m_weightIndex_.size()){
        for(const auto weightTypeWeight : solution.weightMap()){
            const auto weightType = weightTypeWeight.first;
            m_weightIndex_[weightType] = std::vector<size_t>();
            m_weightIndexTwoBtags_[weightType] = std::vector<size_t>();
            m_weightIndexOneBtag_[weightType] = std::vector<size_t>();
            m_weightIndexNoBtags_[weightType] = std::vector<size_t>();
        }
    }
    
    // Set pointers to the specific b-tag multiplicity category
    std::vector<size_t>* v_solutionByCategory(0);
    std::map<KinematicReconstructionSolution::WeightType, std::vector<size_t>>* m_weightIndexByCategory(0);
    const int numberOfBtags = solution.numberOfBtags();
    if(numberOfBtags == 2){
        v_solutionByCategory = &v_solutionTwoBtags_;
        m_weightIndexByCategory = &m_weightIndexTwoBtags_;
    }
    else if(numberOfBtags == 1){
        v_solutionByCategory = &v_solutionOneBtag_;
        m_weightIndexByCategory = &m_weightIndexOneBtag_;
    }
    else if(numberOfBtags == 0){
        v_solutionByCategory = &v_solutionNoBtags_;
        m_weightIndexByCategory = &m_weightIndexNoBtags_;
    }
    else{
        std::cerr<<"ERROR in KinematicReconstructionSolutions::addSolution()! Invalid number of b-tags: "<<numberOfBtags
                 <<"\n...break\n"<<std::endl;
        exit(731);
    }
    
    // Add solution to all solutions, and to specific b-tag multiplicity category
    v_solution_.push_back(solution);
    const size_t solutionIndex = std::distance(v_solution_.begin(), v_solution_.end()) - 1;
    v_solutionByCategory->push_back(solutionIndex);
    
    // Fill for each weight type the indices of solutions, ordered for the weight
    for(const auto weightTypeWeight : solution.weightMap()){
        const KinematicReconstructionSolution::WeightType weightType = weightTypeWeight.first;
        const double weight = weightTypeWeight.second;
        this->insertIndex(solutionIndex, weight, m_weightIndex_.at(weightType));
        this->insertIndexByCategory(*v_solutionByCategory, weight, m_weightIndexByCategory->at(weightType));
    }
}



void KinematicReconstructionSolutions::insertIndex(const size_t solutionIndex,
                                                   const double weight, std::vector<size_t>& v_index)const
{
    if(!v_index.size()){
        v_index.push_back(solutionIndex);
        return;
    }
    
    bool isInserted(false);
    for(std::vector<size_t>::iterator i_index = v_index.begin(); i_index != v_index.end(); ++i_index){
        const double& iWeight = v_solution_.at(*i_index).weight();
        if(iWeight < weight){
            v_index.insert(i_index, solutionIndex);
            isInserted = true;
            break;
        }
    }
    if(!isInserted) v_index.push_back(solutionIndex);
}



void KinematicReconstructionSolutions::insertIndexByCategory(const std::vector<size_t>& v_solutionIndex,
                                                             const double weight, std::vector<size_t>& v_solutionIndexByCategory)const
{
    const size_t solutionIndexByCategory = std::distance(v_solutionIndex.begin(), v_solutionIndex.end()) - 1;
    
    if(!v_solutionIndexByCategory.size()){
        v_solutionIndexByCategory.push_back(solutionIndexByCategory);
        return;
    }
    
    bool isInserted(false);
    for(std::vector<size_t>::iterator i_index = v_solutionIndexByCategory.begin(); i_index != v_solutionIndexByCategory.end(); ++i_index){
        const double& iWeight = v_solution_.at(v_solutionIndex.at(*i_index)).weight();
        if(iWeight < weight){
            v_solutionIndexByCategory.insert(i_index, solutionIndexByCategory);
            isInserted = true;
            break;
        }
    }
    if(!isInserted) v_solutionIndexByCategory.push_back(solutionIndexByCategory);
}



const KinematicReconstructionSolution& KinematicReconstructionSolutions::solution(const KinematicReconstructionSolution::WeightType weightType,
                                                                                  const size_t solutionNumber)const
{
    const size_t index = m_weightIndex_.at(weightType).at(solutionNumber);
    return v_solution_.at(index);
}



const KinematicReconstructionSolution& KinematicReconstructionSolutions::solutionTwoBtags(const KinematicReconstructionSolution::WeightType weightType,
                                                                                          const size_t solutionNumber)const
{
    const size_t index = m_weightIndexTwoBtags_.at(weightType).at(solutionNumber);
    return v_solution_.at(v_solutionTwoBtags_.at(index));
}



const KinematicReconstructionSolution& KinematicReconstructionSolutions::solutionOneBtag(const KinematicReconstructionSolution::WeightType weightType,
                                                                                         const size_t solutionNumber)const
{
    const size_t index = m_weightIndexOneBtag_.at(weightType).at(solutionNumber);
    return v_solution_.at(v_solutionOneBtag_.at(index));
}



const KinematicReconstructionSolution& KinematicReconstructionSolutions::solutionNoBtags(const KinematicReconstructionSolution::WeightType weightType,
                                                                                         const size_t solutionNumber)const
{
    const size_t index = m_weightIndexNoBtags_.at(weightType).at(solutionNumber);
    return v_solution_.at(v_solutionNoBtags_.at(index));
}






