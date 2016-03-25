#ifndef KinematicReconstruction_MeanSol_h
#define KinematicReconstruction_MeanSol_h

#include <vector>

#include <TLorentzVector.h>


#include "classesFwd.h"




class KinematicReconstruction_MeanSol{
    
public:
    
    KinematicReconstruction_MeanSol(const double& topm);
    ~KinematicReconstruction_MeanSol();  
    void add(const TLorentzVector& top, const TLorentzVector& topbar, const TLorentzVector& n, const TLorentzVector& nbar, const double& weight, const double& mbl_weight);
    void add(const TLorentzVector& top, const TLorentzVector& topbar, const TLorentzVector& n, const TLorentzVector& nbar, const double& weight);
    
    void getMeanVect(TLorentzVector& lv, const std::vector<TLorentzVector>& vlv, const double& mass)const;
    void getMeanSol(TLorentzVector& top, TLorentzVector& topbar, TLorentzVector& n, TLorentzVector& nbar)const;
    double getSumWeight()const;
    int getNsol()const;
    void clear();
    
    
    // FIXME: right now this is only using getMeanSol(), since the whole class is based on TLorentzVector
    /// Get Lorentz vectors of (anti-)top and (anti-)neutrino for mean solution
    void meanSolution(LV& top, LV& antiTop, LV& neutrino, LV& antiNeutrino)const;
    
    
    
private:
    
    std::vector<TLorentzVector> v_top_;
    std::vector<TLorentzVector> v_topbar_;
    std::vector<TLorentzVector> v_n_;
    std::vector<TLorentzVector> v_nbar_;
    
    std::vector<double> v_weight_;
    double sum_weight_;  
    double max_sum_weight_;
    
    const double mass_top_;
};



#endif


