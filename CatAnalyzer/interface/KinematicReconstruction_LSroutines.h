#ifndef KinematicReconstruction_LSroutines_h
#define KinematicReconstruction_LSroutines_h

#include <vector>

#include <TLorentzVector.h>

class TH1F;
class TF1;
class TF2;

#include "classesFwd.h"





class KinematicReconstruction_LSroutines{

public:
    KinematicReconstruction_LSroutines();
    KinematicReconstruction_LSroutines(const double& mass_Wp, const double& mass_Wm);
    KinematicReconstruction_LSroutines(const double& mass_top, const double& mass_topbar, 
                                       const double& mass_b, const double& mass_bbar, 
                                       const double& mass_w, const double& mass_wbar, 
                                       const double& mass_al, const double& mass_l);
    ~KinematicReconstruction_LSroutines();
    
    void fDelete()const;
    
    
    void ini(const double& mass_Wp, const double& mass_Wm);
    
    void setConstraints(const TLorentzVector& LV_al, 
                        const TLorentzVector& LV_l, 
                        const TLorentzVector& LV_b, 
                        const TLorentzVector& LV_bbar,
                        const double& missPx,
                        const double& missPy
                       );
    
    void setConstraints(const LV& LV_al, 
                        const LV& LV_l, 
                        const LV& LV_b, 
                        const LV& LV_bbar,
                        const double& missPx,
                        const double& missPy
                       );
    
    int getNsol()const;
    
    struct TopSolution{
        double dR;
        double dN;
        TLorentzVector top;
        TLorentzVector topbar;
        TLorentzVector neutrino;
        TLorentzVector neutrinobar;
        TLorentzVector wp;
        TLorentzVector wm;
        double x1;
        double x2;
        double mtt;
        double weight;
        
    };
    const std::vector<TopSolution>* getTtSol()const;
    
    void setTrueInfo(const LV& LV_Top,const LV& LV_AntiTop,const LV& LV_Neutrino, const LV& LV_AntiNeutrino);
    void sortBy(std::string ch);
    void print()const;
    
private:
    void filldR();
    void filldN();
    void swapTopSol(TopSolution& sol1, TopSolution& sol2)const;
    void sortTopSol(std::vector<TopSolution>& v)const;
    void doAll();
    void topRec(const double& px_neutrino);
    void findCoeff(double* const koeficienty);
    void quartic_equation(const double& h0, const double& h1, const double& h2, const double& h3, const double& h4, std::vector<double>& v)const;
    void cubic_equation(const double& a, const double& b, const double& c, const double& d, std::vector<double> &v)const;
    void quadratic_equation(const double& a, const double& b, const double& c, std::vector<double>& v)const;
    void linear_equation(const double& a, const double& b, std::vector<double>& v)const;
    int sign(const long double& ld)const;
    double landau2D(const double& x, const double& y)const;
    
    
    //Utility Methods
    double sqr(const double& x)const;
    void swap(double& realone, double& realtwo)const;
    
    int nSol_;
    double coeffs_[5];
    std::vector<double> vect_pxv_;
    std::vector<TopSolution> ttSol_;
    TLorentzVector al_;
    TLorentzVector l_;
    TLorentzVector b_;
    TLorentzVector bbar_;
    TLorentzVector top_;
    TLorentzVector topbar_;
    TLorentzVector neutrino_;
    TLorentzVector neutrinobar_;
    TLorentzVector w_;
    TLorentzVector wbar_;
    TLorentzVector tt_;
    
    TLorentzVector true_top_;
    TLorentzVector true_topbar_;
    TLorentzVector true_neutrino_;
    TLorentzVector true_neutrinobar_;
    
    double px_miss_;
    double py_miss_;
    
    double mt_;
    double mtbar_;
    double mb_;
    double mbbar_;
    double mw_;
    double mwbar_;
    double ml_;
    double mal_;
    double mv_;
    double mav_;
    
    double a1_,a2_,a3_,a4_;
    double b1_,b2_,b3_,b4_;
    double c22_,c21_,c20_,c11_,c10_,c00_;
    double d22_,d21_,d20_,d11_,d10_,d00_;
    double d0_,d1_,d2_;
    double c0_,c1_,c2_;
    
};


#endif
