#include <vector>
#include <iostream>

#include <TLorentzVector.h>
#include <TH1F.h>
#include <TF1.h>
#include <TF2.h>
#include <TAttLine.h>
#include <TFile.h>
#include <TMath.h>

#include "CATTools/CatAnalyzer/interface/KinematicReconstruction_LSroutines.h"
#include "CATTools/CatAnalyzer/interface/classes.h"




constexpr double TopMASS = 172.5;








KinematicReconstruction_LSroutines::KinematicReconstruction_LSroutines()
{
    mt_    = TopMASS;
    mtbar_ = TopMASS;
    mb_    = 4.8;
    mbbar_ = 4.8;
    mw_    = 80.4;
    mwbar_ = 80.4; 
    ml_    = 0.;
    mal_    = 0.; 
    mv_=0;
    mav_=0;
}



KinematicReconstruction_LSroutines::KinematicReconstruction_LSroutines(const double& mass_Wp, const double& mass_Wm)
{
    mt_    = TopMASS;
    mtbar_ = TopMASS;
    mb_    = 4.8;
    mbbar_ = 4.8;
    mw_    = mass_Wp;
    mwbar_ = mass_Wm; 
    ml_    = 0.;
    mal_    = 0.;
    mv_=0;
    mav_=0;
}



KinematicReconstruction_LSroutines::KinematicReconstruction_LSroutines(const double& mass_top, const double& mass_topbar, 
                                                                       const double& mass_b, const double& mass_bbar, 
                                                                       const double& mass_Wp, const double& mass_Wm, 
                                                                       const double& mass_al, const double& mass_l)
{
    mt_    = mass_top;
    mtbar_ = mass_topbar;
    mb_    = mass_b;
    mbbar_ = mass_bbar;
    mw_    = mass_Wp;
    mwbar_ = mass_Wm;
    mal_=mass_al;
    ml_=mass_l;
    mv_=0;
    mav_=0;
}



KinematicReconstruction_LSroutines::~KinematicReconstruction_LSroutines()
{
    
}



void KinematicReconstruction_LSroutines::fDelete()const
{
    delete this;
}



void KinematicReconstruction_LSroutines::ini(const double& mass_Wp, const double& mass_Wm)
{
    mt_    = TopMASS;
    mtbar_ = TopMASS;
    mb_    = 4.8;
    mbbar_ = 4.8;
    mw_    = mass_Wp;
    mwbar_ = mass_Wm; 
    ml_    = 0.;
    mal_    = 0.;
    mv_=0;
    mav_=0;
}



// FIXME: following two functions do almost the same thing --> use the one in the other, such that implementation exists only once
void KinematicReconstruction_LSroutines::setConstraints(const TLorentzVector& LV_al, const TLorentzVector& LV_l, const TLorentzVector& LV_b, const TLorentzVector& LV_bbar, const double& missPx, const double& missPy)
{
    l_  = LV_l;
    al_ = LV_al;
    b_  = LV_b;
    bbar_ = LV_bbar;
    px_miss_ = missPx;
    py_miss_ = missPy;
    this->doAll();
}



void KinematicReconstruction_LSroutines::setConstraints(const LV& LV_al, const LV& LV_l, const LV& LV_b, const LV& LV_bbar, const double& missPx, const double& missPy)
{
    TLorentzVector temp_al(LV_al.Px(),LV_al.Py(),LV_al.Pz(),LV_al.E());
    TLorentzVector temp_l(LV_l.Px(),LV_l.Py(),LV_l.Pz(),LV_l.E());
    TLorentzVector temp_b(LV_b.Px(),LV_b.Py(),LV_b.Pz(),LV_b.E());
    TLorentzVector temp_bbar(LV_bbar.Px(),LV_bbar.Py(),LV_bbar.Pz(),LV_bbar.E());
    l_  = temp_l;
    al_ = temp_al;
    b_  = temp_b;
    bbar_ = temp_bbar;
    px_miss_ = missPx;
    py_miss_ = missPy;
    this->doAll();
}



int KinematicReconstruction_LSroutines::getNsol()const
{
    return nSol_;
}




const std::vector<KinematicReconstruction_LSroutines::TopSolution>* KinematicReconstruction_LSroutines::getTtSol()const
{
    return &ttSol_;
}



void KinematicReconstruction_LSroutines::setTrueInfo(const LV& LV_Top, const LV& LV_AntiTop,const LV& LV_Neutrino,const LV& LV_AntiNeutrino)
{
    true_top_         = TLorentzVector(LV_Top.Px(),LV_Top.Py(),LV_Top.Pz(),LV_Top.E());  
    true_topbar_      = TLorentzVector(LV_AntiTop.Px(),LV_AntiTop.Py(),LV_AntiTop.Pz(),LV_AntiTop.E());  
    true_neutrino_    = TLorentzVector(LV_Neutrino.Px(),LV_Neutrino.Py(),LV_Neutrino.Pz(),LV_Neutrino.E());  
    true_neutrinobar_ = TLorentzVector(LV_AntiNeutrino.Px(),LV_AntiNeutrino.Py(),LV_AntiNeutrino.Pz(),LV_AntiNeutrino.E());  
    filldR();
    filldN();
}



void KinematicReconstruction_LSroutines::print()const
{
    for(int i=0;i<(int)ttSol_.size();++i)
    {
        printf("\nSol: %d:   weight: %f dN: %f\n",i+1,ttSol_[i].weight,ttSol_[i].dN);
        ttSol_[i].top.Print();
        ttSol_[i].topbar.Print();
        ttSol_[i].neutrino.Print();
        ttSol_[i].neutrinobar.Print();
    }
}



double KinematicReconstruction_LSroutines::landau2D(const double& xv, const double& xvbar)const
{
    return 30.641*TMath::Landau(xv,57.941,22.344,0)*TMath::Landau(xvbar,57.533,22.232,0);   // top_mass = 172.5 GeV  CME = 7 TeV
}



void KinematicReconstruction_LSroutines::doAll()
{
    
        this->findCoeff(coeffs_);
        this->quartic_equation(coeffs_[0],coeffs_[1],coeffs_[2],coeffs_[3],coeffs_[4],vect_pxv_);
        nSol_=vect_pxv_[0];
        
            for(int i=1;i<=nSol_;++i)
            {
                this->topRec(vect_pxv_[i]);
                TopSolution TS_temp;
                
                TS_temp.top = top_;
                TS_temp.topbar = topbar_;
                TS_temp.wp = w_;
                TS_temp.wm = wbar_;
                TS_temp.neutrino = neutrino_;
                TS_temp.neutrinobar = neutrinobar_;
                
                //proton Energy [GeV]  
                double protonE = 4000; //FIXME:  set as global variable
                TS_temp.x1 = (top_.E()+topbar_.E()+top_.Pz()+topbar_.Pz())/(2*protonE);
                TS_temp.x2 = (top_.E()+topbar_.E()-top_.Pz()-topbar_.Pz())/(2*protonE);
                TS_temp.mtt = tt_.M();
                
                TS_temp.weight = 1.0/tt_.M();
                //TS_temp.weight=Landau2D(neutrino_.E(),neutrinobar_.E());
                
                ttSol_.push_back(TS_temp);

            }
        nSol_=ttSol_.size();
        if(nSol_>0) this->sortTopSol(ttSol_);
}



void KinematicReconstruction_LSroutines::filldR()
{
    for(int i=0;i<nSol_;++i)
            {
                ttSol_[i].dR = sqrt(pow(ttSol_[i].top.DeltaR(true_top_),2)+pow(ttSol_[i].topbar.DeltaR(true_topbar_),2));   
            }
}



void KinematicReconstruction_LSroutines::filldN()
{
    for(int i=0;i<nSol_;++i)
       {
           ttSol_[i].dN = sqrt(pow((ttSol_[i].neutrino.Px()-true_neutrino_.Px()),2)+pow((ttSol_[i].neutrino.Py()-true_neutrino_.Py()),2)+pow((ttSol_[i].neutrino.Pz()-true_neutrino_.Pz()),2)+pow((ttSol_[i].neutrinobar.Px()-true_neutrinobar_.Px()),2)+pow((ttSol_[i].neutrinobar.Py()-true_neutrinobar_.Py()),2)+pow((ttSol_[i].neutrinobar.Pz()-true_neutrinobar_.Pz()),2));
       }
}



void KinematicReconstruction_LSroutines::swapTopSol(KinematicReconstruction_LSroutines::TopSolution& sol1, KinematicReconstruction_LSroutines::TopSolution& sol2)const
{
    KinematicReconstruction_LSroutines::TopSolution aux = sol1;
    sol1 = sol2;
    sol2 = aux;
}



void KinematicReconstruction_LSroutines::sortBy(std::string ch)
{
    if(ch=="dR"&&ttSol_.size()>0)
    {
     for(uint i=0;i<ttSol_.size()-1;++i)
        {
            if(ttSol_[i].dR > ttSol_[i+1].dR){ swapTopSol(ttSol_[i],ttSol_[i+1]);i=-1;}
        }   
    }
    if(ch=="dN"&&ttSol_.size()>0)
    {
     for(uint i=0;i<ttSol_.size()-1;++i)
        {
            if(ttSol_[i].dN > ttSol_[i+1].dN){ swapTopSol(ttSol_[i],ttSol_[i+1]);i=-1; }
        }   
    }
    if(ch=="dRN"&&ttSol_.size()>0)
    {
     for(uint i=0;i<ttSol_.size()-1;++i)
        {
            if(ttSol_[i].dN*ttSol_[i].dR > ttSol_[i+1].dN*ttSol_[i+1].dR){ swapTopSol(ttSol_[i],ttSol_[i+1]);i=-1; }
        }
    }
}



void KinematicReconstruction_LSroutines::sortTopSol(std::vector<KinematicReconstruction_LSroutines::TopSolution>& v)const
{
    //std::vector< KinematicReconstruction_LSroutines::TopSolution > result;
    for(uint i=0;i<v.size()-1;++i)
    {
      if(v[i].weight < v[i+1].weight){this->swapTopSol(v[i],v[i+1]);i=-1;}
    }
    
    //v.swap(result);
}



void KinematicReconstruction_LSroutines::topRec(const double& px_neutrino)
{
    double pxp, pyp, pzp, pup, pvp, pwp;
    
    d0_ = d00_;
    d1_ = d11_+d10_*px_neutrino;
    d2_ = d22_+d21_*px_neutrino+d20_*px_neutrino*px_neutrino;
    
    c0_ = c00_;
    c1_ = c11_+c10_*px_neutrino;
    c2_ = c22_+c21_*px_neutrino+c20_*px_neutrino*px_neutrino;
    
    
    pup = px_neutrino;
    pvp = (c0_*d2_-c2_*d0_)/(c1_*d0_-c0_*d1_);
    pwp = (-1)*(a1_+a2_*pup+a3_*pvp)/a4_;
    
    pxp = px_miss_-pup;   
    pyp = py_miss_-pvp;
    pzp = (-1)*(b1_+b2_*pxp+b3_*pyp)/b4_;
    
    neutrinobar_.SetXYZM(pxp, pyp, pzp, mav_);
    neutrino_.SetXYZM(pup, pvp, pwp, mv_);
        
    top_ = b_ + al_ + neutrino_;
    topbar_ = bbar_ + l_ + neutrinobar_; 
    tt_=top_+topbar_;
    w_ = al_ + neutrino_;
    wbar_ = l_ + neutrinobar_;
}



void KinematicReconstruction_LSroutines::findCoeff(double* const koeficienty)
{
    a1_ = ((b_.E()+al_.E())*(mw_*mw_-mal_*mal_-mv_*mv_)-al_.E()*(mt_*mt_-mb_*mb_-mal_*mal_-mv_*mv_)+2*b_.E()*al_.E()*al_.E()-2*al_.E()*(al_.Vect().Dot(b_.Vect())))/(2*al_.E()*(b_.E()+al_.E()));
    a2_ = 2*(b_.E()*al_.Px()-al_.E()*b_.Px())/(2*al_.E()*(b_.E()+al_.E()));
    a3_ = 2*(b_.E()*al_.Py()-al_.E()*b_.Py())/(2*al_.E()*(b_.E()+al_.E()));
    a4_ = 2*(b_.E()*al_.Pz()-al_.E()*b_.Pz())/(2*al_.E()*(b_.E()+al_.E()));
        
    //printf("Koefs ai: %f %f %f %f\n",a1_,a2_,a3_,a4_);//printout
    
    b1_ = ((bbar_.E()+l_.E())*(mwbar_*mwbar_-ml_*ml_-mav_*mav_)-l_.E()*(mtbar_*mtbar_-mbbar_*mbbar_-ml_*ml_-mav_*mav_)+2*bbar_.E()*l_.E()*l_.E()-2*l_.E()*(l_.Vect().Dot(bbar_.Vect())))/(2*l_.E()*(bbar_.E()+l_.E()));
    b2_ = 2*(bbar_.E()*l_.Px()-l_.E()*bbar_.Px())/(2*l_.E()*(bbar_.E()+l_.E()));
    b3_ = 2*(bbar_.E()*l_.Py()-l_.E()*bbar_.Py())/(2*l_.E()*(bbar_.E()+l_.E()));
    b4_ = 2*(bbar_.E()*l_.Pz()-l_.E()*bbar_.Pz())/(2*l_.E()*(bbar_.E()+l_.E()));
    
    //printf("Koefs bi: %f %f %f %f\n",b1_,b2_,b3_,b4_);//printout
    
    c22_ = (sqr((mw_*mw_-mal_*mal_-mv_*mv_))-4*(sqr(al_.E())-sqr(al_.Pz()))*sqr(a1_/a4_)-4*(mw_*mw_-mal_*mal_-mv_*mv_)*al_.Pz()*(a1_/a4_))/sqr(2*(b_.E()+al_.E())); 
    c21_ = (4*(mw_*mw_-mal_*mal_-mv_*mv_)*(al_.Px()-al_.Pz()*(a2_/a4_))-8*(sqr(al_.E())-sqr(al_.Pz()))*(a1_*a2_/sqr(a4_))-8*al_.Px()*al_.Pz()*(a1_/a4_))/sqr(2*(b_.E()+al_.E())); 
    c20_ = (-4*(sqr(al_.E())-sqr(al_.Px()))-4*(sqr(al_.E())-sqr(al_.Pz()))*sqr(a2_/a4_)-8*al_.Px()*al_.Pz()*(a2_/a4_))/sqr(2*(b_.E()+al_.E())); 
    c11_ = (4*(mw_*mw_-mal_*mal_-mv_*mv_)*(al_.Py()-al_.Pz()*(a3_/a4_))-8*(sqr(al_.E())-sqr(al_.Pz()))*(a1_*a3_/sqr(a4_))-8*al_.Py()*al_.Pz()*(a1_/a4_))/sqr(2*(b_.E()+al_.E())); 
    c10_ = (-8*(sqr(al_.E())-sqr(al_.Pz()))*(a2_*a3_/sqr(a4_)) + 8*al_.Px()*al_.Py() - 8*al_.Px()*al_.Pz()*(a3_/a4_) - 8*al_.Py()*al_.Pz()*(a2_/a4_))/sqr(2*(b_.E()+al_.E()));
    c00_ = (-4*(sqr(al_.E())-sqr(al_.Py())) -4*(sqr(al_.E())-sqr(al_.Pz()))*sqr(a3_/a4_)-8*al_.Py()*al_.Pz()*(a3_/a4_))/sqr(2*(b_.E()+al_.E()));
    
    // printf("Koefs ci: %f %f %f %f %f %f\n",c22_,c21_/c22_,c11_/c22_,c20_/c22_,c10_/c22_,c00_/c22_);//printout
    
    
    
    double D22,D21,D20,D11,D10,D00;
    D22 = (sqr((mwbar_*mwbar_-ml_*ml_-mav_*mav_))-4*(sqr(l_.E())-sqr(l_.Pz()))*sqr(b1_/b4_)-4*(mwbar_*mwbar_-ml_*ml_-mav_*mav_)*l_.Pz()*(b1_/b4_))/sqr(2*(bbar_.E()+l_.E())); 
    D21 = (4*(mwbar_*mwbar_-ml_*ml_-mav_*mav_)*(l_.Px()-l_.Pz()*(b2_/b4_))-8*(sqr(l_.E())-sqr(l_.Pz()))*(b1_*b2_/sqr(b4_))-8*l_.Px()*l_.Pz()*(b1_/b4_))/sqr(2*(bbar_.E()+l_.E())); 
    D20 = (-4*(sqr(l_.E())-sqr(l_.Px()))-4*(sqr(l_.E())-sqr(l_.Pz()))*sqr(b2_/b4_)-8*l_.Px()*l_.Pz()*(b2_/b4_))/sqr(2*(bbar_.E()+l_.E())); 
    D11 = (4*(mwbar_*mwbar_-ml_*ml_-mav_*mav_)*(l_.Py()-l_.Pz()*(b3_/b4_))-8*(sqr(l_.E())-sqr(l_.Pz()))*(b1_*b3_/sqr(b4_))-8*l_.Py()*l_.Pz()*(b1_/b4_))/sqr(2*(bbar_.E()+l_.E())); 
    D10 = (-8*(sqr(l_.E())-sqr(l_.Pz()))*(b2_*b3_/sqr(b4_)) + 8*l_.Px()*l_.Py() - 8*l_.Px()*l_.Pz()*(b3_/b4_) - 8*l_.Py()*l_.Pz()*(b2_/b4_))/sqr(2*(bbar_.E()+l_.E()));
    D00  = (-4*(sqr(l_.E())-sqr(l_.Py())) -4*(sqr(l_.E())-sqr(l_.Pz()))*sqr(b3_/b4_)-8*l_.Py()*l_.Pz()*(b3_/b4_))/sqr(2*(bbar_.E()+l_.E()));
    
    //printf("Koefs di_: %f %f %f %f %f %f\n",D22,D21/D22,D11/D22,D20/D22,D10/D22,D00/D22);//printout
    
    
    d22_ = D22+sqr(px_miss_)*D20+sqr(py_miss_)*D00+px_miss_*py_miss_*D10+px_miss_*D21+py_miss_*D11;
    d21_ = -D21-2*px_miss_*D20-py_miss_*D10;
    d20_ = D20;
    d11_ = -D11-2*py_miss_*D00-px_miss_*D10;
    d10_ = D10;
    d00_  = D00;
    
    //printf("Koefs di: %f %f %f %f %f %f\n",d22_, d21_/d22_, d11_/d22_, d20_/d22_, d10_/d22_, d00_/d22_);//printout
    
    
    koeficienty[4] = sqr(c00_)*sqr(d22_)+c11_*d22_*(c11_*d00_-c00_*d11_)+c00_*c22_*(sqr(d11_)-2*d00_*d22_)+c22_*d00_*(c22_*d00_-c11_*d11_);
    koeficienty[3] = c00_*d21_*(2*c00_*d22_-c11_*d11_)+c00_*d11_*(2*c22_*d10_+c21_*d11_)+c22_*d00_*(2*c21_*d00_-c11_*d10_)-c00_*d22_*(c11_*d10_+c10_*d11_)-2*c00_*d00_*(c22_*d21_+c21_*d22_)-d00_*d11_*(c11_*c21_+c10_*c22_)+c11_*d00_*(c11_*d21_+2*c10_*d22_);
    koeficienty[2] = sqr(c00_)*(2*d22_*d20_+sqr(d21_))-c00_*d21_*(c11_*d10_+c10_*d11_)+c11_*d20_*(c11_*d00_-c00_*d11_)+c00_*d10_*(c22_*d10_-c10_*d22_)+c00_*d11_*(2*c21_*d10_+c20_*d11_)+(2*c22_*c20_+sqr(c21_))*sqr(d00_)-2*c00_*d00_*(c22_*d20_+c21_*d21_+c20_*d22_)+c10_*d00_*(2*c11_*d21_+c10_*d22_)-d00_*d10_*(c11_*c21_+c10_*c22_)-d00_*d11_*(c11_*c20_+c10_*c21_);
    // koeficienty[1] = c00_*d21_*(2*c00_*d20_-c10_*d10_)-c00_*d20_*(c11_*d10_+c10_*d11_)+c00_*d10_*(c21_*d10_+2*c20_*d11_)-2*c00_*d00_*(c21_*d20_+c20_*d21_)+c10_*d00_*(2*c11_*d20_+c10_*d21_)-c20_*d00_*(2*c21_*d00_-c10_*d11_)-d00_*d10_*(c11_*c20_+c10_*c21_);
    koeficienty[1] = c00_*d21_*(2*c00_*d20_-c10_*d10_)-c00_*d20_*(c11_*d10_+c10_*d11_)+c00_*d10_*(c21_*d10_+2*c20_*d11_)-2*c00_*d00_*(c21_*d20_+c20_*d21_)+c10_*d00_*(2*c11_*d20_+c10_*d21_)+c20_*d00_*(2*c21_*d00_-c10_*d11_)-d00_*d10_*(c11_*c20_+c10_*c21_);
    koeficienty[0] = sqr(c00_)*sqr(d20_)+c10_*d20_*(c10_*d00_-c00_*d10_)+c20_*d10_*(c00_*d10_-c10_*d00_)+c20_*d00_*(c20_*d00_-2*c00_*d20_);
    //printf("Koefs_in_f: %15.15f %f %f %f %f\n",koeficienty[0],koeficienty[1],koeficienty[2],koeficienty[3],koeficienty[4]); //printout
}



double KinematicReconstruction_LSroutines::sqr(const double& x)const
{
    return (x*x);
}



void KinematicReconstruction_LSroutines::swap(double& realone, double& realtwo)const
{
    if(realtwo < realone){
        double aux = realtwo;
        realtwo = realone;
        realone = aux;
    }
}



int KinematicReconstruction_LSroutines::sign(const long double& ld)const
{
    if(std::abs(ld)<0.0000000000001) return 0;
    return (ld>0)?1:-1;
}



void KinematicReconstruction_LSroutines::quartic_equation(const double& h0, const double& h1, const double& h2, const double& h3, const double& h4, std::vector<double>& v)const
{
     std::vector<double> result;
    
    //printf("Koefs_in_f: %f %f %f %f %f\n",h0,h1,h2,h3,h4); //printout
    //printf("Koefs_norm_in_f: %f %f %f %f %f\n",h0/h0,h1/h0,h2/h0,h3/h0,h4/h0); //printout


        if(sign(a4_)==0||sign(b4_)==0)
        {
            result.push_back(0);
            v.swap(result);    
        }
        else
        {
               //printf("else1\n"); //printout
            if(sign(h0)==0)
            {
                this->cubic_equation(h1,h2,h3,h4,result);
                v.swap(result);
            }
           else
            {
                //printf("else2\n"); //printout
                if(sign(h4)==0)
                {
                    this->cubic_equation(h0,h1,h2,h3,result);
                    result[0]=result[0]+1;
                    result.push_back(0);
                    v.swap(result);
                }
                else
                {
                   //printf("else3\n"); //printout
                    
                     double H1=h1/h0;
                     double H2=h2/h0;
                     double H3=h3/h0;
                     double H4=h4/h0;
                   double K1 = H2 -3*sqr(H1)/8;
                   double K2 = H3 + H1*sqr(H1)/8-H1*H2/2;
                   double K3 = H4-3*sqr(sqr(H1))/256+sqr(H1)*H2/16-H1*H3/4;
                    //printf("Koefs Ki: %f %f %10.10f\n",K1,K2,K3);//printout
                    if(sign(K3)==0)
                    {
                       this->cubic_equation(1,0,K1,K2,result);
                       for(int i=1;i<=result[0];++i)
                       {
                           result[i]=result[i]-H1/4;
                       }
                        result[0]=result[0]+1;
                        result.push_back(-H1/4);
                        v.swap(result);
                       
                    }
                    else
                    {
                        //printf("else4\n"); //printout
                        std::vector<double> result_t12;
                        
                        std::vector<double> result_t1;
                            result_t1.push_back(0);
                        std::vector<double> result_t2;
                            result_t2.push_back(0);
                        
                        this->cubic_equation(1,2*K1,(K1*K1-4*K3),(-1)*K2*K2,result_t12); 
                        
                        //std::cout << "hehehe:  " << result_t12[0]  <<  std::endl; //printout
                        
                        
                        for(int i=1;i<=result_t12[0];++i)
                        {
                            //std::cout << "heh:  " << result_t12[i]  << std::endl; //printout
                            
                            if(result_t12[i]>=0)
                            {
                                result_t1[0]=result_t1[0]+2;
                                result_t1.push_back(sqrt(result_t12[i]));
                                result_t1.push_back((-1)*sqrt(result_t12[i]));
                                result_t2[0]=result_t2[0]+2;
                                result_t2.push_back((K1+result_t12[i]-K2/sqrt(result_t12[i]))/2);
                                result_t2.push_back((K1+result_t12[i]+K2/sqrt(result_t12[i]))/2);                                
                            }
                        }  
                        
                        //std::cout  << std::endl;
                        
                        
                        std::vector<double> pre_result1;

                        result.push_back(0);
                        for(int i=1;i<=result_t1[0];++i)
                        {
                            //std::cout << "quadric_equation:   " << i << " " << result_t1[i] << " " << result_t2[i] << std::endl; //printout
                             
                             this->quadratic_equation(1,result_t1[i],result_t2[i],pre_result1);
                             
                             for(int j=1;j<=pre_result1[0];++j)
                             {
                                // if(pre_result1[0]==2)std::cout << "quadric_equation:   " << i << " " << pre_result1[1] << " " << pre_result1[2] << std::endl; //printout

                                 int flag=1;
                                for(int r=1;r<=result.at(0);++r)
                                {
                                 if(fabs(result.at(r)-pre_result1[j])<0.02)flag=0;
                                 //printf("Result-result: %10.10f  \n",result[r]-pre_result1[j]);
                                }
                                if(flag)
                                {
                                    result.at(0)=result.at(0)+1;
                                    result.push_back(pre_result1[j]);
                                }
                             }
                            pre_result1.clear();                          
                        }                          
                       for(int k=1;k<=result.at(0);++k)
                       {
                           
                           //printf("Result: %f   %f \n",H1/4,h1/4); //printout
                           result.at(k)=result.at(k)-H1/4;
                           
                       }
                       v.swap(result);   
                    }
                }
            }
        }
}



void KinematicReconstruction_LSroutines::cubic_equation(const double& a, const double& b, const double& c, const double& d, std::vector<double>& v)const
{
        
    std::vector<double> result;
    if(a==0)
    {
        this->quadratic_equation(b,c,d,result);
        v.swap(result);
    }
    else
    {
       double s1 = b/a;
       double s2 = c/a;
       double s3 = d/a;
       
       double q = (s1*s1-3*s2)/9;
       double q3 = q*q*q;
       double r = (2*s1*s1*s1-9*s1*s2+27*s3)/54;
       double r2 = r*r;
       double S = r2-q3;
       if(sign(S)<0)
       {
           
           result.push_back(3);
           double F = acos(r/sqrt(q3));
           result.push_back(-2*sqrt(fabs(q))*cos(F/3)-s1/3);
           result.push_back(-2*sqrt(fabs(q))*cos((F+2*TMath::Pi())/3)-s1/3);
           result.push_back(-2*sqrt(fabs(q))*cos((F-2*TMath::Pi())/3)-s1/3);  
           v.swap(result);
           
       }
       else 
       {
           if(sign(S)==0)
           {
                long double A = r+sqrt(fabs(r2-q3));
                A = A<0 ? pow(fabs(A),(long double)1.0/3) : -pow(fabs(A),(long double)1.0/3);
                long double B = sign(A) == 0 ? 0 : q/A; 
                result.push_back(2);
                result.push_back(A+B-s1/3);
                result.push_back(-0.5*(A+B)-s1/3);  //!!!
                v.swap(result);
           }
           else
           {
               long double A = r+sqrt(fabs(r2-q3));
               A = A<0 ? pow(fabs(A),(long double)1.0/3) : -pow(fabs(A),(long double)1.0/3);
               long double B = sign(A) == 0 ? 0 : q/A; 
               result.push_back(1);
               result.push_back(A+B-s1/3);
               v.swap(result);
           }
       }
       
    }
}



void KinematicReconstruction_LSroutines::quadratic_equation(const double& a, const double& b, const double& c, std::vector<double>& v)const
{
     std::vector<double> result;
     //printf("a: %10.10f\n",a);//printout
    if(a==0)
    {
        this->linear_equation(b,c,result);
        v.swap(result);
    }
    else
    {
        double D = b*b-4*a*c;
        //printf("D: %10.10f\n",D);//printout
        if(this->sign(D)<0)
        {
            result.push_back(0);
            v.swap(result);
        }
        else 
        {
            if(sign(D)==0)
            {
                result.push_back(1);
                result.push_back((-1)*b/(2*a));
                v.swap(result);
            }
            else
            {
                result.push_back(2);
                result.push_back((-b-sqrt(D))/(2*a));
                result.push_back((-b+sqrt(D))/(2*a));
                v.swap(result);
            }
        }
    }
}



void KinematicReconstruction_LSroutines::linear_equation(const double& a, const double& b, std::vector<double>& v)const
{
    std::vector<double> result;
    if(a==0)
    {
        result.push_back(0);
        v.swap(result);
    }
    else
    {
        result.push_back(1);
        result.push_back((-1)*(b/a));
        v.swap(result);
    }
}


