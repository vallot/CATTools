#ifndef LeptonWEIGHT_H
#define LeptonWEIGHT_H

#include <vector>
#include "TLorentzVector.h"
#include <iostream>
using namespace std;

class LeptonWeight
{

  public:
    LeptonWeight(){};
    ~LeptonWeight(){};

  enum LeptonType {
      Muon,
      Electron
  };


  float SF(float pt, float eta, LeptonType type){

    float scale = 1.0;
    if(type == Muon){
      scale = 1.0;
      /* //for 2012
      if( pt < 30 ){
        if( fabs(eta) >= 0.0 && fabs(eta) < 0.9 ) scale = 0.980;  
        if( fabs(eta) >= 0.9 && fabs(eta) < 1.2 ) scale = 0.972;  
        if( fabs(eta) >= 1.2 && fabs(eta) < 2.1 ) scale = 0.996;
        if( fabs(eta) >= 2.1 && fabs(eta) < 2.5 ) scale = 1.019;
      }
      if( pt >= 30 && pt < 40 ){
        if( fabs(eta) >= 0.0 && fabs(eta) < 0.9 ) scale = 0.988;
        if( fabs(eta) >= 0.9 && fabs(eta) < 1.2 ) scale = 0.994;
        if( fabs(eta) >= 1.2 && fabs(eta) < 2.1 ) scale = 0.982;
        if( fabs(eta) >= 2.1 && fabs(eta) < 2.5 ) scale = 1.018;
      }
      if( pt >= 40 && pt < 50 ){
        if( fabs(eta) >= 0.0 && fabs(eta) < 0.9 ) scale = 0.982;
        if( fabs(eta) >= 0.9 && fabs(eta) < 1.2 ) scale = 0.979;
        if( fabs(eta) >= 1.2 && fabs(eta) < 2.1 ) scale = 0.986;
        if( fabs(eta) >= 2.1 && fabs(eta) < 2.5 ) scale = 1.000;
      }
      if( pt >= 50 ){
        if( fabs(eta) >= 0.0 && fabs(eta) < 0.9 ) scale = 0.992;
        if( fabs(eta) >= 0.9 && fabs(eta) < 1.2 ) scale = 0.988;
        if( fabs(eta) >= 1.2 && fabs(eta) < 2.1 ) scale = 1.000;
        if( fabs(eta) >= 2.1 && fabs(eta) < 2.5 ) scale = 1.024;
      }*/
    }else if(type == Electron){
      /* // for 2012
      if( pt < 30 ){
        if( fabs(eta) >= 0.0 && fabs(eta) < 0.8 )   scale = 0.971;      
        if( fabs(eta) >= 0.8 && fabs(eta) < 1.479 ) scale = 0.962;      
        if( fabs(eta) >= 1.479 && fabs(eta) < 2.5 ) scale = 0.921;
      }
      if( pt >= 30 && pt < 40 ){
        if( fabs(eta) >= 0.0 && fabs(eta) < 0.8 )   scale = 0.942;
        if( fabs(eta) >= 0.8 && fabs(eta) < 1.479 ) scale = 0.956;
        if( fabs(eta) >= 1.479 && fabs(eta) < 2.5 ) scale = 0.921;
      }
      if( pt >= 40 && pt < 50 ){
        if( fabs(eta) >= 0.0 && fabs(eta) < 0.8 )   scale = 0.967;
        if( fabs(eta) >= 0.8 && fabs(eta) < 1.479 ) scale = 0.961;
        if( fabs(eta) >= 1.479 && fabs(eta) < 2.5 ) scale = 0.959;
      }
      if( pt >= 50 ){
        if( fabs(eta) >= 0.0 && fabs(eta) < 0.8 )   scale = 0.964;
        if( fabs(eta) >= 0.8 && fabs(eta) < 1.479 ) scale = 0.963;
        if( fabs(eta) >= 1.479 && fabs(eta) < 2.5 ) scale = 0.963;
      }*/
    }else{
      return scale = 1;
    }

    return scale;

  }

  
};

#endif

