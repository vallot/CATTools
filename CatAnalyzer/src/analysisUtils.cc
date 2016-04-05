#include <iostream>
#include <cstdlib>
#include <cmath>
#include <cstdio>

#include <TLorentzVector.h>

#include "CATTools/CatAnalyzer/interface/analysisUtils.h"
#include "CATTools/CatAnalyzer/interface/classes.h"






// --- Several conversion functions -------------------------------------------------------------------------------------

void common::LVtod4(const LV& lv, double* d)
{
    d[0] = lv.E();
    d[1] = lv.Px();
    d[2] = lv.Py();
    d[3] = lv.Pz();
}



std::string common::d2s(const double& d)
{
    char result[100];
    if(std::abs(d) < 5) {
        std::sprintf(result, "%.3f", d);
        std::string s = std::string(result);
        while (s.length() > 0 && s[s.length()-1] == '0') s.erase(s.end()-1);
        if (s.length() > 0 && s[s.length()-1] == '.') s.erase(s.end()-1);
        return s;
    }
    else {
        std::sprintf(result, "%.0f", d);
        return std::string(result);
    }
}



const TLorentzVector common::LVtoTLV(const LV& lv)
{
    return TLorentzVector(lv.X(), lv.Y(), lv.Z(), lv.T());
}



const LV common::TLVtoLV(const TLorentzVector& lv)
{
    LV result; 
    result.SetXYZT(lv.X(), lv.Y(), lv.Z(), lv.T());
    return result;
}







// --- Functions concerning the treatment of indices of vectors (for working with data stored in nTuple branches) -------------

std::vector<int> common::mergeIndices(const std::vector<int>& v_index1, const std::vector<int>& v_index2, const bool allowOverlap)
{
    std::vector<int> result(v_index1);
    
    for(const int index : v_index2){
        // Check if index is already contained in first vector
        if(std::find(v_index1.begin(), v_index1.end(), index) != v_index1.end()){
            if(allowOverlap) continue;
            else{
                std::cerr<<"ERROR in common::mergeIndices()! Two collections should be merged, but they overlap (not allowed)\n...break\n"<<std::endl;
                exit(92);
            }
        }
        else result.push_back(index);
    }
    
    return result;
}



std::vector<double> common::parametersLV(const VLV& v_lv, const common::LVParameter& parameter)
{
    std::vector<double> v_variable;
    for(const LV& lv : v_lv){
        if(parameter == LVpt) v_variable.push_back(lv.pt());
        else if (parameter == LVeta) v_variable.push_back(lv.eta());
        else{
            std::cerr<<"Error in parametersLV! Lorentz vector parameter is not implemented\n...break\n";
            exit(638);
        }
    }
    return v_variable;
}



void common::orderIndices(int& index1, int& index2, const VLV& v_lv, const common::LVParameter& parameter, const bool absoluteValue)
{
    const std::vector<double> v_variable = parametersLV(v_lv, parameter);
    orderIndices(index1, index2, v_variable, absoluteValue);
}



void common::orderIndices(std::vector<int>& v_index, const VLV& v_lv, const common::LVParameter& parameter, const bool absoluteValue)
{
    const std::vector<double> v_variable = parametersLV(v_lv, parameter);
    orderIndices(v_index, v_variable, absoluteValue);
}



void common::selectIndicesLeptonIso(std::vector<int>& v_index, const VLV& v_lvLep, const VLV& v_lvJet)
{
    std::vector<int> result;
    for(const int index : v_index){
        bool isReject = 0;
        for(const LV& jet : v_lvJet){
            TLorentzVector jetTLV = common::LVtoTLV(jet);
            TLorentzVector lepTLV = common::LVtoTLV(v_lvLep.at(index));
            double dR = jetTLV.DeltaR(lepTLV);
            double relPt = (lepTLV.Vect().Cross(jetTLV.Vect()).Mag())/(jetTLV.Vect().Mag());
            if(dR < 0.5 && relPt < 15)isReject=1;
        }
        if(!isReject)result.push_back(index);
    }
    
    v_index.clear();
    v_index = result;
}



void common::selectIndices(std::vector<int>& v_index, const VLV& v_lv, const common::LVParameter& parameter,
                          const double cutValue, const bool lowerThreshold)
{
    const std::vector<double> v_variable = parametersLV(v_lv, parameter);
    selectIndices(v_index, v_variable, cutValue, lowerThreshold);
}

void common::selectIndices(std::vector<int>& v_index, const std::vector<int>& v_variable0, const std::vector<double>& v_variable,
                          const double cutValue, const bool lowerThreshold)
{
    std::vector<int> result;
    for(const int index : v_index){
        const int& variable0(v_variable0.at(index));
        if(std::abs(variable0)==11) {
        const double& variable(v_variable.at(index));
            if(lowerThreshold){
                if(variable<cutValue) continue;
            }
            else{
                if(variable>cutValue) continue;
            }
        }
        result.push_back(index);
    }
    v_index.clear();
    v_index = result;

}


int common::extremumIndex(const VLV& v_lv, const common::LVParameter& parameter, const bool maximumValue)
{
    const std::vector<double> v_variable = parametersLV(v_lv, parameter);
    return extremumIndex(v_variable, maximumValue);
}



void common::orderLV(LV& lv1, LV& lv2, const LV& inputLv1, const LV& inputLv2, const common::LVParameter& parameter, const bool absoluteValue)
{
    double variable1;
    double variable2;
    if(parameter == LVpt){
        variable1 = inputLv1.pt();
        variable2 = inputLv2.pt();
    }
    else if(parameter == LVeta){
        variable1 = inputLv1.eta();
        variable2 = inputLv2.eta();
    }
    else{
        std::cerr<<"Error in orderLV! Lorentz vector parameter is not implemented\n...break\n";
        exit(639);
    }
    
    if(absoluteValue){
        variable1 = std::abs(variable1);
        variable2 = std::abs(variable2);
    }
    
    if (variable1 > variable2) {
        lv1 = inputLv1;
        lv2 = inputLv2;
    } else {
        lv1 = inputLv2;
        lv2 = inputLv1;
    }
}









