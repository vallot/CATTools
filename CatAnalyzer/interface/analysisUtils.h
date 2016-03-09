#ifndef common_analysisUtils_h
#define common_analysisUtils_h

#include <vector>
#include <cmath>
#include <string>

class TLorentzVector;

#include "classesFwd.h"





// --- Several conversion functions -------------------------------------------------------------------------------------

namespace common{
    /// Convert LorentzVector to an array[E, px, py, pz]
    void LVtod4(const LV& lv, double* d);
    
    /// Convert double to string (smart number of digits)
    std::string d2s(const double& d);
    
    
    /// Conversion from TLorentzVector to our LV type (ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiM4D<double> >)
    const LV TLVtoLV(const TLorentzVector& lv);
    
    
    /// Conversion from our LV type to TLorentzVector (ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiM4D<double> >)
    const TLorentzVector LVtoTLV(const LV& lv);
}





// --- Functions concerning the treatment of indices of vectors (for working with data stored in nTuple branches) -------------

namespace common{
    
    /// Setting up a vector with all indices corresonding to its size N {0, 1, 2, ..., N-1}
    template<class T> std::vector<int> initialiseIndices(const std::vector<T>& v_variable);
    
    
    
    /// Template function to order two indices by comparison operator for the two elements of a vector
    /// corresponding to the indices of a given variable
    /// Result is (if absoluteValue==false): variable.at(index1) > variable.at(index2)
    /// Result is (if absoluteValue==true): |variable.at(index1)| > |variable.at(index2)|
    template<class T> void orderIndices(int& index1, int& index2, const std::vector<T>& v_variable, const bool absoluteValue =false);
    
    /// Template function to order a vector of indices by comparison operator for the two elements of a vector
    /// corresponding to the indices of a given variable
    /// Result is (if absoluteValue==false): variable.at(index1) > variable.at(index2) > ... > variable.at(indexN)
    /// Result is (if absoluteValue==true): |variable.at(index1)| > |variable.at(index2)| > ... > |variable.at(indexN)|
    template<class T> void orderIndices(std::vector<int>& v_index, const std::vector<T>& v_variable, const bool absoluteValue =false);
    
    /// Template function to select from a vector of indices those whose corresponding value in the vector
    /// survive the cut, i.e. variable.at(index) > cutValue (if lowerThreshold==true) or
    /// variable.at(index) < cutValue (if lowerThreshold==false)
    template<class T> void selectIndices(std::vector<int>& v_index, const std::vector<T>& v_variable,
                                         const T& cutValue, const bool lowerThreshold =true);
    
    
    
    /// Template function returning the index of a vector whose element has the most extreme value (in case of several, will be the first one found)
    /// Result is (if maximumValue==true): Index of element with maximum value
    /// Result is (if maximumValue==false): Index of element with minimum value
    /// If the vector is empty, return -1
    template <class T> int extremumIndex(const std::vector<T>& v_variable, const bool maximumValue =true);
    
    
    
    /// Merge the indices of two collections in one
    /// If allowOverlap=false, it is checked that there is no duplicated index, else an error is thrown
    /// If allowOverlap=true, both collections can have identical indices, the merged collection contains it then once
    std::vector<int> mergeIndices(const std::vector<int>& v_index1, const std::vector<int>& v_index2, const bool allowOverlap =false);
    
    
    
    /// Possible parameters of a lorentz vector to be used in following tools
    enum LVParameter{LVpt, LVeta};
    
    /// Extract a specific value for each entry of a vector of LVs
    std::vector<double> parametersLV(const VLV& v_lv, const LVParameter& parameter);
    
    /// Function to order two indices by comparison operator for the parameter of the LVs of a vector
    /// corresponding to the indices of a given variable
    /// Result is (if absoluteValue==false): LV.at(index1).parameter > LV.at(index2).parameter
    /// Result is (if absoluteValue==true): |LV.at(index1).parameter| > |LV.at(index2).parameter|
    void orderIndices(int& index1, int& index2, const VLV& v_lv, const LVParameter& parameter, const bool absoluteValue =false);
    
    /// Function to order a vector of indices by comparison operator for the parameter of the lorentz vectors of a vector
    /// corresponding to the indices of a given variable
    /// Result is (if absoluteValue==false): LV.at(index1).parameter > LV.at(index2).parameter > ... > LV.at(indexN).parameter
    /// Result is (if absoluteValue==true): |LV.at(index1).parameter| > |LV.at(index2).parameter| > ... > |LV.at(indexN).parameter|
    void orderIndices(std::vector<int>& v_index, const VLV& v_lv, const LVParameter& parameter, const bool absoluteValue =false);
    
    /// Function to select from a vector of indices those whose corresponding value of a LV parameter
    /// in the vector survive the cut
    void selectIndices(std::vector<int>& v_index, const VLV& v_lv, const LVParameter& parameter,
                                         const double cutValue, const bool lowerThreshold =true);
    
    /// Function to select from a vector of indices those whose corresponding value of a LV parameter
    /// in the vector survive the cut REMOVING ELECTRONS WITH D0>0.04
    void selectIndices(std::vector<int>& v_index, const std::vector<int>& v_variable0, const std::vector<double>& v_variable, 
                                         const double cutValue, const bool lowerThreshold =true);
    
    /// Function to select from a vector of lepton indices those whose fulfill iso-selection criteria
    void selectIndicesLeptonIso(std::vector<int>& v_index, const VLV& v_lvLep, const VLV& v_lvJet);
    
    /// Function returning the index of a vector whose parameter of the LVs has the most extreme value (in case of several, will be the first one found)
    /// Result is (if maximumValue==true): Index of element with maximum value of LV parameter
    /// Result is (if maximumValue==false): Index of element with minimum value of LV parameter
    /// If the vector is empty, return -1
    int extremumIndex(const VLV& v_lv, const LVParameter& parameter, const bool maximumValue =true);
    
    
    
    /// Function to order two Lorentz vectors by comparison operator for the parameter of the LVs
    /// Result is (if absoluteValue==false): lv1.parameter > lv2.parameter
    /// Result is (if absoluteValue==true): |lv1.parameter| > |lv2.parameter|
    void orderLV(LV& lv1, LV& lv2, const LV& inputLv1, const LV& inputLv2, const LVParameter& parameter, const bool absoluteValue =false);
}










template<class T> std::vector<int> common::initialiseIndices(const std::vector<T>& v_variable)
{
    std::vector<int> result;
    for(size_t i=0; i<v_variable.size(); ++i) result.push_back(i);
    return result;
}



template<class T> void common::orderIndices(int& index1, int& index2, const std::vector<T>& v_variable, const bool absoluteValue)
{
    const T& variable1 = absoluteValue ? std::abs(v_variable.at(index1)) : v_variable.at(index1);
    const T& variable2 = absoluteValue ? std::abs(v_variable.at(index2)) : v_variable.at(index2);
    if(variable1<variable2){
        const int tmpIndex2 = index1;
        index1 = index2;
        index2 = tmpIndex2;
    }
}



template<class T> void common::orderIndices(std::vector<int>& v_index, const std::vector<T>& v_variable, const bool absoluteValue)
{
    std::vector<int> v_indexResult;
    for(std::vector<int>::const_iterator i_index2 = v_index.begin(); i_index2 != v_index.end(); ++i_index2){
        const int index2(*i_index2);
        if(i_index2 == v_index.begin()){
            v_indexResult.push_back(index2);
            continue;
        }
        const T& variable2 = absoluteValue ? std::abs(v_variable.at(index2)) : v_variable.at(index2);
        
        bool isInserted(false);
        for(std::vector<int>::iterator i_index1 = v_indexResult.begin(); i_index1 != v_indexResult.end(); ++i_index1){
            const T& variable1 = absoluteValue ? std::abs(v_variable.at(*i_index1)) : v_variable.at(*i_index1);
            if(variable1<variable2){
                v_indexResult.insert(i_index1, index2);
                isInserted = true;
                break;
            }
        }
        if(!isInserted) v_indexResult.push_back(index2);
    }
    v_index = v_indexResult;
}



template<class T> void common::selectIndices(std::vector<int>& v_index, const std::vector<T>& v_variable,
                                            const T& cutValue, const bool lowerThreshold)
{
    std::vector<int> result;
    for(const int index : v_index){
        const T& variable(v_variable.at(index));
        if(lowerThreshold){
            if(variable<cutValue) continue;
        }
        else{
            if(variable>cutValue) continue;
        }
        result.push_back(index);
    }
    v_index.clear();
    v_index = result;
}



template<class T> int common::extremumIndex(const std::vector<T>& v_variable, const bool maximumValue)
{
    if(v_variable.size() == 0) return -1;
    
    int result(0);
    T extremum(v_variable.at(0));
    
    for(size_t i = 1; i < v_variable.size(); ++i){
        if(maximumValue){
            if(extremum < v_variable.at(i)){
                extremum = v_variable.at(i);
                result = i;
            }
        }
        else{
            if(extremum > v_variable.at(i)){
                extremum = v_variable.at(i);
                result = i;
            }
        }
    }
    
    return result;
}







#endif







