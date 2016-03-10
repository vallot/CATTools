#ifndef common_utils_h
#define common_utils_h

#include <string>
#include <functional>
#include <vector>

class TString;




namespace common{

    /// return CMSSW_BASE environment variable as string, with error checking
    const std::string CMSSW_BASE();

    /// Return the path where relevant input data (e.g. histograms for scale factors) is stored
    const std::string DATA_PATH_COMMON();
    
    
    
    /**
     * Helper function to create a function which checks if a string is in the passed vector of string
     * (Check for exact agreement)
     * 
     * @param v_string a vector of allowed strings (TString)
     * @return a function taking a std::string and returning a bool
     */
    #ifndef __CINT__
    std::function<bool(const std::string& s)> makeStringCheck(const std::vector<TString> v_string);
    #endif
    
    /**
     * Helper function to create a function which checks if a string is in the passed vector of string
     * (Check only if the beginning of the string is contained in the vector of string)
     * 
     * @param v_string a vector of allowed strings (TString)
     * @return a function taking a std::string and returning a bool
     */
    #ifndef __CINT__
    std::function<bool(const std::string& s)> makeStringCheckBegin(const std::vector<TString> v_string);
    #endif
}







#endif




