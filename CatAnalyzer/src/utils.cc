#include <iostream>
#include <cstdlib>
#include <algorithm>

#include <TString.h>

#include "CATTools/CatAnalyzer/interface/utils.h"





const std::string common::CMSSW_BASE()
{
    const char* cmssw_base = std::getenv("CMSSW_BASE");
    if (!cmssw_base) {
        std::cerr << "Error! Environmental variable CMSSW_BASE not set!\n"
                  << "Please run cmsenv first.\n"
                  << "When running without CMSSW, you still need this variable so that the\n"
                  << "certain files can be found.\n";
        exit(1);            
    }
    std::string result(cmssw_base);
    return result;
}



const std::string common::DATA_PATH_COMMON()
{
    std::string result(CMSSW_BASE());
    result.append("/src/CATTools/CatAnalyzer/data");
    return result;
}



std::function<bool(const std::string& s)> common::makeStringCheck(const std::vector<TString> v_string)
{
    return [v_string](const std::string& test) -> bool {
        const TString tTest(test);
        return std::find(begin(v_string), end(v_string), tTest) != end(v_string);
    };
}



std::function<bool(const std::string& s)> common::makeStringCheckBegin(const std::vector<TString> v_string)
{
    return [v_string](const std::string& test) -> bool {
        const TString tTest(test);
        for(const auto& string : v_string){
            if(string == ""){
                if(tTest == "") return true;
                else return false;
            }
            else if(tTest.BeginsWith(string)) return true;
        }
        return false;
    };
}





