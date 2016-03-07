#ifndef classes_lorentz_for_top_Fwd_h
#define classes_lorentz_for_top_Fwd_h

#include <vector>





// For LorentzVector - ROOT dictionary

namespace ROOT{
    namespace Math{
        template<typename T> class PtEtaPhiM4D;
        template<class T> class LorentzVector;
    }
}

/// Our Lorentz vector as used in the nTuple
typedef ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiM4D<double> > LV;

/// Vector of our Lorentz vector as used in the nTuple
typedef std::vector<LV> VLV;



#endif



