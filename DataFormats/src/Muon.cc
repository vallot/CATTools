#include "CATTools/DataFormats/interface/Muon.h"

using namespace cat;

/// default constructor
Muon::Muon() {
}

Muon::Muon(const reco::LeafCandidate & aMuon) : Lepton( aMuon ) {
}

/// destructor
Muon::~Muon() {
}

float Muon::scaleFactor(const std::string& name) const {
  if (name == "NUM_TightIDandIPCut_DEN_genTracks_PAR_pt_spliteta_bin1"){
    float abseta = std::abs(this->eta());
    if (this->pt()>20. && this->pt() <= 25.){
      if      ( abseta > 0.0 && abseta <= 0.9) return 0.9752116203308105;
      else if ( abseta > 0.9 && abseta <= 1.2) return 0.9738101959228516;
      else if ( abseta > 1.2 && abseta <= 2.1) return 0.9983288645744324;
      else if ( abseta > 2.1 && abseta <= 2.4) return 0.9877836108207703;
      else return 1.;
    }    
    if (this->pt()>25. && this->pt() <= 30.){
      if      ( abseta > 0.0 && abseta <= 0.9) return 0.9848297238349915;
      else if ( abseta > 0.9 && abseta <= 1.2) return 0.978645384311676;
      else if ( abseta > 1.2 && abseta <= 2.1) return 0.9905462265014648;
      else if ( abseta > 2.1 && abseta <= 2.4) return 0.9802553653717041;
      else return 1.;
    }
    if (this->pt()>30. && this->pt() <= 40.){
      if      ( abseta > 0.0 && abseta <= 0.9) return 0.9861794114112854;
      else if ( abseta > 0.9 && abseta <= 1.2) return 0.9798933267593384;
      else if ( abseta > 1.2 && abseta <= 2.1) return 0.9923668503761292;
      else if ( abseta > 2.1 && abseta <= 2.4) return 0.9785045385360718;
      else return 1.;
    }
    if (this->pt()>40. && this->pt() <= 50.){
      if      ( abseta > 0.0 && abseta <= 0.9) return 0.987443208694458;
      else if ( abseta > 0.9 && abseta <= 1.2) return 0.980233907699585;
      else if ( abseta > 1.2 && abseta <= 2.1) return 0.9927627444267273;
      else if ( abseta > 2.1 && abseta <= 2.4) return 0.9778544902801514;
      else return 1.;
    }
    if (this->pt()>50. && this->pt() <= 60.){
      if      ( abseta > 0.0 && abseta <= 0.9) return 0.9834294319152832;
      else if ( abseta > 0.9 && abseta <= 1.2) return 0.9773300886154175;
      else if ( abseta > 1.2 && abseta <= 2.1) return 0.9886322021484375;
      else if ( abseta > 2.1 && abseta <= 2.4) return 0.9654409885406494;
      else return 1.;
    }
    if (this->pt()>60. && this->pt() <= 120.){
      if      ( abseta > 0.0 && abseta <= 0.9) return 0.9863178730010986;
      else if ( abseta > 0.9 && abseta <= 1.2) return 0.9795225858688354;
      else if ( abseta > 1.2 && abseta <= 2.1) return 0.9950451850891113;
      else if ( abseta > 2.1 && abseta <= 2.4) return 0.9689615368843079;
      else return 1.;
    }
  }

  if (name == "NUM_TightRelIso_DEN_TightID_PAR_pt_spliteta_bin1"){
    float abseta = std::abs(this->eta());
    if (this->pt()>20. && this->pt() <= 25.){
      if      ( abseta > 0.0 && abseta <= 0.9) return 1.0043761730194092;
      else if ( abseta > 0.9 && abseta <= 1.2) return 1.004357933998108;
      else if ( abseta > 1.2 && abseta <= 2.1) return 0.9970762133598328;
      else if ( abseta > 2.1 && abseta <= 2.4) return 0.9957730770111084;
      else return 1.;
    }    
    if (this->pt()>25. && this->pt() <= 30.){
      if      ( abseta > 0.0 && abseta <= 0.9) return 0.9995378255844116;
      else if ( abseta > 0.9 && abseta <= 1.2) return 1.002331256866455;
      else if ( abseta > 1.2 && abseta <= 2.1) return 1.0006532669067383;
      else if ( abseta > 2.1 && abseta <= 2.4) return 0.9939026832580566;
      else return 1.;
    }
    if (this->pt()>30. && this->pt() <= 40.){
      if      ( abseta > 0.0 && abseta <= 0.9) return 1.000901222229004;
      else if ( abseta > 0.9 && abseta <= 1.2) return 1.004658579826355;
      else if ( abseta > 1.2 && abseta <= 2.1) return 1.0023553371429443;
      else if ( abseta > 2.1 && abseta <= 2.4) return 0.997478187084198;
      else return 1.;
    }
    if (this->pt()>40. && this->pt() <= 50.){
      if      ( abseta > 0.0 && abseta <= 0.9) return 0.9986253976821899;
      else if ( abseta > 0.9 && abseta <= 1.2) return 1.0013608932495117;
      else if ( abseta > 1.2 && abseta <= 2.1) return 0.999933660030365;
      else if ( abseta > 2.1 && abseta <= 2.4) return 1.002805233001709;
      else return 1.;
    }
    if (this->pt()>50. && this->pt() <= 60.){
      if      ( abseta > 0.0 && abseta <= 0.9) return 1.0002487897872925;
      else if ( abseta > 0.9 && abseta <= 1.2) return 0.9986217021942139;
      else if ( abseta > 1.2 && abseta <= 2.1) return 1.0002963542938232;
      else if ( abseta > 2.1 && abseta <= 2.4) return 1.0043764114379883;
      else return 1.;
    }
    if (this->pt()>60. && this->pt() <= 120.){
      if      ( abseta > 0.0 && abseta <= 0.9) return 0.9986850023269653;
      else if ( abseta > 0.9 && abseta <= 1.2) return 1.0054655075073242;
      else if ( abseta > 1.2 && abseta <= 2.1) return 1.0004935264587402;
      else if ( abseta > 2.1 && abseta <= 2.4) return 1.0010104179382324;
      else return 1.;
    }
  }
  if (name == "runD_IsoMu20_OR_IsoTkMu20_HLTv4p3_PtEtaBins"){
    float abseta = std::abs(this->eta());
    if (this->pt()>22. && this->pt() <= 25.){
      if      ( abseta > 0.0 && abseta <= 0.9) return 1.0123462677001953;
      else if ( abseta > 0.9 && abseta <= 1.2) return 1.0364969968795776;
      else if ( abseta > 1.2 && abseta <= 2.1) return 1.0033341646194458;
      else if ( abseta > 2.1 && abseta <= 2.4) return 1.0173189640045166;
      else return 1.;
    }    
    if (this->pt()>25. && this->pt() <= 30.){
      if      ( abseta > 0.0 && abseta <= 0.9) return 1.0114353895187378;
      else if ( abseta > 0.9 && abseta <= 1.2) return 1.0190913677215576;
      else if ( abseta > 1.2 && abseta <= 2.1) return 1.0003139972686768;
      else if ( abseta > 2.1 && abseta <= 2.4) return 1.0156738758087158;
      else return 1.;
    }
    if (this->pt()>30. && this->pt() <= 40.){
      if      ( abseta > 0.0 && abseta <= 0.9) return 1.0017163753509521;
      else if ( abseta > 0.9 && abseta <= 1.2) return 1.0062072277069092;
      else if ( abseta > 1.2 && abseta <= 2.1) return 0.9933236837387085;
      else if ( abseta > 2.1 && abseta <= 2.4) return 1.0214488506317139;
      else return 1.;
    }
    if (this->pt()>40. && this->pt() <= 50.){
      if      ( abseta > 0.0 && abseta <= 0.9) return 0.999156653881073;
      else if ( abseta > 0.9 && abseta <= 1.2) return 0.9967776536941528;
      else if ( abseta > 1.2 && abseta <= 2.1) return 0.9920070171356201;
      else if ( abseta > 2.1 && abseta <= 2.4) return 1.0131736993789673;
      else return 1.;
    }
    if (this->pt()>50. && this->pt() <= 60.){
      if      ( abseta > 0.0 && abseta <= 0.9) return 0.9958007335662842;
      else if ( abseta > 0.9 && abseta <= 1.2) return 0.9995957016944885;
      else if ( abseta > 1.2 && abseta <= 2.1) return 0.9910058975219727;
      else if ( abseta > 2.1 && abseta <= 2.4) return 1.0221701860427856;
      else return 1.;
    }
    if (this->pt()>60. && this->pt() <= 120.){
      if      ( abseta > 0.0 && abseta <= 0.9) return 0.9904376864433289;
      else if ( abseta > 0.9 && abseta <= 1.2) return 0.9838433265686035;
      else if ( abseta > 1.2 && abseta <= 2.1) return 0.9919053316116333;
      else if ( abseta > 2.1 && abseta <= 2.4) return 1.0137783288955688;
      else return 1.;
    }
  }
  
   return 1.;
}
