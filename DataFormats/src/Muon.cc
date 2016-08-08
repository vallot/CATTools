#include "CATTools/DataFormats/interface/Muon.h"

using namespace cat;

/// default constructor
Muon::Muon() {}

Muon::Muon(const reco::LeafCandidate & aMuon) :
  Lepton( aMuon ),
  isGlobalMuon_(false),
  isSoftMuon_(false),
  normalizedChi2_(-1),
  ipsig_(-1),
  numberOfValidHits_(0),
  numberOfValidMuonHits_(0),
  numberOfMatchedStations_(0),
  numberOfValidPixelHits_(0),
  trackerLayersWithMeasurement_(0)
{}

/// destructor
Muon::~Muon() {
}

float Muon::scaleFactor(const std::string& name, int sign) const {
  if (name == "NUM_TightIDandIPCut_DEN_genTracks_PAR_pt_spliteta_bin1"){
    float abseta = std::abs(this->eta());
    if (this->pt()>20. && this->pt() <= 25.){
      if      ( abseta > 0.0 && abseta <= 0.9) return 0.9752116203308105 + sign*(0.0030660638813280626+0.01);
      else if ( abseta > 0.9 && abseta <= 1.2) return 0.9738101959228516 + sign*(0.004502934246978295+0.01);
      else if ( abseta > 1.2 && abseta <= 2.1) return 0.9983288645744324 + sign*(0.002331323348626783+0.01);
      else if ( abseta > 2.1 && abseta <= 2.4) return 0.9877836108207703 + sign*(0.004915740433340289+0.01);
      else return 1.;
    }    
    if (this->pt()>25. && this->pt() <= 30.){
      if      ( abseta > 0.0 && abseta <= 0.9) return 0.9848297238349915 + sign*(0.0016307213764927449+0.01);
      else if ( abseta > 0.9 && abseta <= 1.2) return 0.978645384311676 + sign*(0.0027064755458685794+0.01);
      else if ( abseta > 1.2 && abseta <= 2.1) return 0.9905462265014648 + sign*(0.001402578599690647+0.01);
      else if ( abseta > 2.1 && abseta <= 2.4) return 0.9802553653717041 + sign*(0.003173276637083633+0.01);
      else return 1.;
    }
    if (this->pt()>30. && this->pt() <= 40.){
      if      ( abseta > 0.0 && abseta <= 0.9) return 0.9861794114112854 + sign*(0.0006187187412138267+0.01);
      else if ( abseta > 0.9 && abseta <= 1.2) return 0.9798933267593384 + sign*(0.001057081371390319+0.01);
      else if ( abseta > 1.2 && abseta <= 2.1) return 0.9923668503761292 + sign*(0.0005653311393042486+0.01);
      else if ( abseta > 2.1 && abseta <= 2.4) return 0.9785045385360718 + sign*(0.0015542030446523895+0.01);
      else return 1.;
    }
    if (this->pt()>40. && this->pt() <= 50.){
      if      ( abseta > 0.0 && abseta <= 0.9) return 0.987443208694458 + sign*(0.000494159746725046+0.01);
      else if ( abseta > 0.9 && abseta <= 1.2) return 0.980233907699585 + sign*(0.000819615406448897+0.01);
      else if ( abseta > 1.2 && abseta <= 2.1) return 0.9927627444267273 + sign*(0.0004155573807947332+0.01);
      else if ( abseta > 2.1 && abseta <= 2.4) return 0.9778544902801514 + sign*(0.001456799997296391+0.01);
      else return 1.;
    }
    if (this->pt()>50. && this->pt() <= 60.){
      if      ( abseta > 0.0 && abseta <= 0.9) return 0.9834294319152832 + sign*(0.0011818999573518245+0.01);
      else if ( abseta > 0.9 && abseta <= 1.2) return 0.9773300886154175 + sign*(0.001955436343316424+0.01);
      else if ( abseta > 1.2 && abseta <= 2.1) return 0.9886322021484375 + sign*(0.0011254961157344963+0.01);
      else if ( abseta > 2.1 && abseta <= 2.4) return 0.9654409885406494 + sign*(0.003709169009223743+0.01);
      else return 1.;
    }
    if (this->pt()>60. && this->pt() <= 120.){
      if      ( abseta > 0.0 && abseta <= 0.9) return 0.9863178730010986 + sign*(0.002073330940717176+0.01);
      else if ( abseta > 0.9 && abseta <= 1.2) return 0.9795225858688354 + sign*(0.0035622593553725837+0.01);
      else if ( abseta > 1.2 && abseta <= 2.1) return 0.9950451850891113 + sign*(0.002673833447209764+0.01);
      else if ( abseta > 2.1 && abseta <= 2.4) return 0.9689615368843079 + sign*(0.011084748199568817+0.01);
      else return 1.;
    }
  }

  if (name == "NUM_TightRelIso_DEN_TightID_PAR_pt_spliteta_bin1"){
    float abseta = std::abs(this->eta());
    if (this->pt()>20. && this->pt() <= 25.){
      if      ( abseta > 0.0 && abseta <= 0.9) return 1.0043761730194092 + sign*(0.003959090391076143+0.01);
      else if ( abseta > 0.9 && abseta <= 1.2) return 1.004357933998108 + sign*(0.006125539530136138+0.01);
      else if ( abseta > 1.2 && abseta <= 2.1) return 0.9970762133598328 + sign*(0.003109125287470401+0.01);
      else if ( abseta > 2.1 && abseta <= 2.4) return 0.9957730770111084 + sign*(0.006137193387970902+0.01);
      else return 1.;
    }    
    if (this->pt()>25. && this->pt() <= 30.){
      if      ( abseta > 0.0 && abseta <= 0.9) return 0.9995378255844116 + sign*(0.0022512071035640673+0.01);
      else if ( abseta > 0.9 && abseta <= 1.2) return 1.002331256866455 + sign*(0.004003683572512011+0.01);
      else if ( abseta > 1.2 && abseta <= 2.1) return 1.0006532669067383 + sign*(0.002067755362435184+0.01);
      else if ( abseta > 2.1 && abseta <= 2.4) return 0.9939026832580566 + sign*(0.004261971076013437+0.01);
      else return 1.;
    }
    if (this->pt()>30. && this->pt() <= 40.){
      if      ( abseta > 0.0 && abseta <= 0.9) return 1.000901222229004 + sign*(0.0007979481788689052+0.01);
      else if ( abseta > 0.9 && abseta <= 1.2) return 1.004658579826355 + sign*(0.0014502638048416372+0.01);
      else if ( abseta > 1.2 && abseta <= 2.1) return 1.0023553371429443 + sign*(0.0008445520691793605+0.01);
      else if ( abseta > 2.1 && abseta <= 2.4) return 0.997478187084198 + sign*(0.001781225374381486+0.01);
      else return 1.;
    }
    if (this->pt()>40. && this->pt() <= 50.){
      if      ( abseta > 0.0 && abseta <= 0.9) return 0.9986253976821899 + sign*(0.0004518361024064332+0.01);
      else if ( abseta > 0.9 && abseta <= 1.2) return 1.0013608932495117 + sign*(0.0004888604573095644+0.01);
      else if ( abseta > 1.2 && abseta <= 2.1) return 0.999933660030365 + sign*(0.0004309914887707696+0.01);
      else if ( abseta > 2.1 && abseta <= 2.4) return 1.002805233001709 + sign*(0.001100242856214239+0.01);
      else return 1.;
    }
    if (this->pt()>50. && this->pt() <= 60.){
      if      ( abseta > 0.0 && abseta <= 0.9) return 1.0002487897872925 + sign*(0.000772847340102783+0.01);
      else if ( abseta > 0.9 && abseta <= 1.2) return 0.9986217021942139 + sign*(0.0012396364566794034+0.01);
      else if ( abseta > 1.2 && abseta <= 2.1) return 1.0002963542938232 + sign*(0.0007614160360063238+0.01);
      else if ( abseta > 2.1 && abseta <= 2.4) return 1.0043764114379883 + sign*(0.001806526581100641+0.01);
      else return 1.;
    }
    if (this->pt()>60. && this->pt() <= 120.){
      if      ( abseta > 0.0 && abseta <= 0.9) return 0.9986850023269653 + sign*(0.0008907575174433545+0.01);
      else if ( abseta > 0.9 && abseta <= 1.2) return 1.0054655075073242 + sign*(0.001589130019220112+0.01);
      else if ( abseta > 1.2 && abseta <= 2.1) return 1.0004935264587402 + sign*(0.0009382223143922724+0.01);
      else if ( abseta > 2.1 && abseta <= 2.4) return 1.0010104179382324 + sign*(0.0022795762936220253+0.01);
      else return 1.;
    }
  }
  
  if (name == "runD_IsoMu20_OR_IsoTkMu20_HLTv4p3_PtEtaBins"){
    float abseta = std::abs(this->eta());
    if (this->pt()>22. && this->pt() <= 25.){
      if      ( abseta > 0.0 && abseta <= 0.9) return 1.0123462677001953 + sign*(0.003915927567865302+0.005);
      else if ( abseta > 0.9 && abseta <= 1.2) return 1.0364969968795776 + sign*(0.007723853105651518+0.005);
      else if ( abseta > 1.2 && abseta <= 2.1) return 1.0033341646194458 + sign*(0.004321950105020405+0.005);
      else if ( abseta > 2.1 && abseta <= 2.4) return 1.0173189640045166 + sign*(0.011553505049544352+0.005);
      else return 1.;
    }    
    if (this->pt()>25. && this->pt() <= 30.){
      if      ( abseta > 0.0 && abseta <= 0.9) return 1.0114353895187378 + sign*(0.002051091269986238+0.005);
      else if ( abseta > 0.9 && abseta <= 1.2) return 1.0190913677215576 + sign*(0.004311463510911091+0.005);
      else if ( abseta > 1.2 && abseta <= 2.1) return 1.0003139972686768 + sign*(0.002556847344445013+0.005);
      else if ( abseta > 2.1 && abseta <= 2.4) return 1.0156738758087158 + sign*(0.006752866642118899+0.005);
      else return 1.;
    }
    if (this->pt()>30. && this->pt() <= 40.){
      if      ( abseta > 0.0 && abseta <= 0.9) return 1.0017163753509521 + sign*(0.0008555351596684594+0.005);
      else if ( abseta > 0.9 && abseta <= 1.2) return 1.0062072277069092 + sign*(0.0018825015620286413+0.005);
      else if ( abseta > 1.2 && abseta <= 2.1) return 0.9933236837387085 + sign*(0.0011911812478541811+0.005);
      else if ( abseta > 2.1 && abseta <= 2.4) return 1.0214488506317139 + sign*(0.00345789723771596+0.005);
      else return 1.;
    }
    if (this->pt()>40. && this->pt() <= 50.){
      if      ( abseta > 0.0 && abseta <= 0.9) return 0.999156653881073 + sign*(0.0003927590496207436+0.005);
      else if ( abseta > 0.9 && abseta <= 1.2) return 0.9967776536941528 + sign*(0.0013975747759820665+0.005);
      else if ( abseta > 1.2 && abseta <= 2.1) return 0.9920070171356201 + sign*(0.0008800556423714225+0.005);
      else if ( abseta > 2.1 && abseta <= 2.4) return 1.0131736993789673 + sign*(0.002758106472665241+0.005);
      else return 1.;
    }
    if (this->pt()>50. && this->pt() <= 60.){
      if      ( abseta > 0.0 && abseta <= 0.9) return 0.9958007335662842 + sign*(0.0013143767249472538+0.005);
      else if ( abseta > 0.9 && abseta <= 1.2) return 0.9995957016944885 + sign*(0.0028673199626467415+0.005);
      else if ( abseta > 1.2 && abseta <= 2.1) return 0.9910058975219727 + sign*(0.0017550948031006969+0.005);
      else if ( abseta > 2.1 && abseta <= 2.4) return 1.0221701860427856 + sign*(0.005536197016181099+0.005);
      else return 1.;
    }
    if (this->pt()>60. && this->pt() <= 120.){
      if      ( abseta > 0.0 && abseta <= 0.9) return 0.9904376864433289 + sign*(0.0016431593847415083+0.005);
      else if ( abseta > 0.9 && abseta <= 1.2) return 0.9838433265686035 + sign*(0.0035930482614946967+0.005);
      else if ( abseta > 1.2 && abseta <= 2.1) return 0.9919053316116333 + sign*(0.002361883767990562+0.005);
      else if ( abseta > 2.1 && abseta <= 2.4) return 1.0137783288955688 + sign*(0.007529513743816103+0.005);
      else return 1.;
    }
  }
  
   return 1.;
}
