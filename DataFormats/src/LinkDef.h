#ifdef __CINT__
#include "../interface/Run.h"
#include "../interface/Event.h"
#include "../interface/Particle.h"
#include "../interface/MCParticle.h"
#include "../interface/Jet.h"
#include "../interface/GenJet.h"
#include "../interface/PFJet.h"
#include "../interface/Muon.h"
#include "../interface/Electron.h"
#include "../interface/Lepton.h"
#include "../interface/Photon.h"
#include "../interface/MET.h"
#include "../interface/PFMET.h"
#include "../interface/GenEvent.h"
#include "../interface/GenTop.h"
#include "../interface/NPGenEvent.h"
#include "../interface/SpinCorrGen.h"
#include "../interface/Vertex.h"
#include "../interface/HLTInfo.h"
#else
#include "../interface/Run.h"
#include "../interface/Event.h"
#include "../interface/Particle.h"
#include "../interface/MCParticle.h"
#include "../interface/Jet.h"
#include "../interface/GenJet.h"
#include "../interface/PFJet.h"
#include "../interface/Muon.h"
#include "../interface/Electron.h"
#include "../interface/Lepton.h"
#include "../interface/Photon.h"
#include "../interface/MET.h"
#include "../interface/PFMET.h"
#include "../interface/GenEvent.h"
#include "../interface/GenTop.h"
#include "../interface/NPGenEvent.h"
#include "../interface/SpinCorrGen.h"
#include "../interface/Vertex.h"
#include "../interface/HLTInfo.h"
#endif

#ifdef __CINT__
#pragma link off all globals;
#pragma link off all classes;
#pragma link off all functions;

#pragma link C++ class cat::Run;
#pragma link C++ class cat::Event;
#pragma link C++ class cat::Particle;
#pragma link C++ class cat::MCParticle;
#pragma link C++ class cat::Jet;
#pragma link C++ class cat::GenJet;
#pragma link C++ class cat::PFJet;
#pragma link C++ class cat::Muon;
#pragma link C++ class cat::Electron;
#pragma link C++ class cat::Lepton;
#pragma link C++ class cat::Photon;
#pragma link C++ class cat::MET;
#pragma link C++ class cat::PFMET;
#pragma link C++ class cat::GenEvent;
#pragma link C++ class cat::GenTop;
#pragma link C++ class cat::NPGenEvent;
#pragma link C++ class cat::SpinCorrGen;
#pragma link C++ class cat::Vertex;
#pragma link C++ class cat::HLTInfo;

#pragma link C++ struct cat::triggeredObject;

#endif


