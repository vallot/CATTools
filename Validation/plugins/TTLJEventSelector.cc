#include "FWCore/Framework/interface/one/EDFilter.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "CATTools/DataFormats/interface/Muon.h"
#include "CATTools/DataFormats/interface/Electron.h"
#include "CATTools/DataFormats/interface/Jet.h"
#include "CATTools/DataFormats/interface/MET.h"
#include "CATTools/DataFormats/interface/SecVertex.h"

#include "CATTools/CommonTools/interface/TTbarModeDefs.h"
#include "CATTools/CommonTools/interface/ScaleFactorEvaluator.h"
#include "CATTools/CatAnalyzer/interface/TopTriggerSF.h"

#include "DataFormats/Candidate/interface/LeafCandidate.h"
//#include "DataFormats/Candidate/interface/CompositeCandidate.h"
//#include "DataFormats/Candidate/interface/CompositeRefCandidate.h"
#include "DataFormats/Candidate/interface/CompositePtrCandidate.h"

#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "CommonTools/Utils/interface/TFileDirectory.h"
#include "TH1D.h"
#include "TH2F.h"

using namespace std;

namespace cat {

struct ControlPlotsTTLJ
{
  const static int nMaxCutstep;
  typedef TH1D* H1;
  typedef TH2D* H2;

  H1 hCutstep, hCutstepNoweight;
  H2 h2Cutstep, h2CutstepNoweight;

  H1 h0a_vertex_n;

  H1 h0b_vertex_n;
  H1 h0b_met_pt, h0b_met_phi;
  H1 h0b_leptons_n, h0b_leptons_pt, h0b_leptons_eta;
  H1 h0b_jets_n, h0b_jets_pt, h0b_jets_eta, h0b_jets_ht;
  H1 h0b_bjets_n;

  H1 h0c_vertex_n;
  H1 h0c_met_pt, h0c_met_phi;
  H1 h0c_leptons_n;
  H1 h0c_jets_n, h0c_jets_pt, h0c_jets_eta, h0c_jets_ht;
  H1 h0c_bjets_n;

  H1 h1_vertex_n;
  H1 h1_met_pt, h1_met_phi;
  H1 h1_leptons_n;
  H1 h1_lepton1_pt, h1_lepton1_eta, h1_lepton1_phi, h1_lepton1_q;
  H1 h1_jets_n, h1_jets_pt, h1_jets_eta, h1_jets_ht;
  H1 h1_jet1_m, h1_jet1_pt, h1_jet1_eta, h1_jet1_phi, h1_jet1_btag;
  H1 h1_jet2_m, h1_jet2_pt, h1_jet2_eta, h1_jet2_phi, h1_jet2_btag;
  H1 h1_jet3_m, h1_jet3_pt, h1_jet3_eta, h1_jet3_phi, h1_jet3_btag;
  H1 h1_jet4_m, h1_jet4_pt, h1_jet4_eta, h1_jet4_phi, h1_jet4_btag;
  H1 h1_jet5_m, h1_jet5_pt, h1_jet5_eta, h1_jet5_phi, h1_jet5_btag;
  H1 h1_jet6_m, h1_jet6_pt, h1_jet6_eta, h1_jet6_phi, h1_jet6_btag;
  H1 h1_bjets_n;
  H1 h1_event_st;

  H1 h2a_vertex_n;
  H1 h2a_met_pt, h2a_met_phi;
  H1 h2a_leptons_n;
  H1 h2a_lepton1_pt, h2a_lepton1_eta, h2a_lepton1_phi, h2a_lepton1_q;
  H1 h2a_jets_n, h2a_jets_pt, h2a_jets_eta, h2a_jets_ht;
  H1 h2a_jet1_m, h2a_jet1_pt, h2a_jet1_eta, h2a_jet1_phi, h2a_jet1_btag;
  H1 h2a_jet2_m, h2a_jet2_pt, h2a_jet2_eta, h2a_jet2_phi, h2a_jet2_btag;
  H1 h2a_jet3_m, h2a_jet3_pt, h2a_jet3_eta, h2a_jet3_phi, h2a_jet3_btag;
  H1 h2a_jet4_m, h2a_jet4_pt, h2a_jet4_eta, h2a_jet4_phi, h2a_jet4_btag;
  H1 h2a_jet5_m, h2a_jet5_pt, h2a_jet5_eta, h2a_jet5_phi, h2a_jet5_btag;
  H1 h2a_jet6_m, h2a_jet6_pt, h2a_jet6_eta, h2a_jet6_phi, h2a_jet6_btag;
  H1 h2a_bjets_n;
  H1 h2a_event_st;

  H1 h2b_vertex_n;
  H1 h2b_met_pt, h2b_met_phi;
  H1 h2b_leptons_n;
  H1 h2b_lepton1_pt, h2b_lepton1_eta, h2b_lepton1_phi, h2b_lepton1_q;
  H1 h2b_jets_n, h2b_jets_pt, h2b_jets_eta, h2b_jets_ht;
  H1 h2b_jet1_m, h2b_jet1_pt, h2b_jet1_eta, h2b_jet1_phi, h2b_jet1_btag;
  H1 h2b_jet2_m, h2b_jet2_pt, h2b_jet2_eta, h2b_jet2_phi, h2b_jet2_btag;
  H1 h2b_jet3_m, h2b_jet3_pt, h2b_jet3_eta, h2b_jet3_phi, h2b_jet3_btag;
  H1 h2b_jet4_m, h2b_jet4_pt, h2b_jet4_eta, h2b_jet4_phi, h2b_jet4_btag;
  H1 h2b_jet5_m, h2b_jet5_pt, h2b_jet5_eta, h2b_jet5_phi, h2b_jet5_btag;
  H1 h2b_jet6_m, h2b_jet6_pt, h2b_jet6_eta, h2b_jet6_phi, h2b_jet6_btag;
  H1 h2b_bjets_n;
  H1 h2b_event_st;

  H1 h3_vertex_n;
  H1 h3_met_pt, h3_met_phi;
  H1 h3_jets_n, h3_jets_pt, h3_jets_eta, h3_jets_ht;
  H1 h3_jet1_m, h3_jet1_pt, h3_jet1_eta, h3_jet1_phi, h3_jet1_btag;
  H1 h3_jet2_m, h3_jet2_pt, h3_jet2_eta, h3_jet2_phi, h3_jet2_btag;
  H1 h3_jet3_m, h3_jet3_pt, h3_jet3_eta, h3_jet3_phi, h3_jet3_btag;
  H1 h3_jet4_m, h3_jet4_pt, h3_jet4_eta, h3_jet4_phi, h3_jet4_btag;
  H1 h3_jet5_m, h3_jet5_pt, h3_jet5_eta, h3_jet5_phi, h3_jet5_btag;
  H1 h3_jet6_m, h3_jet6_pt, h3_jet6_eta, h3_jet6_phi, h3_jet6_btag;
  H1 h3_bjets_n;
  H1 h3_event_st;

  H1 h4_vertex_n;
  H1 h4_met_pt, h4_met_phi;
  H1 h4_jets_n, h4_jets_pt, h4_jets_eta, h4_jets_ht;
  H1 h4_jet1_m, h4_jet1_pt, h4_jet1_eta, h4_jet1_phi, h4_jet1_btag;
  H1 h4_jet2_m, h4_jet2_pt, h4_jet2_eta, h4_jet2_phi, h4_jet2_btag;
  H1 h4_jet3_m, h4_jet3_pt, h4_jet3_eta, h4_jet3_phi, h4_jet3_btag;
  H1 h4_jet4_m, h4_jet4_pt, h4_jet4_eta, h4_jet4_phi, h4_jet4_btag;
  H1 h4_jet5_m, h4_jet5_pt, h4_jet5_eta, h4_jet5_phi, h4_jet5_btag;
  H1 h4_jet6_m, h4_jet6_pt, h4_jet6_eta, h4_jet6_phi, h4_jet6_btag;
  H1 h4_bjets_n;
  H1 h4_event_st;

  void book(TFileDirectory&& dir)
  {
    const double maxeta = 3;
    const double pi = 3.141592;

    // There are step0a, step0b and step0c cut steps
    // By putting step0a to underflow bin and step0b to -1, step0c to 0,
    // We can start cut steps from 1.
    hCutstep = dir.make<TH1D>("cutstep", "cutstep", nMaxCutstep, -2, nMaxCutstep-2);
    hCutstepNoweight = dir.make<TH1D>("cutstepNoweight", "cutstepNoweight", nMaxCutstep, -2, nMaxCutstep-2);
    h2Cutstep = dir.make<TH2D>("cutcorrelation", "cutcorrelation", nMaxCutstep, -2, nMaxCutstep-2, nMaxCutstep, -2, nMaxCutstep-2);
    h2CutstepNoweight = dir.make<TH2D>("cutcorrelationNoweight", "cutcorrelationNoweight", nMaxCutstep, -2, nMaxCutstep-2, nMaxCutstep, -2, nMaxCutstep-2);

    hCutstep->GetXaxis()->SetBinLabel(1, "S0a all event");
    hCutstep->GetXaxis()->SetBinLabel(2, "S0b Trigger");
    hCutstep->GetXaxis()->SetBinLabel(3, "S0c Event filter");
    hCutstep->GetXaxis()->SetBinLabel(4, "S1  Signal lepton");
    hCutstep->GetXaxis()->SetBinLabel(5, "S2a Veto muon");
    hCutstep->GetXaxis()->SetBinLabel(6, "S2b Veto electron");
    hCutstep->GetXaxis()->SetBinLabel(7, "S3a nJet1");
    hCutstep->GetXaxis()->SetBinLabel(8, "S4  nBJet1");

    hCutstepNoweight->GetXaxis()->SetBinLabel(1, "S0a all event");
    hCutstepNoweight->GetXaxis()->SetBinLabel(2, "S0b Trigger");
    hCutstepNoweight->GetXaxis()->SetBinLabel(3, "S0c Event filter");
    hCutstepNoweight->GetXaxis()->SetBinLabel(4, "S1  Signal lepton");
    hCutstepNoweight->GetXaxis()->SetBinLabel(5, "S2a Veto muon");
    hCutstepNoweight->GetXaxis()->SetBinLabel(6, "S2b Veto electron");
    hCutstepNoweight->GetXaxis()->SetBinLabel(7, "S3a nJet1");
    hCutstepNoweight->GetXaxis()->SetBinLabel(8, "S4  nBJet1");

    h2Cutstep->GetXaxis()->SetBinLabel(1, "S0a all event");
    h2Cutstep->GetXaxis()->SetBinLabel(2, "S0b Trigger");
    h2Cutstep->GetXaxis()->SetBinLabel(3, "S0c Event filter");
    h2Cutstep->GetXaxis()->SetBinLabel(4, "S1  Signal lepton");
    h2Cutstep->GetXaxis()->SetBinLabel(5, "S2a Veto muon");
    h2Cutstep->GetXaxis()->SetBinLabel(6, "S2b Veto electron");
    h2Cutstep->GetXaxis()->SetBinLabel(7, "S3a nJet1");
    h2Cutstep->GetXaxis()->SetBinLabel(8, "S4  nBJet1");

    h2Cutstep->GetYaxis()->SetBinLabel(1, "S0a all event");
    h2Cutstep->GetYaxis()->SetBinLabel(2, "S0b Trigger");
    h2Cutstep->GetYaxis()->SetBinLabel(3, "S0c Event filter");
    h2Cutstep->GetYaxis()->SetBinLabel(4, "S1  Signal lepton");
    h2Cutstep->GetYaxis()->SetBinLabel(5, "S2a Veto muon");
    h2Cutstep->GetYaxis()->SetBinLabel(6, "S2b Veto electron");
    h2Cutstep->GetYaxis()->SetBinLabel(7, "S3a nJet1");
    h2Cutstep->GetYaxis()->SetBinLabel(8, "S4  nBJet1");

    h2CutstepNoweight->GetXaxis()->SetBinLabel(1, "S0a all event");
    h2CutstepNoweight->GetXaxis()->SetBinLabel(2, "S0b Trigger");
    h2CutstepNoweight->GetXaxis()->SetBinLabel(3, "S0c Event filter");
    h2CutstepNoweight->GetXaxis()->SetBinLabel(4, "S1  Signal lepton");
    h2CutstepNoweight->GetXaxis()->SetBinLabel(5, "S2a Veto muon");
    h2CutstepNoweight->GetXaxis()->SetBinLabel(6, "S2b Veto electron");
    h2CutstepNoweight->GetXaxis()->SetBinLabel(7, "S3a nJet1");
    h2CutstepNoweight->GetXaxis()->SetBinLabel(8, "S4  nBJet1");

    h2CutstepNoweight->GetYaxis()->SetBinLabel(1, "S0a all event");
    h2CutstepNoweight->GetYaxis()->SetBinLabel(2, "S0b Trigger");
    h2CutstepNoweight->GetYaxis()->SetBinLabel(3, "S0c Event filter");
    h2CutstepNoweight->GetYaxis()->SetBinLabel(4, "S1  Signal lepton");
    h2CutstepNoweight->GetYaxis()->SetBinLabel(5, "S2a Veto muon");
    h2CutstepNoweight->GetYaxis()->SetBinLabel(6, "S2b Veto electron");
    h2CutstepNoweight->GetYaxis()->SetBinLabel(7, "S3a nJet1");
    h2CutstepNoweight->GetYaxis()->SetBinLabel(8, "S4  nBJet1");

    auto subdir = dir.mkdir("step0a");
    h0a_vertex_n = subdir.make<TH1D>("vertex_n", "vertex_n", 100, 0, 100);

    subdir = dir.mkdir("step0b");
    h0b_vertex_n = subdir.make<TH1D>("vertex_n", "vertex_n", 100, 0, 100);
    h0b_met_pt = subdir.make<TH1D>("met_pt", "met_pt", 1000, 0, 1000);
    h0b_met_phi = subdir.make<TH1D>("met_phi", "met_phi", 100, -pi, pi);
    h0b_leptons_n = subdir.make<TH1D>("leptons_n", "leptons_n", 10, 0, 10);
    h0b_jets_n = subdir.make<TH1D>("jets_n", "jets_n", 10, 0, 10);
    h0b_jets_pt  = subdir.make<TH1D>("jets_pt", "jets_pt", 1000, 0, 1000);
    h0b_jets_eta = subdir.make<TH1D>("jets_eta", "jets_eta", 100, -maxeta, maxeta);
    h0b_jets_ht = subdir.make<TH1D>("jets_ht", "jets_ht", 1000, 0, 1000);
    h0b_bjets_n = subdir.make<TH1D>("bjets_n", "bjets_n", 10, 0, 10);

    subdir = dir.mkdir("step0c");
    h0c_vertex_n = subdir.make<TH1D>("vertex_n", "vertex_n", 100, 0, 100);
    h0c_met_pt = subdir.make<TH1D>("met_pt", "met_pt", 1000, 0, 1000);
    h0c_met_phi = subdir.make<TH1D>("met_phi", "met_phi", 100, -pi, pi);
    h0c_leptons_n = subdir.make<TH1D>("leptons_n", "leptons_n", 10, 0, 10);
    h0c_jets_n = subdir.make<TH1D>("jets_n", "jets_n", 10, 0, 10);
    h0c_jets_pt  = subdir.make<TH1D>("jets_pt", "jets_pt", 1000, 0, 1000);
    h0c_jets_eta = subdir.make<TH1D>("jets_eta", "jets_eta", 100, -maxeta, maxeta);
    h0c_jets_ht = subdir.make<TH1D>("jets_ht", "jets_ht", 1000, 0, 1000);
    h0c_bjets_n = subdir.make<TH1D>("bjets_n", "bjets_n", 10, 0, 10);

    subdir = dir.mkdir("step1");
    h1_vertex_n = subdir.make<TH1D>("vertex_n", "vertex_n", 100, 0, 100);
    h1_met_pt = subdir.make<TH1D>("met_pt", "met_pt", 1000, 0, 1000);
    h1_met_phi = subdir.make<TH1D>("met_phi", "met_phi", 100, -pi, pi);
    h1_leptons_n = subdir.make<TH1D>("leptons_n", "leptons_n", 10, 0, 10);

    h1_lepton1_pt  = subdir.make<TH1D>("lepton1_pt", "lepton1_pt", 1000, 0, 1000);
    h1_lepton1_eta = subdir.make<TH1D>("lepton1_eta", "lepton1_eta", 100, -maxeta, maxeta);
    h1_lepton1_phi = subdir.make<TH1D>("lepton1_phi", "lepton1_phi", 100, -pi, pi);
    h1_lepton1_q   = subdir.make<TH1D>("lepton1_q", "lepton1_q", 3, -1.5, 1.5);

    h1_jets_n = subdir.make<TH1D>("jets_n", "jets_n", 10, 0, 10);
    h1_jets_pt  = subdir.make<TH1D>("jets_pt", "jets_pt", 1000, 0, 1000);
    h1_jets_eta = subdir.make<TH1D>("jets_eta", "jets_eta", 100, -maxeta, maxeta);
    h1_jets_ht = subdir.make<TH1D>("jets_ht", "jets_ht", 1000, 0, 1000);

    h1_jet1_m   = subdir.make<TH1D>("jet1_m", "jet1_m", 500, 0, 500);
    h1_jet1_pt  = subdir.make<TH1D>("jet1_pt", "jet1_pt", 1000, 0, 1000);
    h1_jet1_eta = subdir.make<TH1D>("jet1_eta", "jet1_eta", 100, -maxeta, maxeta);
    h1_jet1_phi = subdir.make<TH1D>("jet1_phi", "jet1_phi", 100, -pi, pi);
    h1_jet1_btag = subdir.make<TH1D>("jet1_btag", "jet1_btag", 100, 0, 1);

    h1_jet2_m   = subdir.make<TH1D>("jet2_m", "jet2_m", 500, 0, 500);
    h1_jet2_pt  = subdir.make<TH1D>("jet2_pt", "jet2_pt", 1000, 0, 1000);
    h1_jet2_eta = subdir.make<TH1D>("jet2_eta", "jet2_eta", 100, -maxeta, maxeta);
    h1_jet2_phi = subdir.make<TH1D>("jet2_phi", "jet2_phi", 100, -pi, pi);
    h1_jet2_btag = subdir.make<TH1D>("jet2_btag", "jet2_btag", 100, 0, 1);

    h1_jet3_m   = subdir.make<TH1D>("jet3_m", "jet3_m", 500, 0, 500);
    h1_jet3_pt  = subdir.make<TH1D>("jet3_pt", "jet3_pt", 1000, 0, 1000);
    h1_jet3_eta = subdir.make<TH1D>("jet3_eta", "jet3_eta", 100, -maxeta, maxeta);
    h1_jet3_phi = subdir.make<TH1D>("jet3_phi", "jet3_phi", 100, -pi, pi);
    h1_jet3_btag = subdir.make<TH1D>("jet3_btag", "jet3_btag", 100, 0, 1);

    h1_jet4_m   = subdir.make<TH1D>("jet4_m", "jet4_m", 500, 0, 500);
    h1_jet4_pt  = subdir.make<TH1D>("jet4_pt", "jet4_pt", 1000, 0, 1000);
    h1_jet4_eta = subdir.make<TH1D>("jet4_eta", "jet4_eta", 100, -maxeta, maxeta);
    h1_jet4_phi = subdir.make<TH1D>("jet4_phi", "jet4_phi", 100, -pi, pi);
    h1_jet4_btag = subdir.make<TH1D>("jet4_btag", "jet4_btag", 100, 0, 1);

    h1_jet5_m   = subdir.make<TH1D>("jet5_m", "jet5_m", 500, 0, 500);
    h1_jet5_pt  = subdir.make<TH1D>("jet5_pt", "jet5_pt", 1000, 0, 1000);
    h1_jet5_eta = subdir.make<TH1D>("jet5_eta", "jet5_eta", 100, -maxeta, maxeta);
    h1_jet5_phi = subdir.make<TH1D>("jet5_phi", "jet5_phi", 100, -pi, pi);
    h1_jet5_btag = subdir.make<TH1D>("jet5_btag", "jet5_btag", 100, 0, 1);

    h1_jet6_m   = subdir.make<TH1D>("jet6_m", "jet6_m", 500, 0, 500);
    h1_jet6_pt  = subdir.make<TH1D>("jet6_pt", "jet6_pt", 1000, 0, 1000);
    h1_jet6_eta = subdir.make<TH1D>("jet6_eta", "jet6_eta", 100, -maxeta, maxeta);
    h1_jet6_phi = subdir.make<TH1D>("jet6_phi", "jet6_phi", 100, -pi, pi);
    h1_jet6_btag = subdir.make<TH1D>("jet6_btag", "jet6_btag", 100, 0, 1);

    h1_bjets_n = subdir.make<TH1D>("bjets_n", "bjets_n", 10, 0, 10);

    h1_event_st = subdir.make<TH1D>("event_st", "event_st", 1000, 0, 1000);

    subdir = dir.mkdir("step2a");
    h2a_vertex_n = subdir.make<TH1D>("vertex_n", "vertex_n", 100, 0, 100);
    h2a_met_pt = subdir.make<TH1D>("met_pt", "met_pt", 1000, 0, 1000);
    h2a_met_phi = subdir.make<TH1D>("met_phi", "met_phi", 100, -pi, pi);
    h2a_leptons_n = subdir.make<TH1D>("leptons_n", "leptons_n", 10, 0, 10);

    h2a_lepton1_pt  = subdir.make<TH1D>("lepton1_pt", "lepton1_pt", 1000, 0, 1000);
    h2a_lepton1_eta = subdir.make<TH1D>("lepton1_eta", "lepton1_eta", 100, -maxeta, maxeta);
    h2a_lepton1_phi = subdir.make<TH1D>("lepton1_phi", "lepton1_phi", 100, -pi, pi);
    h2a_lepton1_q   = subdir.make<TH1D>("lepton1_q", "lepton1_q", 3, -1.5, 1.5);

    h2a_jets_n = subdir.make<TH1D>("jets_n", "jets_n", 10, 0, 10);
    h2a_jets_pt  = subdir.make<TH1D>("jets_pt", "jets_pt", 1000, 0, 1000);
    h2a_jets_eta = subdir.make<TH1D>("jets_eta", "jets_eta", 100, -maxeta, maxeta);
    h2a_jets_ht = subdir.make<TH1D>("jets_ht", "jets_ht", 1000, 0, 1000);

    h2a_jet1_m   = subdir.make<TH1D>("jet1_m", "jet1_m", 500, 0, 500);
    h2a_jet1_pt  = subdir.make<TH1D>("jet1_pt", "jet1_pt", 1000, 0, 1000);
    h2a_jet1_eta = subdir.make<TH1D>("jet1_eta", "jet1_eta", 100, -maxeta, maxeta);
    h2a_jet1_phi = subdir.make<TH1D>("jet1_phi", "jet1_phi", 100, -pi, pi);
    h2a_jet1_btag = subdir.make<TH1D>("jet1_btag", "jet1_btag", 100, 0, 1);

    h2a_jet2_m   = subdir.make<TH1D>("jet2_m", "jet2_m", 500, 0, 500);
    h2a_jet2_pt  = subdir.make<TH1D>("jet2_pt", "jet2_pt", 1000, 0, 1000);
    h2a_jet2_eta = subdir.make<TH1D>("jet2_eta", "jet2_eta", 100, -maxeta, maxeta);
    h2a_jet2_phi = subdir.make<TH1D>("jet2_phi", "jet2_phi", 100, -pi, pi);
    h2a_jet2_btag = subdir.make<TH1D>("jet2_btag", "jet2_btag", 100, 0, 1);

    h2a_jet3_m   = subdir.make<TH1D>("jet3_m", "jet3_m", 500, 0, 500);
    h2a_jet3_pt  = subdir.make<TH1D>("jet3_pt", "jet3_pt", 1000, 0, 1000);
    h2a_jet3_eta = subdir.make<TH1D>("jet3_eta", "jet3_eta", 100, -maxeta, maxeta);
    h2a_jet3_phi = subdir.make<TH1D>("jet3_phi", "jet3_phi", 100, -pi, pi);
    h2a_jet3_btag = subdir.make<TH1D>("jet3_btag", "jet3_btag", 100, 0, 1);

    h2a_jet4_m   = subdir.make<TH1D>("jet4_m", "jet4_m", 500, 0, 500);
    h2a_jet4_pt  = subdir.make<TH1D>("jet4_pt", "jet4_pt", 1000, 0, 1000);
    h2a_jet4_eta = subdir.make<TH1D>("jet4_eta", "jet4_eta", 100, -maxeta, maxeta);
    h2a_jet4_phi = subdir.make<TH1D>("jet4_phi", "jet4_phi", 100, -pi, pi);
    h2a_jet4_btag = subdir.make<TH1D>("jet4_btag", "jet4_btag", 100, 0, 1);

    h2a_jet5_m   = subdir.make<TH1D>("jet5_m", "jet5_m", 500, 0, 500);
    h2a_jet5_pt  = subdir.make<TH1D>("jet5_pt", "jet5_pt", 1000, 0, 1000);
    h2a_jet5_eta = subdir.make<TH1D>("jet5_eta", "jet5_eta", 100, -maxeta, maxeta);
    h2a_jet5_phi = subdir.make<TH1D>("jet5_phi", "jet5_phi", 100, -pi, pi);
    h2a_jet5_btag = subdir.make<TH1D>("jet5_btag", "jet5_btag", 100, 0, 1);

    h2a_jet6_m   = subdir.make<TH1D>("jet6_m", "jet6_m", 500, 0, 500);
    h2a_jet6_pt  = subdir.make<TH1D>("jet6_pt", "jet6_pt", 1000, 0, 1000);
    h2a_jet6_eta = subdir.make<TH1D>("jet6_eta", "jet6_eta", 100, -maxeta, maxeta);
    h2a_jet6_phi = subdir.make<TH1D>("jet6_phi", "jet6_phi", 100, -pi, pi);
    h2a_jet6_btag = subdir.make<TH1D>("jet6_btag", "jet6_btag", 100, 0, 1);

    h2a_bjets_n = subdir.make<TH1D>("bjets_n", "bjets_n", 10, 0, 10);

    h2a_event_st = subdir.make<TH1D>("event_st", "event_st", 1000, 0, 1000);

    subdir = dir.mkdir("step2b");
    h2b_vertex_n = subdir.make<TH1D>("vertex_n", "vertex_n", 100, 0, 100);
    h2b_met_pt = subdir.make<TH1D>("met_pt", "met_pt", 1000, 0, 1000);
    h2b_met_phi = subdir.make<TH1D>("met_phi", "met_phi", 100, -pi, pi);
    h2b_leptons_n = subdir.make<TH1D>("leptons_n", "leptons_n", 10, 0, 10);

    h2b_lepton1_pt  = subdir.make<TH1D>("lepton1_pt", "lepton1_pt", 1000, 0, 1000);
    h2b_lepton1_eta = subdir.make<TH1D>("lepton1_eta", "lepton1_eta", 100, -maxeta, maxeta);
    h2b_lepton1_phi = subdir.make<TH1D>("lepton1_phi", "lepton1_phi", 100, -pi, pi);
    h2b_lepton1_q   = subdir.make<TH1D>("lepton1_q", "lepton1_q", 3, -1.5, 1.5);

    h2b_jets_n = subdir.make<TH1D>("jets_n", "jets_n", 10, 0, 10);
    h2b_jets_pt  = subdir.make<TH1D>("jets_pt", "jets_pt", 1000, 0, 1000);
    h2b_jets_eta = subdir.make<TH1D>("jets_eta", "jets_eta", 100, -maxeta, maxeta);
    h2b_jets_ht = subdir.make<TH1D>("jets_ht", "jets_ht", 1000, 0, 1000);

    h2b_jet1_m   = subdir.make<TH1D>("jet1_m", "jet1_m", 500, 0, 500);
    h2b_jet1_pt  = subdir.make<TH1D>("jet1_pt", "jet1_pt", 1000, 0, 1000);
    h2b_jet1_eta = subdir.make<TH1D>("jet1_eta", "jet1_eta", 100, -maxeta, maxeta);
    h2b_jet1_phi = subdir.make<TH1D>("jet1_phi", "jet1_phi", 100, -pi, pi);
    h2b_jet1_btag = subdir.make<TH1D>("jet1_btag", "jet1_btag", 100, 0, 1);

    h2b_jet2_m   = subdir.make<TH1D>("jet2_m", "jet2_m", 500, 0, 500);
    h2b_jet2_pt  = subdir.make<TH1D>("jet2_pt", "jet2_pt", 1000, 0, 1000);
    h2b_jet2_eta = subdir.make<TH1D>("jet2_eta", "jet2_eta", 100, -maxeta, maxeta);
    h2b_jet2_phi = subdir.make<TH1D>("jet2_phi", "jet2_phi", 100, -pi, pi);
    h2b_jet2_btag = subdir.make<TH1D>("jet2_btag", "jet2_btag", 100, 0, 1);

    h2b_jet3_m   = subdir.make<TH1D>("jet3_m", "jet3_m", 500, 0, 500);
    h2b_jet3_pt  = subdir.make<TH1D>("jet3_pt", "jet3_pt", 1000, 0, 1000);
    h2b_jet3_eta = subdir.make<TH1D>("jet3_eta", "jet3_eta", 100, -maxeta, maxeta);
    h2b_jet3_phi = subdir.make<TH1D>("jet3_phi", "jet3_phi", 100, -pi, pi);
    h2b_jet3_btag = subdir.make<TH1D>("jet3_btag", "jet3_btag", 100, 0, 1);

    h2b_jet4_m   = subdir.make<TH1D>("jet4_m", "jet4_m", 500, 0, 500);
    h2b_jet4_pt  = subdir.make<TH1D>("jet4_pt", "jet4_pt", 1000, 0, 1000);
    h2b_jet4_eta = subdir.make<TH1D>("jet4_eta", "jet4_eta", 100, -maxeta, maxeta);
    h2b_jet4_phi = subdir.make<TH1D>("jet4_phi", "jet4_phi", 100, -pi, pi);
    h2b_jet4_btag = subdir.make<TH1D>("jet4_btag", "jet4_btag", 100, 0, 1);

    h2b_jet5_m   = subdir.make<TH1D>("jet5_m", "jet5_m", 500, 0, 500);
    h2b_jet5_pt  = subdir.make<TH1D>("jet5_pt", "jet5_pt", 1000, 0, 1000);
    h2b_jet5_eta = subdir.make<TH1D>("jet5_eta", "jet5_eta", 100, -maxeta, maxeta);
    h2b_jet5_phi = subdir.make<TH1D>("jet5_phi", "jet5_phi", 100, -pi, pi);
    h2b_jet5_btag = subdir.make<TH1D>("jet5_btag", "jet5_btag", 100, 0, 1);

    h2b_jet6_m   = subdir.make<TH1D>("jet6_m", "jet6_m", 500, 0, 500);
    h2b_jet6_pt  = subdir.make<TH1D>("jet6_pt", "jet6_pt", 1000, 0, 1000);
    h2b_jet6_eta = subdir.make<TH1D>("jet6_eta", "jet6_eta", 100, -maxeta, maxeta);
    h2b_jet6_phi = subdir.make<TH1D>("jet6_phi", "jet6_phi", 100, -pi, pi);
    h2b_jet6_btag = subdir.make<TH1D>("jet6_btag", "jet6_btag", 100, 0, 1);

    h2b_bjets_n = subdir.make<TH1D>("bjets_n", "bjets_n", 10, 0, 10);

    h2b_event_st = subdir.make<TH1D>("event_st", "event_st", 1000, 0, 1000);

    subdir = dir.mkdir("step3");
    h3_vertex_n = subdir.make<TH1D>("vertex_n", "vertex_n", 100, 0, 100);
    h3_met_pt = subdir.make<TH1D>("met_pt", "met_pt", 1000, 0, 1000);
    h3_met_phi = subdir.make<TH1D>("met_phi", "met_phi", 100, -pi, pi);

    h3_jets_n = subdir.make<TH1D>("jets_n", "jets_n", 10, 0, 10);
    h3_jets_pt  = subdir.make<TH1D>("jets_pt", "jets_pt", 1000, 0, 1000);
    h3_jets_eta = subdir.make<TH1D>("jets_eta", "jets_eta", 100, -maxeta, maxeta);
    h3_jets_ht = subdir.make<TH1D>("jets_ht", "jets_ht", 1000, 0, 1000);

    h3_jet1_m   = subdir.make<TH1D>("jet1_m", "jet1_m", 500, 0, 500);
    h3_jet1_pt  = subdir.make<TH1D>("jet1_pt", "jet1_pt", 1000, 0, 1000);
    h3_jet1_eta = subdir.make<TH1D>("jet1_eta", "jet1_eta", 100, -maxeta, maxeta);
    h3_jet1_phi = subdir.make<TH1D>("jet1_phi", "jet1_phi", 100, -pi, pi);
    h3_jet1_btag = subdir.make<TH1D>("jet1_btag", "jet1_btag", 100, 0, 1);

    h3_jet2_m   = subdir.make<TH1D>("jet2_m", "jet2_m", 500, 0, 500);
    h3_jet2_pt  = subdir.make<TH1D>("jet2_pt", "jet2_pt", 1000, 0, 1000);
    h3_jet2_eta = subdir.make<TH1D>("jet2_eta", "jet2_eta", 100, -maxeta, maxeta);
    h3_jet2_phi = subdir.make<TH1D>("jet2_phi", "jet2_phi", 100, -pi, pi);
    h3_jet2_btag = subdir.make<TH1D>("jet2_btag", "jet2_btag", 100, 0, 1);

    h3_jet3_m   = subdir.make<TH1D>("jet3_m", "jet3_m", 500, 0, 500);
    h3_jet3_pt  = subdir.make<TH1D>("jet3_pt", "jet3_pt", 1000, 0, 1000);
    h3_jet3_eta = subdir.make<TH1D>("jet3_eta", "jet3_eta", 100, -maxeta, maxeta);
    h3_jet3_phi = subdir.make<TH1D>("jet3_phi", "jet3_phi", 100, -pi, pi);
    h3_jet3_btag = subdir.make<TH1D>("jet3_btag", "jet3_btag", 100, 0, 1);

    h3_jet4_m   = subdir.make<TH1D>("jet4_m", "jet4_m", 500, 0, 500);
    h3_jet4_pt  = subdir.make<TH1D>("jet4_pt", "jet4_pt", 1000, 0, 1000);
    h3_jet4_eta = subdir.make<TH1D>("jet4_eta", "jet4_eta", 100, -maxeta, maxeta);
    h3_jet4_phi = subdir.make<TH1D>("jet4_phi", "jet4_phi", 100, -pi, pi);
    h3_jet4_btag = subdir.make<TH1D>("jet4_btag", "jet4_btag", 100, 0, 1);

    h3_jet5_m   = subdir.make<TH1D>("jet5_m", "jet5_m", 500, 0, 500);
    h3_jet5_pt  = subdir.make<TH1D>("jet5_pt", "jet5_pt", 1000, 0, 1000);
    h3_jet5_eta = subdir.make<TH1D>("jet5_eta", "jet5_eta", 100, -maxeta, maxeta);
    h3_jet5_phi = subdir.make<TH1D>("jet5_phi", "jet5_phi", 100, -pi, pi);
    h3_jet5_btag = subdir.make<TH1D>("jet5_btag", "jet5_btag", 100, 0, 1);

    h3_jet6_m   = subdir.make<TH1D>("jet6_m", "jet6_m", 500, 0, 500);
    h3_jet6_pt  = subdir.make<TH1D>("jet6_pt", "jet6_pt", 1000, 0, 1000);
    h3_jet6_eta = subdir.make<TH1D>("jet6_eta", "jet6_eta", 100, -maxeta, maxeta);
    h3_jet6_phi = subdir.make<TH1D>("jet6_phi", "jet6_phi", 100, -pi, pi);
    h3_jet6_btag = subdir.make<TH1D>("jet6_btag", "jet6_btag", 100, 0, 1);

    h3_bjets_n = subdir.make<TH1D>("bjets_n", "bjets_n", 10, 0, 10);

    h3_event_st = subdir.make<TH1D>("event_st", "event_st", 1000, 0, 1000);

    subdir = dir.mkdir("step4");
    h4_vertex_n = subdir.make<TH1D>("vertex_n", "vertex_n", 100, 0, 100);
    h4_met_pt = subdir.make<TH1D>("met_pt", "met_pt", 1000, 0, 1000);
    h4_met_phi = subdir.make<TH1D>("met_phi", "met_phi", 100, -pi, pi);

    h4_jets_n = subdir.make<TH1D>("jets_n", "jets_n", 10, 0, 10);
    h4_jets_pt  = subdir.make<TH1D>("jets_pt", "jets_pt", 1000, 0, 1000);
    h4_jets_eta = subdir.make<TH1D>("jets_eta", "jets_eta", 100, -maxeta, maxeta);
    h4_jets_ht = subdir.make<TH1D>("jets_ht", "jets_ht", 1000, 0, 1000);

    h4_jet1_m   = subdir.make<TH1D>("jet1_m", "jet1_m", 500, 0, 500);
    h4_jet1_pt  = subdir.make<TH1D>("jet1_pt", "jet1_pt", 1000, 0, 1000);
    h4_jet1_eta = subdir.make<TH1D>("jet1_eta", "jet1_eta", 100, -maxeta, maxeta);
    h4_jet1_phi = subdir.make<TH1D>("jet1_phi", "jet1_phi", 100, -pi, pi);
    h4_jet1_btag = subdir.make<TH1D>("jet1_btag", "jet1_btag", 100, 0, 1);

    h4_jet2_m   = subdir.make<TH1D>("jet2_m", "jet2_m", 500, 0, 500);
    h4_jet2_pt  = subdir.make<TH1D>("jet2_pt", "jet2_pt", 1000, 0, 1000);
    h4_jet2_eta = subdir.make<TH1D>("jet2_eta", "jet2_eta", 100, -maxeta, maxeta);
    h4_jet2_phi = subdir.make<TH1D>("jet2_phi", "jet2_phi", 100, -pi, pi);
    h4_jet2_btag = subdir.make<TH1D>("jet2_btag", "jet2_btag", 100, 0, 1);

    h4_jet3_m   = subdir.make<TH1D>("jet3_m", "jet3_m", 500, 0, 500);
    h4_jet3_pt  = subdir.make<TH1D>("jet3_pt", "jet3_pt", 1000, 0, 1000);
    h4_jet3_eta = subdir.make<TH1D>("jet3_eta", "jet3_eta", 100, -maxeta, maxeta);
    h4_jet3_phi = subdir.make<TH1D>("jet3_phi", "jet3_phi", 100, -pi, pi);
    h4_jet3_btag = subdir.make<TH1D>("jet3_btag", "jet3_btag", 100, 0, 1);

    h4_jet4_m   = subdir.make<TH1D>("jet4_m", "jet4_m", 500, 0, 500);
    h4_jet4_pt  = subdir.make<TH1D>("jet4_pt", "jet4_pt", 1000, 0, 1000);
    h4_jet4_eta = subdir.make<TH1D>("jet4_eta", "jet4_eta", 100, -maxeta, maxeta);
    h4_jet4_phi = subdir.make<TH1D>("jet4_phi", "jet4_phi", 100, -pi, pi);
    h4_jet4_btag = subdir.make<TH1D>("jet4_btag", "jet4_btag", 100, 0, 1);

    h4_jet5_m   = subdir.make<TH1D>("jet5_m", "jet5_m", 500, 0, 500);
    h4_jet5_pt  = subdir.make<TH1D>("jet5_pt", "jet5_pt", 1000, 0, 1000);
    h4_jet5_eta = subdir.make<TH1D>("jet5_eta", "jet5_eta", 100, -maxeta, maxeta);
    h4_jet5_phi = subdir.make<TH1D>("jet5_phi", "jet5_phi", 100, -pi, pi);
    h4_jet5_btag = subdir.make<TH1D>("jet5_btag", "jet5_btag", 100, 0, 1);

    h4_jet6_m   = subdir.make<TH1D>("jet6_m", "jet6_m", 500, 0, 500);
    h4_jet6_pt  = subdir.make<TH1D>("jet6_pt", "jet6_pt", 1000, 0, 1000);
    h4_jet6_eta = subdir.make<TH1D>("jet6_eta", "jet6_eta", 100, -maxeta, maxeta);
    h4_jet6_phi = subdir.make<TH1D>("jet6_phi", "jet6_phi", 100, -pi, pi);
    h4_jet6_btag = subdir.make<TH1D>("jet6_btag", "jet6_btag", 100, 0, 1);

    h4_bjets_n = subdir.make<TH1D>("bjets_n", "bjets_n", 10, 0, 10);

    h4_event_st = subdir.make<TH1D>("event_st", "event_st", 1000, 0, 1000);
  };
};
const int ControlPlotsTTLJ::nMaxCutstep = 8; // 5+3

class TTLJEventSelector : public edm::one::EDFilter<edm::one::SharedResources>
{
public:
  TTLJEventSelector(const edm::ParameterSet& pset);
  bool filter(edm::Event& event, const edm::EventSetup&) override;
  ~TTLJEventSelector();

private:
  typedef std::vector<float> vfloat;
  typedef std::vector<double> vdouble;
  edm::EDGetTokenT<float> pileupWeightToken_, genWeightToken_;
  edm::EDGetTokenT<vfloat> genWeightsToken_;
  int genWeightIndex_;

  edm::EDGetTokenT<cat::MuonCollection> muonToken_;
  edm::EDGetTokenT<cat::ElectronCollection> electronToken_;
  edm::EDGetTokenT<cat::JetCollection> jetToken_;
  edm::EDGetTokenT<cat::METCollection> metToken_;

  edm::EDGetTokenT<int> recoFilterToken_;
  edm::EDGetTokenT<int> trigElToken_, trigMuToken_;
  edm::EDGetTokenT<int> nVertexToken_;

  std::vector<edm::EDGetTokenT<float> > extWeightTokensF_;
  std::vector<edm::EDGetTokenT<double> > extWeightTokensD_;

private:
  TH1D* h_weight, * h_pileupWeight, * h_genWeight;
  ControlPlotsTTLJ h_el, h_mu;

private:
  double shiftedMuonPt(const cat::Muon& mu) { return mu.pt()+muonScale_*mu.shiftedEn(); }
  double shiftedElectronPt(const cat::Electron& el) { return el.pt()+electronScale_*el.shiftedEn(); }
  double shiftedLepPt(const reco::Candidate& cand)
  {
    auto muonP = dynamic_cast<const cat::Muon*>(&cand);
    auto electronP = dynamic_cast<const cat::Electron*>(&cand);
    if ( muonP ) return shiftedMuonPt(*muonP);
    else if ( electronP ) return shiftedElectronPt(*electronP);
    return cand.pt();
  }
  double shiftedJetPt(const reco::Candidate& cand)
  {
    const auto jet = dynamic_cast<const cat::Jet&>(cand);
    double pt = jet.pt();
    if      ( jetScale_ == +1 ) pt *= jet.shiftedEnUp();
    else if ( jetScale_ == -1 ) pt *= jet.shiftedEnDown();

    if ( isMC_ and !isSkipJER_ ) pt *= jet.smearedRes(jetResol_);

    return pt;
  }

  bool isGoodMuon(const cat::Muon& mu)
  {
    if ( std::abs(mu.eta()) > 2.1 ) return false;
    if ( shiftedMuonPt(mu) < 26 ) return false;

    if ( mu.relIso(0.4) > 0.15 ) return false;
    if ( !mu.isTightMuon() ) return false;
    return true;
  }
  bool isGoodElectron(const cat::Electron& el)
  {
    if ( std::abs(el.eta()) > 2.4 ) return false;
    if ( shiftedElectronPt(el) < 30 ) return false;

    if ( isMVAElectronSel_ and !el.isTrigMVAValid() ) return false;

    //if ( el.relIso(0.3) >= 0.11 ) return false;
    if ( !el.electronID(elIdName_) ) return false;
    //if ( !el.isPF() or !el.passConversionVeto() ) return false;
    const double scEta = std::abs(el.scEta());
    if ( isEcalCrackVeto_ and scEta > 1.4442 and scEta < 1.566 ) return false;
    return true;
  }
  bool isVetoMuon(const cat::Muon& mu)
  {
    if ( std::abs(mu.eta()) > 2.4 ) return false;
    if ( mu.pt() < 10 ) return false;

    if ( !mu.isLooseMuon() ) return false;
    if ( mu.relIso(0.4) > 0.25 ) return false;
    return true;
  }
  bool isVetoElectron(const cat::Electron& el)
  {
    if ( !el.electronID(elVetoIdName_) ) return false;
    return true;
  }
  bool isBjet(const cat::Jet& jet)
  {
    const double bTag = jet.bDiscriminator(bTagName_);
    if      ( bTagWP_ == BTagWP::CSVL ) return bTag > WP_BTAG_CSVv2L;
    else if ( bTagWP_ == BTagWP::CSVM ) return bTag > WP_BTAG_CSVv2M;
    else if ( bTagWP_ == BTagWP::CSVT ) return bTag > WP_BTAG_CSVv2T;
    return false;
  }

private:
  typedef reco::Candidate::LorentzVector LV;

  // Energy scales
  int muonScale_, electronScale_, jetScale_, jetResol_;
  bool isSkipJER_; // Do not apply JER, needed to remove randomness during the Synchronization

  // Efficiency SF
  ScaleFactorEvaluator muonSF_, electronSF_;
  double muonSFShift_, electronSFShift_;
  double trigSFShift_;

  bool isMC_;
  bool isIgnoreTrig_; // Accept event even if it does not pass HLT. Needed for synchronization
  const int applyFilterAt_;

  // ID variables
  bool isEcalCrackVeto_, isMVAElectronSel_;
  std::string bTagName_;
  std::string elIdName_, elVetoIdName_;
  enum class BTagWP { CSVL, CSVM, CSVT } bTagWP_;

};

}

using namespace cat;

TTLJEventSelector::TTLJEventSelector(const edm::ParameterSet& pset):
  isMC_(pset.getParameter<bool>("isMC")),
  applyFilterAt_(pset.getParameter<int>("applyFilterAt"))
{
  const auto muonSet = pset.getParameter<edm::ParameterSet>("muon");
  muonToken_ = consumes<cat::MuonCollection>(muonSet.getParameter<edm::InputTag>("src"));
  muonScale_ = muonSet.getParameter<int>("scaleDirection");
  if ( isMC_ ) {
    const auto muonSFSet = muonSet.getParameter<edm::ParameterSet>("efficiencySF");
    // FIXME : for muons, eta bins are folded - always double check this with cfg
    muonSF_.set(muonSFSet.getParameter<vdouble>("pt_bins"),
                muonSFSet.getParameter<vdouble>("abseta_bins"),
                muonSFSet.getParameter<vdouble>("values"),
                muonSFSet.getParameter<vdouble>("errors"));
    muonSFShift_ = muonSet.getParameter<int>("efficiencySFDirection");
  }

  const auto electronSet = pset.getParameter<edm::ParameterSet>("electron");
  electronToken_ = consumes<cat::ElectronCollection>(electronSet.getParameter<edm::InputTag>("src"));
  elIdName_ = electronSet.getParameter<string>("idName");
  elVetoIdName_ = electronSet.getParameter<string>("vetoIdName");
  electronScale_ = electronSet.getParameter<int>("scaleDirection");
  if ( isMC_ ) {
    const auto electronSFSet = electronSet.getParameter<edm::ParameterSet>("efficiencySF");
    // FIXME : for electrons, eta bins are NOT folded - always double check this with cfg
    electronSF_.set(electronSFSet.getParameter<vdouble>("pt_bins"),
                    electronSFSet.getParameter<vdouble>("abseta_bins"),
                    electronSFSet.getParameter<vdouble>("values"),
                    electronSFSet.getParameter<vdouble>("errors"));
    electronSFShift_ = electronSet.getParameter<int>("efficiencySFDirection");
  }
  isEcalCrackVeto_ = isMVAElectronSel_ = false;
  if ( elIdName_.substr(0,3) == "mva" ) {
    isMVAElectronSel_ = true;
  }
  else {
    isEcalCrackVeto_ = electronSet.getParameter<bool>("applyEcalCrackVeto");
  }

  const auto jetSet = pset.getParameter<edm::ParameterSet>("jet");
  jetToken_ = consumes<cat::JetCollection>(jetSet.getParameter<edm::InputTag>("src"));
  jetScale_ = jetSet.getParameter<int>("scaleDirection");
  jetResol_ = jetSet.getParameter<int>("resolDirection");
  bTagName_ = jetSet.getParameter<string>("bTagName");
  const auto bTagWPStr = jetSet.getParameter<string>("bTagWP");
  if      ( bTagWPStr == "CSVL" ) bTagWP_ = BTagWP::CSVL;
  else if ( bTagWPStr == "CSVM" ) bTagWP_ = BTagWP::CSVM;
  else if ( bTagWPStr == "CSVT" ) bTagWP_ = BTagWP::CSVT;
  else edm::LogError("TTLJEventSelector") << "Wrong bTagWP parameter " << bTagWPStr;
  isSkipJER_ = jetSet.getParameter<bool>("skipJER");

  const auto metSet = pset.getParameter<edm::ParameterSet>("met");
  metToken_ = consumes<cat::METCollection>(metSet.getParameter<edm::InputTag>("src"));

  const auto vertexSet = pset.getParameter<edm::ParameterSet>("vertex");
  nVertexToken_ = consumes<int>(vertexSet.getParameter<edm::InputTag>("nVertex"));
  pileupWeightToken_ = consumes<float>(vertexSet.getParameter<edm::InputTag>("pileupWeight"));

  const auto filterSet = pset.getParameter<edm::ParameterSet>("filters");
  recoFilterToken_ = consumes<int>(filterSet.getParameter<edm::InputTag>("filterRECO"));
  trigElToken_ = consumes<int>(filterSet.getParameter<edm::InputTag>("trigEL"));
  trigMuToken_ = consumes<int>(filterSet.getParameter<edm::InputTag>("trigMU"));
  isIgnoreTrig_ = filterSet.getParameter<bool>("ignoreTrig");
  trigSFShift_ = filterSet.getParameter<int>("efficiencySFDirection");

  if ( isMC_ ) {
    const auto genWeightSet = pset.getParameter<edm::ParameterSet>("genWeight");

    genWeightIndex_ = genWeightSet.getParameter<int>("index");
    if ( genWeightIndex_ < 0 ) genWeightToken_ = consumes<float>(genWeightSet.getParameter<edm::InputTag>("src"));
    else genWeightsToken_ = consumes<vfloat>(genWeightSet.getParameter<edm::InputTag>("src"));
  }

  // Other weights
  const auto extWeightLabels = pset.getParameter<std::vector<edm::InputTag> >("extWeights");
  for ( auto x : extWeightLabels ) {
    extWeightTokensF_.push_back(consumes<float>(x));
    extWeightTokensD_.push_back(consumes<double>(x));
  }

  // Fill histograms, etc
  usesResource("TFileService");
  edm::Service<TFileService> fs;

  auto doverall = fs->mkdir("overall", "overall");
  h_weight = doverall.make<TH1D>("weight", "weight", 200, -10, 10);
  if ( isMC_ ) {
    h_genWeight = doverall.make<TH1D>("genWeight", "genWeight", 200, -10, 10);
    h_pileupWeight = doverall.make<TH1D>("pileupWeight", "pileupWeight", 200, -10, 10);
  }

  h_el.book(fs->mkdir("el"));
  h_mu.book(fs->mkdir("mu"));

  produces<int>("channel");
  produces<float>("weight");
  produces<float>("met");
  produces<float>("metphi");
  produces<std::vector<cat::Lepton> >("leptons");
  produces<std::vector<cat::Jet> >("jets");
}

bool TTLJEventSelector::filter(edm::Event& event, const edm::EventSetup&)
{
  if ( event.isRealData() ) isMC_ = false;

  // Get physics objects
  edm::Handle<cat::MuonCollection> muonHandle;
  event.getByToken(muonToken_, muonHandle);

  edm::Handle<cat::ElectronCollection> electronHandle;
  event.getByToken(electronToken_, electronHandle);

  edm::Handle<cat::JetCollection> jetHandle;
  event.getByToken(jetToken_, jetHandle);

  edm::Handle<cat::METCollection> metHandle;
  event.getByToken(metToken_, metHandle);
  const auto& metP4 = metHandle->at(0).p4();
  double metDpx = 0, metDpy = 0;

  edm::Handle<int> nVertexHandle;
  event.getByToken(nVertexToken_, nVertexHandle);
  const int nVertex = *nVertexHandle;

  std::auto_ptr<std::vector<cat::Lepton> > out_leptons(new std::vector<cat::Lepton>());
  std::auto_ptr<std::vector<cat::Jet> > out_jets(new std::vector<cat::Jet>());

  // Compute event weight - from generator, pileup, etc
  double weight = 1.0;
  if ( isMC_ ) {
    float genWeight = 1.;
    edm::Handle<float> fHandle;
    edm::Handle<vfloat> vfHandle;

    if ( genWeightIndex_ < 0 ) {
      event.getByToken(genWeightToken_, fHandle);
      genWeight = *fHandle;
    }
    else {
      event.getByToken(genWeightsToken_, vfHandle);
      genWeight = vfHandle->at(genWeightIndex_);
    }

    event.getByToken(pileupWeightToken_, fHandle);
    const float pileupWeight = *fHandle;

    h_genWeight->Fill(genWeight);
    h_pileupWeight->Fill(pileupWeight);
    weight *= genWeight*pileupWeight;
    // NOTE: weight value to be multiplied by lepton SF, etc.
  }

  // Apply all other weights
  for ( auto t : extWeightTokensF_ ) {
    edm::Handle<float> h;
    if ( event.getByToken(t, h) ) weight *= *h;
  }
  for ( auto t : extWeightTokensD_ ) {
    edm::Handle<double> h;
    if ( event.getByToken(t, h) ) weight *= *h;
  }

  // Get event filters and triggers
  edm::Handle<int> trigHandle;
  event.getByToken(recoFilterToken_, trigHandle);
  const int isRECOFilterOK = *trigHandle;

  event.getByToken(trigElToken_, trigHandle);
  const int isTrigEl = *trigHandle;
  event.getByToken(trigMuToken_, trigHandle);
  const int isTrigMu = *trigHandle;

  // Select good leptons
  double leptons_st = 0;
  cat::MuonCollection selMuons, vetoMuons;
  for ( int i=0, n=muonHandle->size(); i<n; ++i ) {
    auto& p = muonHandle->at(i);
    const double pt = shiftedMuonPt(p);
    const double scale = pt/p.pt();

    cat::Muon lep(p);
    lep.setP4(p.p4()*scale);
    const bool isGood = isGoodMuon(lep);
    const bool isVeto = isVetoMuon(lep);
    if ( isGood ) selMuons.push_back(lep);
    else if ( isVeto ) vetoMuons.push_back(lep);

    if ( isGood or isVeto ) {
      leptons_st += pt;
      metDpx += lep.px()-p.px();
      metDpy += lep.py()-p.py();
    }
  }
  cat::ElectronCollection selElectrons, vetoElectrons;
  for ( int i=0, n=electronHandle->size(); i<n; ++i ) {
    auto& p = electronHandle->at(i);
    const double pt = shiftedElectronPt(p);
    const double scale = pt/p.pt();

    cat::Electron lep(p);
    lep.setP4(p.p4()*scale);
    selElectrons.push_back(lep);
    const bool isGood = isGoodElectron(lep);
    const bool isVeto = isVetoElectron(lep);
    if ( isGood ) selElectrons.push_back(lep);
    else if ( isVeto ) vetoElectrons.push_back(lep);

    if ( isGood or isVeto ) {
      leptons_st += pt;
      metDpx += lep.px()-p.px();
      metDpy += lep.py()-p.py();
    }
  }
  std::vector<const cat::Lepton*> selLeptons;
  for ( auto& x : selMuons ) selLeptons.push_back(&x);
  for ( auto& x : selElectrons ) selLeptons.push_back(&x);
  std::sort(selLeptons.begin(), selLeptons.end(),
            [&](const cat::Lepton* a, const cat::Lepton* b){return a->pt() > b->pt();});
  // Copy selLeptons to out_leptons
  if ( !selLeptons.empty() ) out_leptons->push_back(*selLeptons.at(0));
  const int leptons_n = selLeptons.size();
  const cat::Lepton* lepton1 = 0;
  int channel = 0;
  if ( leptons_n >= 1 ) {
    // Set lepton1
    lepton1 = selLeptons.at(0);

    const int pdgId1 = std::abs(lepton1->pdgId());
    // Determine channel
    channel = abs(pdgId1);

    // Apply lepton SF
    if ( channel == 11 ) {
      const auto e1 = dynamic_cast<const cat::Electron*>(lepton1);
      const double w1 = electronSF_(lepton1->pt(), std::abs(e1->scEta()), electronSFShift_);
      weight *= w1;
      if ( !isIgnoreTrig_ ) weight *= isTrigEl;// * computeTrigSF(*lepton1, trigSFShift_);
    }
    else if ( channel == 13 ) {
      const double w1 = muonSF_(lepton1->pt(), std::abs(lepton1->eta()), muonSFShift_);
      weight *= w1;
      if ( !isIgnoreTrig_ ) weight *= isTrigMu;// * computeTrigSF(*lepton1, trigSFShift_);
    }
    else edm::LogError("TTLJEventSelector") << "Strange event with nLepton >=2 but not falling info ee,mumu,emu category";
  }

  // Select good jets
  int bjets_n = 0;
  double jets_ht = 0;
  for ( int i=0, n=jetHandle->size(); i<n; ++i ) {
    auto& p = jetHandle->at(i);
    if ( std::abs(p.eta()) > 2.4 ) continue;
    if ( !p.LooseId() ) continue;

    const double pt = shiftedJetPt(p);
    const double scale = pt/p.pt();
    cat::Jet jet(p);
    jet.setP4(scale*p.p4());

    metDpx += jet.px()-p.px();
    metDpy += jet.py()-p.py();
    if ( pt < 30 ) continue;

    if ( leptons_n >= 1 and deltaR(jet.p4(), out_leptons->at(0).p4()) < 0.4 ) continue;

    out_jets->push_back(jet);
    jets_ht += pt;
    if ( isBjet(p) ) ++bjets_n;
  }
  const int jets_n = out_jets->size();
  std::sort(out_jets->begin(), out_jets->end(),
            [&](const cat::Jet& a, const cat::Jet& b){return a.pt() > b.pt();});

  // Update & calculate met
  const double met_pt = hypot(metP4.px()-metDpx, metP4.py()-metDpy);
  const double met_phi = atan2(metP4.px()-metDpx, metP4.py()-metDpy);

  // Check cut steps and fill histograms
  h_weight->Fill(weight);

  h_el.hCutstep->Fill(-2, weight);
  h_el.hCutstepNoweight->Fill(-2);
  h_el.h0a_vertex_n->Fill(nVertex, weight);

  h_mu.hCutstep->Fill(-2, weight);
  h_mu.hCutstepNoweight->Fill(-2);
  h_mu.h0a_vertex_n->Fill(nVertex, weight);

  // El channel Cutstep 0b with trigger requirements
  int cutstep_el = -2;
  if ( isIgnoreTrig_ or isTrigEl ) {
    ++cutstep_el;
    h_el.hCutstep->Fill(-1, weight);
    h_el.hCutstepNoweight->Fill(-1);
    h_el.h0b_vertex_n->Fill(nVertex, weight);
    h_el.h0b_met_pt->Fill(met_pt, weight);
    h_el.h0b_met_phi->Fill(met_phi, weight);
    h_el.h0b_leptons_n->Fill(leptons_n, weight);
    h_el.h0b_jets_n->Fill(jets_n, weight);
    h_el.h0b_bjets_n->Fill(bjets_n, weight);
    h_el.h0b_jets_ht->Fill(jets_ht, weight);
    for ( auto jet : *out_jets ) {
      h_el.h0b_jets_pt->Fill(jet.pt(), weight);
      h_el.h0b_jets_eta->Fill(jet.eta(), weight);
    }

    // Cutstep 0c with reco filters
    if ( isMC_ or isRECOFilterOK ) {
      ++cutstep_el;
      h_el.hCutstep->Fill(0., weight);
      h_el.hCutstepNoweight->Fill(0.);
      h_el.h0c_vertex_n->Fill(nVertex, weight);
      h_el.h0c_met_pt->Fill(met_pt, weight);
      h_el.h0c_met_phi->Fill(met_phi, weight);
      h_el.h0c_leptons_n->Fill(leptons_n, weight);
      h_el.h0c_jets_n->Fill(jets_n, weight);
      h_el.h0c_bjets_n->Fill(bjets_n, weight);
      h_el.h0c_jets_ht->Fill(jets_ht, weight);
      for ( auto jet : *out_jets ) {
        h_el.h0c_jets_pt->Fill(jet.pt(), weight);
        h_el.h0c_jets_eta->Fill(jet.eta(), weight);
      }
    }
  }
  // Mu channel Cutstep 0b with trigger requirements
  int cutstep_mu = -2;
  if ( isIgnoreTrig_ or isTrigMu ) {
    ++cutstep_mu;
    h_mu.hCutstep->Fill(-1, weight);
    h_mu.hCutstepNoweight->Fill(-1);
    h_mu.h0b_vertex_n->Fill(nVertex, weight);
    h_mu.h0b_met_pt->Fill(met_pt, weight);
    h_mu.h0b_met_phi->Fill(met_phi, weight);
    h_mu.h0b_leptons_n->Fill(leptons_n, weight);
    h_mu.h0b_jets_n->Fill(jets_n, weight);
    h_mu.h0b_bjets_n->Fill(bjets_n, weight);
    h_mu.h0b_jets_ht->Fill(jets_ht, weight);
    for ( auto jet : *out_jets ) {
      h_mu.h0b_jets_pt->Fill(jet.pt(), weight);
      h_mu.h0b_jets_eta->Fill(jet.eta(), weight);
    }

    // Cutstep 0c with reco filters
    if ( isMC_ or isRECOFilterOK ) {
      ++cutstep_mu;
      h_mu.hCutstep->Fill(0., weight);
      h_mu.hCutstepNoweight->Fill(0.);
      h_mu.h0c_vertex_n->Fill(nVertex, weight);
      h_mu.h0c_met_pt->Fill(met_pt, weight);
      h_mu.h0c_met_phi->Fill(met_phi, weight);
      h_mu.h0c_leptons_n->Fill(leptons_n, weight);
      h_mu.h0c_jets_n->Fill(jets_n, weight);
      h_mu.h0c_bjets_n->Fill(bjets_n, weight);
      h_mu.h0c_jets_ht->Fill(jets_ht, weight);
      for ( auto jet : *out_jets ) {
        h_mu.h0c_jets_pt->Fill(jet.pt(), weight);
        h_mu.h0c_jets_eta->Fill(jet.eta(), weight);
      }
    }
  }

  // Check each cut steps
  int cutstep = -1;
  // bitset for the cut steps, fill the results only for events that pass step0a,0b,0c
  std::bitset<ControlPlotsTTLJ::nMaxCutstep-2> cutstepBits(0);
  //for ( auto x : cutstepBits ) x = false;
  if ( (channel == 11 and cutstep_el == 0) or
       (channel == 13 and cutstep_mu == 0) ) {
    // Step1 at least one signal lepton
    if ( leptons_n > 0 ) cutstepBits[0] = true;

    // Step 2a veto any additional muon
    if ( vetoMuons.empty() and (
         (channel == 11 and selMuons.empty()) or
         (channel == 13 and selMuons.size() == 1)) ) cutstepBits[1] = true;

    // Step 2b veto any additional electron
    if ( vetoElectrons.empty() and (
         (channel == 11 and selElectrons.size() == 1) or
         (channel == 13 and selElectrons.empty())) ) cutstepBits[2] = true;

    // Step3a Minimal jet multiplicity
    if ( jets_n >= 1 ) cutstepBits[3] = true;
    // Step4 one b jet
    if ( bjets_n >= 1 ) cutstepBits[4] = true;

    // Set the cut step of this event
    const int nstep = cutstepBits.size();
    for ( cutstep=0; cutstep<nstep; ++cutstep ) {
      if ( !cutstepBits[cutstep] ) break;
    }
  }
  else
  {
    cutstep = std::max(cutstep_el, cutstep_mu); // reset the cut step
  }

  // Cut step is ready. Now proceed to fill histograms
  switch (1) default: {
    int icutstep = 0;

    if ( cutstep <= 0 ) break;
    ++icutstep; // =1

    const auto lepton1P4 = shiftedLepPt(*lepton1)/lepton1->pt()*lepton1->p4();

    if ( channel == 11 ) {
      h_el.hCutstep->Fill(icutstep, weight);
      h_el.hCutstepNoweight->Fill(icutstep);
      h_el.h1_vertex_n->Fill(nVertex, weight);
      h_el.h1_met_pt->Fill(met_pt, weight);
      h_el.h1_met_phi->Fill(met_phi, weight);
      h_el.h1_leptons_n->Fill(leptons_n, weight);
      h_el.h1_lepton1_pt->Fill(lepton1P4.pt(), weight);
      h_el.h1_lepton1_eta->Fill(lepton1->eta(), weight);
      h_el.h1_lepton1_phi->Fill(lepton1->phi(), weight);
      h_el.h1_lepton1_q->Fill(lepton1->charge(), weight);
      h_el.h1_jets_n->Fill(jets_n, weight);
      h_el.h1_jets_ht->Fill(jets_ht, weight);
      for ( auto jet : *out_jets ) {
        h_el.h1_jets_pt->Fill(jet.pt(), weight);
        h_el.h1_jets_eta->Fill(jet.eta(), weight);
      }
      if ( jets_n >= 1 ) {
        const auto& jet = out_jets->at(0);
        h_el.h1_jet1_m->Fill(jet.mass(), weight);
        h_el.h1_jet1_pt->Fill(jet.pt(), weight);
        h_el.h1_jet1_eta->Fill(jet.eta(), weight);
        h_el.h1_jet1_phi->Fill(jet.phi(), weight);
        h_el.h1_jet1_btag->Fill(jet.bDiscriminator(bTagName_), weight);
      }
      if ( jets_n >= 2 ) {
        const auto& jet = out_jets->at(1);
        h_el.h1_jet2_m->Fill(jet.mass(), weight);
        h_el.h1_jet2_pt->Fill(jet.pt(), weight);
        h_el.h1_jet2_eta->Fill(jet.eta(), weight);
        h_el.h1_jet2_phi->Fill(jet.phi(), weight);
        h_el.h1_jet2_btag->Fill(jet.bDiscriminator(bTagName_), weight);
      }
      if ( jets_n >= 3 ) {
        const auto& jet = out_jets->at(2);
        h_el.h1_jet3_m->Fill(jet.mass(), weight);
        h_el.h1_jet3_pt->Fill(jet.pt(), weight);
        h_el.h1_jet3_eta->Fill(jet.eta(), weight);
        h_el.h1_jet3_phi->Fill(jet.phi(), weight);
        h_el.h1_jet3_btag->Fill(jet.bDiscriminator(bTagName_), weight);
      }
      if ( jets_n >= 4 ) {
        const auto& jet = out_jets->at(3);
        h_el.h1_jet4_m->Fill(jet.mass(), weight);
        h_el.h1_jet4_pt->Fill(jet.pt(), weight);
        h_el.h1_jet4_eta->Fill(jet.eta(), weight);
        h_el.h1_jet4_phi->Fill(jet.phi(), weight);
        h_el.h1_jet4_btag->Fill(jet.bDiscriminator(bTagName_), weight);
      }
      if ( jets_n >= 5 ) {
        const auto& jet = out_jets->at(4);
        h_el.h1_jet5_m->Fill(jet.mass(), weight);
        h_el.h1_jet5_pt->Fill(jet.pt(), weight);
        h_el.h1_jet5_eta->Fill(jet.eta(), weight);
        h_el.h1_jet5_phi->Fill(jet.phi(), weight);
        h_el.h1_jet5_btag->Fill(jet.bDiscriminator(bTagName_), weight);
      }
      if ( jets_n >= 6 ) {
        const auto& jet = out_jets->at(5);
        h_el.h1_jet6_m->Fill(jet.mass(), weight);
        h_el.h1_jet6_pt->Fill(jet.pt(), weight);
        h_el.h1_jet6_eta->Fill(jet.eta(), weight);
        h_el.h1_jet6_phi->Fill(jet.phi(), weight);
        h_el.h1_jet6_btag->Fill(jet.bDiscriminator(bTagName_), weight);
      }
      h_el.h1_bjets_n->Fill(bjets_n, weight);

      h_el.h1_event_st->Fill(leptons_st+jets_ht+met_pt, weight);
    }
    else if ( channel == 13 ) {
      h_mu.hCutstep->Fill(icutstep, weight);
      h_mu.hCutstepNoweight->Fill(icutstep);
      h_mu.h1_vertex_n->Fill(nVertex, weight);
      h_mu.h1_met_pt->Fill(met_pt, weight);
      h_mu.h1_met_phi->Fill(met_phi, weight);
      h_mu.h1_leptons_n->Fill(leptons_n, weight);
      h_mu.h1_lepton1_pt->Fill(lepton1P4.pt(), weight);
      h_mu.h1_lepton1_eta->Fill(lepton1->eta(), weight);
      h_mu.h1_lepton1_phi->Fill(lepton1->phi(), weight);
      h_mu.h1_lepton1_q->Fill(lepton1->charge(), weight);
      h_mu.h1_jets_n->Fill(jets_n, weight);
      h_mu.h1_jets_ht->Fill(jets_ht, weight);
      for ( auto jet : *out_jets ) {
        h_mu.h1_jets_pt->Fill(jet.pt(), weight);
        h_mu.h1_jets_eta->Fill(jet.eta(), weight);
      }
      if ( jets_n >= 1 ) {
        const auto& jet = out_jets->at(0);
        h_mu.h1_jet1_m->Fill(jet.mass(), weight);
        h_mu.h1_jet1_pt->Fill(jet.pt(), weight);
        h_mu.h1_jet1_eta->Fill(jet.eta(), weight);
        h_mu.h1_jet1_phi->Fill(jet.phi(), weight);
        h_mu.h1_jet1_btag->Fill(jet.bDiscriminator(bTagName_), weight);
      }
      if ( jets_n >= 2 ) {
        const auto& jet = out_jets->at(1);
        h_mu.h1_jet2_m->Fill(jet.mass(), weight);
        h_mu.h1_jet2_pt->Fill(jet.pt(), weight);
        h_mu.h1_jet2_eta->Fill(jet.eta(), weight);
        h_mu.h1_jet2_phi->Fill(jet.phi(), weight);
        h_mu.h1_jet2_btag->Fill(jet.bDiscriminator(bTagName_), weight);
      }
      if ( jets_n >= 3 ) {
        const auto& jet = out_jets->at(2);
        h_mu.h1_jet3_m->Fill(jet.mass(), weight);
        h_mu.h1_jet3_pt->Fill(jet.pt(), weight);
        h_mu.h1_jet3_eta->Fill(jet.eta(), weight);
        h_mu.h1_jet3_phi->Fill(jet.phi(), weight);
        h_mu.h1_jet3_btag->Fill(jet.bDiscriminator(bTagName_), weight);
      }
      if ( jets_n >= 4 ) {
        const auto& jet = out_jets->at(3);
        h_mu.h1_jet4_m->Fill(jet.mass(), weight);
        h_mu.h1_jet4_pt->Fill(jet.pt(), weight);
        h_mu.h1_jet4_eta->Fill(jet.eta(), weight);
        h_mu.h1_jet4_phi->Fill(jet.phi(), weight);
        h_mu.h1_jet4_btag->Fill(jet.bDiscriminator(bTagName_), weight);
      }
      if ( jets_n >= 5 ) {
        const auto& jet = out_jets->at(4);
        h_mu.h1_jet5_m->Fill(jet.mass(), weight);
        h_mu.h1_jet5_pt->Fill(jet.pt(), weight);
        h_mu.h1_jet5_eta->Fill(jet.eta(), weight);
        h_mu.h1_jet5_phi->Fill(jet.phi(), weight);
        h_mu.h1_jet5_btag->Fill(jet.bDiscriminator(bTagName_), weight);
      }
      if ( jets_n >= 6 ) {
        const auto& jet = out_jets->at(5);
        h_mu.h1_jet6_m->Fill(jet.mass(), weight);
        h_mu.h1_jet6_pt->Fill(jet.pt(), weight);
        h_mu.h1_jet6_eta->Fill(jet.eta(), weight);
        h_mu.h1_jet6_phi->Fill(jet.phi(), weight);
        h_mu.h1_jet6_btag->Fill(jet.bDiscriminator(bTagName_), weight);
      }
      h_mu.h1_bjets_n->Fill(bjets_n, weight);

      h_mu.h1_event_st->Fill(leptons_st+jets_ht+met_pt, weight);
    }

    // Finalize cut step 1 and start cut step2
    if ( cutstep <= 1 ) break;
    ++icutstep; // =2

    if ( channel == 11 ) {
      h_el.hCutstep->Fill(icutstep, weight);
      h_el.hCutstepNoweight->Fill(icutstep);
      h_el.h2a_vertex_n->Fill(nVertex, weight);
      h_el.h2a_met_pt->Fill(met_pt, weight);
      h_el.h2a_met_phi->Fill(met_phi, weight);
      h_el.h2a_leptons_n->Fill(leptons_n, weight);
      h_el.h2a_lepton1_pt->Fill(lepton1P4.pt(), weight);
      h_el.h2a_lepton1_eta->Fill(lepton1->eta(), weight);
      h_el.h2a_lepton1_phi->Fill(lepton1->phi(), weight);
      h_el.h2a_lepton1_q->Fill(lepton1->charge(), weight);
      h_el.h2a_jets_n->Fill(jets_n, weight);
      h_el.h2a_jets_ht->Fill(jets_ht, weight);
      for ( auto jet : *out_jets ) {
        h_el.h2a_jets_pt->Fill(jet.pt(), weight);
        h_el.h2a_jets_eta->Fill(jet.eta(), weight);
      }
      if ( jets_n >= 1 ) {
        const auto& jet = out_jets->at(0);
        h_el.h2a_jet1_m->Fill(jet.mass(), weight);
        h_el.h2a_jet1_pt->Fill(jet.pt(), weight);
        h_el.h2a_jet1_eta->Fill(jet.eta(), weight);
        h_el.h2a_jet1_phi->Fill(jet.phi(), weight);
        h_el.h2a_jet1_btag->Fill(jet.bDiscriminator(bTagName_), weight);
      }
      if ( jets_n >= 2 ) {
        const auto& jet = out_jets->at(1);
        h_el.h2a_jet2_m->Fill(jet.mass(), weight);
        h_el.h2a_jet2_pt->Fill(jet.pt(), weight);
        h_el.h2a_jet2_eta->Fill(jet.eta(), weight);
        h_el.h2a_jet2_phi->Fill(jet.phi(), weight);
        h_el.h2a_jet2_btag->Fill(jet.bDiscriminator(bTagName_), weight);
      }
      if ( jets_n >= 3 ) {
        const auto& jet = out_jets->at(2);
        h_el.h2a_jet3_m->Fill(jet.mass(), weight);
        h_el.h2a_jet3_pt->Fill(jet.pt(), weight);
        h_el.h2a_jet3_eta->Fill(jet.eta(), weight);
        h_el.h2a_jet3_phi->Fill(jet.phi(), weight);
        h_el.h2a_jet3_btag->Fill(jet.bDiscriminator(bTagName_), weight);
      }
      if ( jets_n >= 4 ) {
        const auto& jet = out_jets->at(3);
        h_el.h2a_jet4_m->Fill(jet.mass(), weight);
        h_el.h2a_jet4_pt->Fill(jet.pt(), weight);
        h_el.h2a_jet4_eta->Fill(jet.eta(), weight);
        h_el.h2a_jet4_phi->Fill(jet.phi(), weight);
        h_el.h2a_jet4_btag->Fill(jet.bDiscriminator(bTagName_), weight);
      }
      if ( jets_n >= 5 ) {
        const auto& jet = out_jets->at(4);
        h_el.h2a_jet5_m->Fill(jet.mass(), weight);
        h_el.h2a_jet5_pt->Fill(jet.pt(), weight);
        h_el.h2a_jet5_eta->Fill(jet.eta(), weight);
        h_el.h2a_jet5_phi->Fill(jet.phi(), weight);
        h_el.h2a_jet5_btag->Fill(jet.bDiscriminator(bTagName_), weight);
      }
      if ( jets_n >= 6 ) {
        const auto& jet = out_jets->at(5);
        h_el.h2a_jet6_m->Fill(jet.mass(), weight);
        h_el.h2a_jet6_pt->Fill(jet.pt(), weight);
        h_el.h2a_jet6_eta->Fill(jet.eta(), weight);
        h_el.h2a_jet6_phi->Fill(jet.phi(), weight);
        h_el.h2a_jet6_btag->Fill(jet.bDiscriminator(bTagName_), weight);
      }
      h_el.h2a_bjets_n->Fill(bjets_n, weight);

      h_el.h2a_event_st->Fill(leptons_st+jets_ht+met_pt, weight);
    }
    else if ( channel == 13 ) {
      h_mu.hCutstep->Fill(icutstep, weight);
      h_mu.hCutstepNoweight->Fill(icutstep);
      h_mu.h2a_vertex_n->Fill(nVertex, weight);
      h_mu.h2a_met_pt->Fill(met_pt, weight);
      h_mu.h2a_met_phi->Fill(met_phi, weight);
      h_mu.h2a_leptons_n->Fill(leptons_n, weight);
      h_mu.h2a_lepton1_pt->Fill(lepton1P4.pt(), weight);
      h_mu.h2a_lepton1_eta->Fill(lepton1->eta(), weight);
      h_mu.h2a_lepton1_phi->Fill(lepton1->phi(), weight);
      h_mu.h2a_lepton1_q->Fill(lepton1->charge(), weight);
      h_mu.h2a_jets_n->Fill(jets_n, weight);
      h_mu.h2a_jets_ht->Fill(jets_ht, weight);
      for ( auto jet : *out_jets ) {
        h_mu.h2a_jets_pt->Fill(jet.pt(), weight);
        h_mu.h2a_jets_eta->Fill(jet.eta(), weight);
      }
      if ( jets_n >= 1 ) {
        const auto& jet = out_jets->at(0);
        h_mu.h2a_jet1_m->Fill(jet.mass(), weight);
        h_mu.h2a_jet1_pt->Fill(jet.pt(), weight);
        h_mu.h2a_jet1_eta->Fill(jet.eta(), weight);
        h_mu.h2a_jet1_phi->Fill(jet.phi(), weight);
        h_mu.h2a_jet1_btag->Fill(jet.bDiscriminator(bTagName_), weight);
      }
      if ( jets_n >= 2 ) {
        const auto& jet = out_jets->at(1);
        h_mu.h2a_jet2_m->Fill(jet.mass(), weight);
        h_mu.h2a_jet2_pt->Fill(jet.pt(), weight);
        h_mu.h2a_jet2_eta->Fill(jet.eta(), weight);
        h_mu.h2a_jet2_phi->Fill(jet.phi(), weight);
        h_mu.h2a_jet2_btag->Fill(jet.bDiscriminator(bTagName_), weight);
      }
      if ( jets_n >= 3 ) {
        const auto& jet = out_jets->at(2);
        h_mu.h2a_jet3_m->Fill(jet.mass(), weight);
        h_mu.h2a_jet3_pt->Fill(jet.pt(), weight);
        h_mu.h2a_jet3_eta->Fill(jet.eta(), weight);
        h_mu.h2a_jet3_phi->Fill(jet.phi(), weight);
        h_mu.h2a_jet3_btag->Fill(jet.bDiscriminator(bTagName_), weight);
      }
      if ( jets_n >= 4 ) {
        const auto& jet = out_jets->at(3);
        h_mu.h2a_jet4_m->Fill(jet.mass(), weight);
        h_mu.h2a_jet4_pt->Fill(jet.pt(), weight);
        h_mu.h2a_jet4_eta->Fill(jet.eta(), weight);
        h_mu.h2a_jet4_phi->Fill(jet.phi(), weight);
        h_mu.h2a_jet4_btag->Fill(jet.bDiscriminator(bTagName_), weight);
      }
      if ( jets_n >= 5 ) {
        const auto& jet = out_jets->at(4);
        h_mu.h2a_jet5_m->Fill(jet.mass(), weight);
        h_mu.h2a_jet5_pt->Fill(jet.pt(), weight);
        h_mu.h2a_jet5_eta->Fill(jet.eta(), weight);
        h_mu.h2a_jet5_phi->Fill(jet.phi(), weight);
        h_mu.h2a_jet5_btag->Fill(jet.bDiscriminator(bTagName_), weight);
      }
      if ( jets_n >= 6 ) {
        const auto& jet = out_jets->at(5);
        h_mu.h2a_jet6_m->Fill(jet.mass(), weight);
        h_mu.h2a_jet6_pt->Fill(jet.pt(), weight);
        h_mu.h2a_jet6_eta->Fill(jet.eta(), weight);
        h_mu.h2a_jet6_phi->Fill(jet.phi(), weight);
        h_mu.h2a_jet6_btag->Fill(jet.bDiscriminator(bTagName_), weight);
      }
      h_mu.h2a_bjets_n->Fill(bjets_n, weight);

      h_mu.h2a_event_st->Fill(leptons_st+jets_ht+met_pt, weight);
    }

    // Finalize cut step 2 and start cut step3
    if ( cutstep <= 2 ) break;
    ++icutstep; // =3

    if ( channel == 11 ) {
      h_el.hCutstep->Fill(icutstep, weight);
      h_el.hCutstepNoweight->Fill(icutstep);
      h_el.h2b_vertex_n->Fill(nVertex, weight);
      h_el.h2b_met_pt->Fill(met_pt, weight);
      h_el.h2b_met_phi->Fill(met_phi, weight);
      h_el.h2b_jets_n->Fill(jets_n, weight);
      h_el.h2b_jets_ht->Fill(jets_ht, weight);
      for ( auto jet : *out_jets ) {
        h_el.h2b_jets_pt->Fill(jet.pt(), weight);
        h_el.h2b_jets_eta->Fill(jet.eta(), weight);
      }
      if ( jets_n >= 1 ) {
        const auto& jet = out_jets->at(0);
        h_el.h2b_jet1_m->Fill(jet.mass(), weight);
        h_el.h2b_jet1_pt->Fill(jet.pt(), weight);
        h_el.h2b_jet1_eta->Fill(jet.eta(), weight);
        h_el.h2b_jet1_phi->Fill(jet.phi(), weight);
        h_el.h2b_jet1_btag->Fill(jet.bDiscriminator(bTagName_), weight);
      }
      if ( jets_n >= 2 ) {
        const auto& jet = out_jets->at(1);
        h_el.h2b_jet2_m->Fill(jet.mass(), weight);
        h_el.h2b_jet2_pt->Fill(jet.pt(), weight);
        h_el.h2b_jet2_eta->Fill(jet.eta(), weight);
        h_el.h2b_jet2_phi->Fill(jet.phi(), weight);
        h_el.h2b_jet2_btag->Fill(jet.bDiscriminator(bTagName_), weight);
      }
      if ( jets_n >= 3 ) {
        const auto& jet = out_jets->at(2);
        h_el.h2b_jet3_m->Fill(jet.mass(), weight);
        h_el.h2b_jet3_pt->Fill(jet.pt(), weight);
        h_el.h2b_jet3_eta->Fill(jet.eta(), weight);
        h_el.h2b_jet3_phi->Fill(jet.phi(), weight);
        h_el.h2b_jet3_btag->Fill(jet.bDiscriminator(bTagName_), weight);
      }
      if ( jets_n >= 4 ) {
        const auto& jet = out_jets->at(3);
        h_el.h2b_jet4_m->Fill(jet.mass(), weight);
        h_el.h2b_jet4_pt->Fill(jet.pt(), weight);
        h_el.h2b_jet4_eta->Fill(jet.eta(), weight);
        h_el.h2b_jet4_phi->Fill(jet.phi(), weight);
        h_el.h2b_jet4_btag->Fill(jet.bDiscriminator(bTagName_), weight);
      }
      if ( jets_n >= 5 ) {
        const auto& jet = out_jets->at(4);
        h_el.h2b_jet5_m->Fill(jet.mass(), weight);
        h_el.h2b_jet5_pt->Fill(jet.pt(), weight);
        h_el.h2b_jet5_eta->Fill(jet.eta(), weight);
        h_el.h2b_jet5_phi->Fill(jet.phi(), weight);
        h_el.h2b_jet5_btag->Fill(jet.bDiscriminator(bTagName_), weight);
      }
      if ( jets_n >= 6 ) {
        const auto& jet = out_jets->at(5);
        h_el.h2b_jet6_m->Fill(jet.mass(), weight);
        h_el.h2b_jet6_pt->Fill(jet.pt(), weight);
        h_el.h2b_jet6_eta->Fill(jet.eta(), weight);
        h_el.h2b_jet6_phi->Fill(jet.phi(), weight);
        h_el.h2b_jet6_btag->Fill(jet.bDiscriminator(bTagName_), weight);
      }
      h_el.h2b_bjets_n->Fill(bjets_n, weight);

      h_el.h2b_event_st->Fill(leptons_st+jets_ht+met_pt, weight);
    }
    else if ( channel == 13 ) {
      h_mu.hCutstep->Fill(icutstep, weight);
      h_mu.hCutstepNoweight->Fill(icutstep);
      h_mu.h2b_vertex_n->Fill(nVertex, weight);
      h_mu.h2b_met_pt->Fill(met_pt, weight);
      h_mu.h2b_met_phi->Fill(met_phi, weight);
      h_mu.h2b_jets_n->Fill(jets_n, weight);
      h_mu.h2b_jets_ht->Fill(jets_ht, weight);
      for ( auto jet : *out_jets ) {
        h_mu.h2b_jets_pt->Fill(jet.pt(), weight);
        h_mu.h2b_jets_eta->Fill(jet.eta(), weight);
      }
      if ( jets_n >= 1 ) {
        const auto& jet = out_jets->at(0);
        h_mu.h2b_jet1_m->Fill(jet.mass(), weight);
        h_mu.h2b_jet1_pt->Fill(jet.pt(), weight);
        h_mu.h2b_jet1_eta->Fill(jet.eta(), weight);
        h_mu.h2b_jet1_phi->Fill(jet.phi(), weight);
        h_mu.h2b_jet1_btag->Fill(jet.bDiscriminator(bTagName_), weight);
      }
      if ( jets_n >= 2 ) {
        const auto& jet = out_jets->at(1);
        h_mu.h2b_jet2_m->Fill(jet.mass(), weight);
        h_mu.h2b_jet2_pt->Fill(jet.pt(), weight);
        h_mu.h2b_jet2_eta->Fill(jet.eta(), weight);
        h_mu.h2b_jet2_phi->Fill(jet.phi(), weight);
        h_mu.h2b_jet2_btag->Fill(jet.bDiscriminator(bTagName_), weight);
      }
      if ( jets_n >= 3 ) {
        const auto& jet = out_jets->at(2);
        h_mu.h2b_jet3_m->Fill(jet.mass(), weight);
        h_mu.h2b_jet3_pt->Fill(jet.pt(), weight);
        h_mu.h2b_jet3_eta->Fill(jet.eta(), weight);
        h_mu.h2b_jet3_phi->Fill(jet.phi(), weight);
        h_mu.h2b_jet3_btag->Fill(jet.bDiscriminator(bTagName_), weight);
      }
      if ( jets_n >= 4 ) {
        const auto& jet = out_jets->at(3);
        h_mu.h2b_jet4_m->Fill(jet.mass(), weight);
        h_mu.h2b_jet4_pt->Fill(jet.pt(), weight);
        h_mu.h2b_jet4_eta->Fill(jet.eta(), weight);
        h_mu.h2b_jet4_phi->Fill(jet.phi(), weight);
        h_mu.h2b_jet4_btag->Fill(jet.bDiscriminator(bTagName_), weight);
      }
      if ( jets_n >= 5 ) {
        const auto& jet = out_jets->at(4);
        h_mu.h2b_jet5_m->Fill(jet.mass(), weight);
        h_mu.h2b_jet5_pt->Fill(jet.pt(), weight);
        h_mu.h2b_jet5_eta->Fill(jet.eta(), weight);
        h_mu.h2b_jet5_phi->Fill(jet.phi(), weight);
        h_mu.h2b_jet5_btag->Fill(jet.bDiscriminator(bTagName_), weight);
      }
      if ( jets_n >= 6 ) {
        const auto& jet = out_jets->at(5);
        h_mu.h2b_jet6_m->Fill(jet.mass(), weight);
        h_mu.h2b_jet6_pt->Fill(jet.pt(), weight);
        h_mu.h2b_jet6_eta->Fill(jet.eta(), weight);
        h_mu.h2b_jet6_phi->Fill(jet.phi(), weight);
        h_mu.h2b_jet6_btag->Fill(jet.bDiscriminator(bTagName_), weight);
      }
      h_mu.h2b_bjets_n->Fill(bjets_n, weight);

      h_mu.h2b_event_st->Fill(leptons_st+jets_ht+met_pt, weight);

    }

    // Finalize cut step 3 and start cut step 4
    if ( cutstep <= 3 ) break;
    ++icutstep; // =4

    const auto& jet1 = out_jets->at(0);

    if ( channel == 11 ) {
      h_el.hCutstep->Fill(icutstep, weight);
      h_el.hCutstepNoweight->Fill(icutstep);
      h_el.h3_vertex_n->Fill(nVertex, weight);
      h_el.h3_met_pt->Fill(met_pt, weight);
      h_el.h3_met_phi->Fill(met_phi, weight);
      h_el.h3_jets_n->Fill(jets_n, weight);
      h_el.h3_jets_ht->Fill(jets_ht, weight);
      for ( auto jet : *out_jets ) {
        h_el.h3_jets_pt->Fill(jet.pt(), weight);
        h_el.h3_jets_eta->Fill(jet.eta(), weight);
      }
      h_el.h3_jet1_m->Fill(jet1.mass(), weight);
      h_el.h3_jet1_pt->Fill(jet1.pt(), weight);
      h_el.h3_jet1_eta->Fill(jet1.eta(), weight);
      h_el.h3_jet1_phi->Fill(jet1.phi(), weight);
      h_el.h3_jet1_btag->Fill(jet1.bDiscriminator(bTagName_), weight);
      if ( jets_n >= 2 ) {
        const auto& jet = out_jets->at(1);
        h_el.h3_jet2_m->Fill(jet.mass(), weight);
        h_el.h3_jet2_pt->Fill(jet.pt(), weight);
        h_el.h3_jet2_eta->Fill(jet.eta(), weight);
        h_el.h3_jet2_phi->Fill(jet.phi(), weight);
        h_el.h3_jet2_btag->Fill(jet.bDiscriminator(bTagName_), weight);
      }
      if ( jets_n >= 3 ) {
        const auto& jet = out_jets->at(2);
        h_el.h3_jet3_m->Fill(jet.mass(), weight);
        h_el.h3_jet3_pt->Fill(jet.pt(), weight);
        h_el.h3_jet3_eta->Fill(jet.eta(), weight);
        h_el.h3_jet3_phi->Fill(jet.phi(), weight);
        h_el.h3_jet3_btag->Fill(jet.bDiscriminator(bTagName_), weight);
      }
      if ( jets_n >= 4 ) {
        const auto& jet = out_jets->at(3);
        h_el.h3_jet4_m->Fill(jet.mass(), weight);
        h_el.h3_jet4_pt->Fill(jet.pt(), weight);
        h_el.h3_jet4_eta->Fill(jet.eta(), weight);
        h_el.h3_jet4_phi->Fill(jet.phi(), weight);
        h_el.h3_jet4_btag->Fill(jet.bDiscriminator(bTagName_), weight);
      }
      if ( jets_n >= 5 ) {
        const auto& jet = out_jets->at(4);
        h_el.h3_jet5_m->Fill(jet.mass(), weight);
        h_el.h3_jet5_pt->Fill(jet.pt(), weight);
        h_el.h3_jet5_eta->Fill(jet.eta(), weight);
        h_el.h3_jet5_phi->Fill(jet.phi(), weight);
        h_el.h3_jet5_btag->Fill(jet.bDiscriminator(bTagName_), weight);
      }
      if ( jets_n >= 6 ) {
        const auto& jet = out_jets->at(5);
        h_el.h3_jet6_m->Fill(jet.mass(), weight);
        h_el.h3_jet6_pt->Fill(jet.pt(), weight);
        h_el.h3_jet6_eta->Fill(jet.eta(), weight);
        h_el.h3_jet6_phi->Fill(jet.phi(), weight);
        h_el.h3_jet6_btag->Fill(jet.bDiscriminator(bTagName_), weight);
      }
      h_el.h3_bjets_n->Fill(bjets_n, weight);

      h_el.h3_event_st->Fill(leptons_st+jets_ht+met_pt, weight);
    }
    else if ( channel == 13 ) {
      h_mu.hCutstep->Fill(icutstep, weight);
      h_mu.hCutstepNoweight->Fill(icutstep);
      h_mu.h3_vertex_n->Fill(nVertex, weight);
      h_mu.h3_met_pt->Fill(met_pt, weight);
      h_mu.h3_met_phi->Fill(met_phi, weight);
      h_mu.h3_jets_n->Fill(jets_n, weight);
      h_mu.h3_jets_ht->Fill(jets_ht, weight);
      for ( auto jet : *out_jets ) {
        h_mu.h3_jets_pt->Fill(jet.pt(), weight);
        h_mu.h3_jets_eta->Fill(jet.eta(), weight);
      }
      h_mu.h3_jet1_m->Fill(jet1.mass(), weight);
      h_mu.h3_jet1_pt->Fill(jet1.pt(), weight);
      h_mu.h3_jet1_eta->Fill(jet1.eta(), weight);
      h_mu.h3_jet1_phi->Fill(jet1.phi(), weight);
      h_mu.h3_jet1_btag->Fill(jet1.bDiscriminator(bTagName_), weight);
      if ( jets_n >= 2 ) {
        const auto& jet = out_jets->at(1);
        h_mu.h3_jet2_m->Fill(jet.mass(), weight);
        h_mu.h3_jet2_pt->Fill(jet.pt(), weight);
        h_mu.h3_jet2_eta->Fill(jet.eta(), weight);
        h_mu.h3_jet2_phi->Fill(jet.phi(), weight);
        h_mu.h3_jet2_btag->Fill(jet.bDiscriminator(bTagName_), weight);
      }
      if ( jets_n >= 3 ) {
        const auto& jet = out_jets->at(2);
        h_mu.h3_jet3_m->Fill(jet.mass(), weight);
        h_mu.h3_jet3_pt->Fill(jet.pt(), weight);
        h_mu.h3_jet3_eta->Fill(jet.eta(), weight);
        h_mu.h3_jet3_phi->Fill(jet.phi(), weight);
        h_mu.h3_jet3_btag->Fill(jet.bDiscriminator(bTagName_), weight);
      }
      if ( jets_n >= 4 ) {
        const auto& jet = out_jets->at(3);
        h_mu.h3_jet4_m->Fill(jet.mass(), weight);
        h_mu.h3_jet4_pt->Fill(jet.pt(), weight);
        h_mu.h3_jet4_eta->Fill(jet.eta(), weight);
        h_mu.h3_jet4_phi->Fill(jet.phi(), weight);
        h_mu.h3_jet4_btag->Fill(jet.bDiscriminator(bTagName_), weight);
      }
      if ( jets_n >= 5 ) {
        const auto& jet = out_jets->at(4);
        h_mu.h3_jet5_m->Fill(jet.mass(), weight);
        h_mu.h3_jet5_pt->Fill(jet.pt(), weight);
        h_mu.h3_jet5_eta->Fill(jet.eta(), weight);
        h_mu.h3_jet5_phi->Fill(jet.phi(), weight);
        h_mu.h3_jet5_btag->Fill(jet.bDiscriminator(bTagName_), weight);
      }
      if ( jets_n >= 6 ) {
        const auto& jet = out_jets->at(5);
        h_mu.h3_jet6_m->Fill(jet.mass(), weight);
        h_mu.h3_jet6_pt->Fill(jet.pt(), weight);
        h_mu.h3_jet6_eta->Fill(jet.eta(), weight);
        h_mu.h3_jet6_phi->Fill(jet.phi(), weight);
        h_mu.h3_jet6_btag->Fill(jet.bDiscriminator(bTagName_), weight);
      }
      h_mu.h3_bjets_n->Fill(bjets_n, weight);

      h_mu.h3_event_st->Fill(leptons_st+jets_ht+met_pt, weight);

    }

    // Finalize cut step 4 and start next cut step
    if ( cutstep <= 4 ) break;
    ++icutstep; // =5

    if ( channel == 11 ) {
      h_el.hCutstep->Fill(icutstep, weight);
      h_el.hCutstepNoweight->Fill(icutstep);
      h_el.h4_vertex_n->Fill(nVertex, weight);
      h_el.h4_met_pt->Fill(met_pt, weight);
      h_el.h4_met_phi->Fill(met_phi, weight);
      h_el.h4_jets_n->Fill(jets_n, weight);
      h_el.h4_jets_ht->Fill(jets_ht, weight);
      for ( auto jet : *out_jets ) {
        h_el.h4_jets_pt->Fill(jet.pt(), weight);
        h_el.h4_jets_eta->Fill(jet.eta(), weight);
      }
      h_el.h4_jet1_m->Fill(jet1.mass(), weight);
      h_el.h4_jet1_pt->Fill(jet1.pt(), weight);
      h_el.h4_jet1_eta->Fill(jet1.eta(), weight);
      h_el.h4_jet1_phi->Fill(jet1.phi(), weight);
      h_el.h4_jet1_btag->Fill(jet1.bDiscriminator(bTagName_), weight);
      if ( jets_n >= 2 ) {
        const auto& jet = out_jets->at(1);
        h_el.h4_jet2_m->Fill(jet.mass(), weight);
        h_el.h4_jet2_pt->Fill(jet.pt(), weight);
        h_el.h4_jet2_eta->Fill(jet.eta(), weight);
        h_el.h4_jet2_phi->Fill(jet.phi(), weight);
        h_el.h4_jet2_btag->Fill(jet.bDiscriminator(bTagName_), weight);
      }
      if ( jets_n >= 3 ) {
        const auto& jet = out_jets->at(2);
        h_el.h4_jet3_m->Fill(jet.mass(), weight);
        h_el.h4_jet3_pt->Fill(jet.pt(), weight);
        h_el.h4_jet3_eta->Fill(jet.eta(), weight);
        h_el.h4_jet3_phi->Fill(jet.phi(), weight);
        h_el.h4_jet3_btag->Fill(jet.bDiscriminator(bTagName_), weight);
      }
      if ( jets_n >= 4 ) {
        const auto& jet = out_jets->at(3);
        h_el.h4_jet4_m->Fill(jet.mass(), weight);
        h_el.h4_jet4_pt->Fill(jet.pt(), weight);
        h_el.h4_jet4_eta->Fill(jet.eta(), weight);
        h_el.h4_jet4_phi->Fill(jet.phi(), weight);
        h_el.h4_jet4_btag->Fill(jet.bDiscriminator(bTagName_), weight);
      }
      if ( jets_n >= 5 ) {
        const auto& jet = out_jets->at(4);
        h_el.h4_jet5_m->Fill(jet.mass(), weight);
        h_el.h4_jet5_pt->Fill(jet.pt(), weight);
        h_el.h4_jet5_eta->Fill(jet.eta(), weight);
        h_el.h4_jet5_phi->Fill(jet.phi(), weight);
        h_el.h4_jet5_btag->Fill(jet.bDiscriminator(bTagName_), weight);
      }
      if ( jets_n >= 6 ) {
        const auto& jet = out_jets->at(5);
        h_el.h4_jet6_m->Fill(jet.mass(), weight);
        h_el.h4_jet6_pt->Fill(jet.pt(), weight);
        h_el.h4_jet6_eta->Fill(jet.eta(), weight);
        h_el.h4_jet6_phi->Fill(jet.phi(), weight);
        h_el.h4_jet6_btag->Fill(jet.bDiscriminator(bTagName_), weight);
      }
      h_el.h4_bjets_n->Fill(bjets_n, weight);

      h_el.h4_event_st->Fill(leptons_st+jets_ht+met_pt, weight);
    }
    else if ( channel == 13 ) {
      h_mu.hCutstep->Fill(icutstep, weight);
      h_mu.hCutstepNoweight->Fill(icutstep);
      h_mu.h4_vertex_n->Fill(nVertex, weight);
      h_mu.h4_met_pt->Fill(met_pt, weight);
      h_mu.h4_met_phi->Fill(met_phi, weight);
      h_mu.h4_jets_n->Fill(jets_n, weight);
      h_mu.h4_jets_ht->Fill(jets_ht, weight);
      for ( auto jet : *out_jets ) {
        h_mu.h4_jets_pt->Fill(jet.pt(), weight);
        h_mu.h4_jets_eta->Fill(jet.eta(), weight);
      }
      h_mu.h4_jet1_m->Fill(jet1.mass(), weight);
      h_mu.h4_jet1_pt->Fill(jet1.pt(), weight);
      h_mu.h4_jet1_eta->Fill(jet1.eta(), weight);
      h_mu.h4_jet1_phi->Fill(jet1.phi(), weight);
      h_mu.h4_jet1_btag->Fill(jet1.bDiscriminator(bTagName_), weight);
      if ( jets_n >= 2 ) {
        const auto& jet = out_jets->at(1);
        h_mu.h4_jet2_m->Fill(jet.mass(), weight);
        h_mu.h4_jet2_pt->Fill(jet.pt(), weight);
        h_mu.h4_jet2_eta->Fill(jet.eta(), weight);
        h_mu.h4_jet2_phi->Fill(jet.phi(), weight);
        h_mu.h4_jet2_btag->Fill(jet.bDiscriminator(bTagName_), weight);
      }
      if ( jets_n >= 3 ) {
        const auto& jet = out_jets->at(2);
        h_mu.h4_jet3_m->Fill(jet.mass(), weight);
        h_mu.h4_jet3_pt->Fill(jet.pt(), weight);
        h_mu.h4_jet3_eta->Fill(jet.eta(), weight);
        h_mu.h4_jet3_phi->Fill(jet.phi(), weight);
        h_mu.h4_jet3_btag->Fill(jet.bDiscriminator(bTagName_), weight);
      }
      if ( jets_n >= 4 ) {
        const auto& jet = out_jets->at(3);
        h_mu.h4_jet4_m->Fill(jet.mass(), weight);
        h_mu.h4_jet4_pt->Fill(jet.pt(), weight);
        h_mu.h4_jet4_eta->Fill(jet.eta(), weight);
        h_mu.h4_jet4_phi->Fill(jet.phi(), weight);
        h_mu.h4_jet4_btag->Fill(jet.bDiscriminator(bTagName_), weight);
      }
      if ( jets_n >= 5 ) {
        const auto& jet = out_jets->at(4);
        h_mu.h4_jet5_m->Fill(jet.mass(), weight);
        h_mu.h4_jet5_pt->Fill(jet.pt(), weight);
        h_mu.h4_jet5_eta->Fill(jet.eta(), weight);
        h_mu.h4_jet5_phi->Fill(jet.phi(), weight);
        h_mu.h4_jet5_btag->Fill(jet.bDiscriminator(bTagName_), weight);
      }
      if ( jets_n >= 6 ) {
        const auto& jet = out_jets->at(5);
        h_mu.h4_jet6_m->Fill(jet.mass(), weight);
        h_mu.h4_jet6_pt->Fill(jet.pt(), weight);
        h_mu.h4_jet6_eta->Fill(jet.eta(), weight);
        h_mu.h4_jet6_phi->Fill(jet.phi(), weight);
        h_mu.h4_jet6_btag->Fill(jet.bDiscriminator(bTagName_), weight);
      }
      h_mu.h4_bjets_n->Fill(bjets_n, weight);

      h_mu.h4_event_st->Fill(leptons_st+jets_ht+met_pt, weight);

    }


  } // switch(1)

  // Cutsomized cutflow without z-veto cut to be used in DY estimation and other studies
  for ( int istep=1, nstep=cutstepBits.size(); istep<=nstep; ++istep )
  {
    if ( istep != 2 and !cutstepBits[istep-1] ) break;
  }

  // Fill cut flow 2D plot
  for ( int istep=1, nstep=cutstepBits.size(); istep<=nstep; ++istep )
  {
    const bool res1 = cutstepBits[istep-1];

    // Fill diagonal terms
    if      ( channel == 11 ) {
      h_el.h2Cutstep->Fill(istep, istep, res1*weight);
      h_el.h2CutstepNoweight->Fill(istep, istep, res1);
    }
    else if ( channel == 13 ) {
      h_mu.h2Cutstep->Fill(istep, istep, res1*weight);
      h_mu.h2CutstepNoweight->Fill(istep, istep, res1);
    }

    // Fill correlations and anti-correlations
    for ( int jstep=1; jstep<istep; ++jstep ) {
      const bool res2 = cutstepBits[jstep-1];
      const int result = res1 && res2;
      const int aresult = res1 && !res2;
      if      ( channel == 11 ) {
        h_el.h2Cutstep->Fill(istep, jstep, result*weight);
        h_el.h2CutstepNoweight->Fill(istep, jstep, result);
        h_el.h2Cutstep->Fill(jstep, istep, aresult*weight);
        h_el.h2CutstepNoweight->Fill(jstep, istep, aresult);
      }
      else if ( channel == 13 ) {
        h_mu.h2Cutstep->Fill(istep, jstep, result*weight);
        h_mu.h2CutstepNoweight->Fill(istep, jstep, result);
        h_mu.h2Cutstep->Fill(jstep, istep, aresult*weight);
        h_mu.h2CutstepNoweight->Fill(jstep, istep, aresult);
      }
    }
  }

  event.put(std::auto_ptr<int>(new int((int)channel)), "channel");
  event.put(std::auto_ptr<float>(new float(weight)), "weight");
  event.put(std::auto_ptr<float>(new float(metP4.pt())), "met");
  event.put(std::auto_ptr<float>(new float(metP4.phi())), "metphi");
  event.put(out_leptons, "leptons");
  event.put(out_jets, "jets");

  // Apply filter at the given step.
  if ( cutstep >= applyFilterAt_ ) return true;

  return false;
}

TTLJEventSelector::~TTLJEventSelector()
{
  if ( h_el.hCutstepNoweight ) {
    cout << "---- cut flows without weight ----\n";
    cout << "Step\tel\tmu\n";
    const int n = h_el.hCutstepNoweight->GetNbinsX();
    for ( int i=1; i<=n; ++i ) {
      const string name(h_el.hCutstepNoweight->GetXaxis()->GetBinLabel(i));
      if ( name.empty() ) break;
      cout << name.substr(0, name.find(' '));
      cout << '\t' << h_el.hCutstepNoweight->GetBinContent(i);
      cout << '\t' << h_mu.hCutstepNoweight->GetBinContent(i) << '\n';
    }
    cout << "-----------------------------------\n";
  }
}

#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(TTLJEventSelector);

