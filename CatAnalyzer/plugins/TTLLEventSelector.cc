#include "FWCore/Framework/interface/one/EDFilter.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "CATTools/DataFormats/interface/Muon.h"
#include "CATTools/DataFormats/interface/Electron.h"
#include "CATTools/DataFormats/interface/Jet.h"
#include "CATTools/DataFormats/interface/MET.h"
#include "CATTools/DataFormats/interface/SecVertex.h"

#include "DataFormats/Candidate/interface/LeafCandidate.h"
//#include "DataFormats/Candidate/interface/CompositeCandidate.h"
//#include "DataFormats/Candidate/interface/CompositeRefCandidate.h"
#include "DataFormats/Candidate/interface/CompositePtrCandidate.h"

#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "CommonTools/Utils/interface/TFileDirectory.h"
#include "TH1F.h"
#include "TH2F.h"

using namespace std;

namespace cat {

struct ControlPlots
{
  typedef TH1F* H1;
  typedef TH2F* H2;

  H1 hCutstep, hCutstepNoweight;

  H1 h0a_vertex_n;

  H1 h0b_vertex_n;
  H1 h0b_met_pt, h0b_met_phi;
  H1 h0b_leptons_n;
  H1 h0b_jets_n, h0b_jets_ht;
  H1 h0b_bjets_n, h0b_bjets_ht;

  H1 h0c_vertex_n;
  H1 h0c_met_pt, h0c_met_phi;
  H1 h0c_leptons_n;
  H1 h0c_jets_n, h0c_jets_ht;
  H1 h0c_bjets_n, h0c_bjets_ht;

  H1 h1_vertex_n;
  H1 h1_met_pt, h1_met_phi;
  H1 h1_leptons_n;
  H1 h1_lepton1_pt, h1_lepton1_eta, h1_lepton1_phi, h1_lepton1_q;
  H1 h1_lepton2_pt, h1_lepton2_eta, h1_lepton2_phi, h1_lepton2_q;
  H1 h1_z_m, h1_z_pt, h1_z_eta, h1_z_phi;
  H1 h1_jets_n, h1_jets_ht;
  H1 h1_jet1_m, h1_jet1_pt, h1_jet1_eta, h1_jet1_phi, h1_jet1_btag;
  H1 h1_jet2_m, h1_jet2_pt, h1_jet2_eta, h1_jet2_phi, h1_jet2_btag;
  H1 h1_jet3_m, h1_jet3_pt, h1_jet3_eta, h1_jet3_phi, h1_jet3_btag;
  H1 h1_jet4_m, h1_jet4_pt, h1_jet4_eta, h1_jet4_phi, h1_jet4_btag;
  H1 h1_bjets_n, h1_bjets_ht;
  H1 h1_event_st;

  H1 h2_vertex_n;
  H1 h2_met_pt, h2_met_phi;
  H1 h2_leptons_n;
  H1 h2_lepton1_pt, h2_lepton1_eta, h2_lepton1_phi, h2_lepton1_q;
  H1 h2_lepton2_pt, h2_lepton2_eta, h2_lepton2_phi, h2_lepton2_q;
  H1 h2_z_m, h2_z_pt, h2_z_eta, h2_z_phi;
  H1 h2_jets_n, h2_jets_ht;
  H1 h2_jet1_m, h2_jet1_pt, h2_jet1_eta, h2_jet1_phi, h2_jet1_btag;
  H1 h2_jet2_m, h2_jet2_pt, h2_jet2_eta, h2_jet2_phi, h2_jet2_btag;
  H1 h2_jet3_m, h2_jet3_pt, h2_jet3_eta, h2_jet3_phi, h2_jet3_btag;
  H1 h2_jet4_m, h2_jet4_pt, h2_jet4_eta, h2_jet4_phi, h2_jet4_btag;
  H1 h2_bjets_n, h2_bjets_ht;
  H1 h2_event_st;

  H1 h3_vertex_n;
  H1 h3_met_pt, h3_met_phi;
  H1 h3_z_m, h3_z_pt, h3_z_eta, h3_z_phi;
  H1 h3_jets_n, h3_jets_ht;
  H1 h3_jet1_m, h3_jet1_pt, h3_jet1_eta, h3_jet1_phi, h3_jet1_btag;
  H1 h3_jet2_m, h3_jet2_pt, h3_jet2_eta, h3_jet2_phi, h3_jet2_btag;
  H1 h3_jet3_m, h3_jet3_pt, h3_jet3_eta, h3_jet3_phi, h3_jet3_btag;
  H1 h3_jet4_m, h3_jet4_pt, h3_jet4_eta, h3_jet4_phi, h3_jet4_btag;
  H1 h3_bjets_n, h3_bjets_ht;
  H1 h3_event_st;

  H1 h4_vertex_n;
  H1 h4_met_pt, h4_met_phi;
  H1 h4_z_m, h4_z_pt, h4_z_eta, h4_z_phi;
  H1 h4_jets_n, h4_jets_ht;
  H1 h4_jet1_m, h4_jet1_pt, h4_jet1_eta, h4_jet1_phi, h4_jet1_btag;
  H1 h4_jet2_m, h4_jet2_pt, h4_jet2_eta, h4_jet2_phi, h4_jet2_btag;
  H1 h4_jet3_m, h4_jet3_pt, h4_jet3_eta, h4_jet3_phi, h4_jet3_btag;
  H1 h4_jet4_m, h4_jet4_pt, h4_jet4_eta, h4_jet4_phi, h4_jet4_btag;
  H1 h4_bjets_n, h4_bjets_ht;
  H1 h4_event_st;

  H1 h5_vertex_n;
  H1 h5_met_pt, h5_met_phi;
  H1 h5_z_m, h5_z_pt, h5_z_eta, h5_z_phi;
  H1 h5_jets_n, h5_jets_ht;
  H1 h5_jet1_m, h5_jet1_pt, h5_jet1_eta, h5_jet1_phi, h5_jet1_btag;
  H1 h5_jet2_m, h5_jet2_pt, h5_jet2_eta, h5_jet2_phi, h5_jet2_btag;
  H1 h5_jet3_m, h5_jet3_pt, h5_jet3_eta, h5_jet3_phi, h5_jet3_btag;
  H1 h5_jet4_m, h5_jet4_pt, h5_jet4_eta, h5_jet4_phi, h5_jet4_btag;
  H1 h5_bjets_n, h5_bjets_ht;
  H1 h5_event_st;

  void book(TFileDirectory&& dir)
  {
    const double maxeta = 3;
    const double pi = 3.141592;

    // There are step0a, step0b and step0c cut steps
    // By putting step0a to underflow bin and step0b to -1, step0c to 0,
    // We can start cut steps from 1.
    hCutstep = dir.make<TH1F>("cutstep", "cutstep", 10, -2, 10);
    hCutstepNoweight = dir.make<TH1F>("cutstepNoweight", "cutstepNoweight", 10, -2, 10);

    hCutstep->GetXaxis()->SetBinLabel(1, "S0a all event");
    hCutstep->GetXaxis()->SetBinLabel(2, "S0b Trigger");
    hCutstep->GetXaxis()->SetBinLabel(3, "S0c Event filter");
    hCutstep->GetXaxis()->SetBinLabel(4, "S1 Dilepton");
    hCutstep->GetXaxis()->SetBinLabel(5, "S2 Z veto");
    hCutstep->GetXaxis()->SetBinLabel(6, "S3 nJet2");
    hCutstep->GetXaxis()->SetBinLabel(7, "S4 MET40");
    hCutstep->GetXaxis()->SetBinLabel(8, "S5 nBJet1");
    hCutstep->GetXaxis()->SetBinLabel(9, "S6 nBJet2");

    hCutstepNoweight->GetXaxis()->SetBinLabel(1, "S0a all event");
    hCutstepNoweight->GetXaxis()->SetBinLabel(2, "S0b Trigger");
    hCutstepNoweight->GetXaxis()->SetBinLabel(3, "S0c Event filter");
    hCutstepNoweight->GetXaxis()->SetBinLabel(4, "S1 Dilepton");
    hCutstepNoweight->GetXaxis()->SetBinLabel(5, "S2 Z veto");
    hCutstepNoweight->GetXaxis()->SetBinLabel(6, "S3 nJet2");
    hCutstepNoweight->GetXaxis()->SetBinLabel(7, "S4 MET40");
    hCutstepNoweight->GetXaxis()->SetBinLabel(8, "S5 nBJet1");
    hCutstepNoweight->GetXaxis()->SetBinLabel(9, "S6 nBJet2");

    auto subdir = dir.mkdir("step0a");
    h0a_vertex_n = subdir.make<TH1F>("vertex_n", "vertex_n", 100, 0, 100);

    subdir = dir.mkdir("step0b");
    h0b_vertex_n = subdir.make<TH1F>("vertex_n", "vertex_n", 100, 0, 100);
    h0b_met_pt = subdir.make<TH1F>("met_pt", "met_pt", 1000, 0, 1000);
    h0b_met_phi = subdir.make<TH1F>("met_phi", "met_phi", 100, -pi, pi);
    h0b_leptons_n = subdir.make<TH1F>("leptons_n", "leptons_n", 10, 0, 10);
    h0b_jets_n = subdir.make<TH1F>("jets_n", "jets_n", 10, 0, 10);
    h0b_jets_ht = subdir.make<TH1F>("jets_ht", "jets_ht", 1000, 0, 1000);
    h0b_bjets_n = subdir.make<TH1F>("bjets_n", "bjets_n", 10, 0, 10);
    h0b_bjets_ht = subdir.make<TH1F>("bjets_ht", "bjets_ht", 1000, 0, 1000);

    subdir = dir.mkdir("step0c");
    h0c_vertex_n = subdir.make<TH1F>("vertex_n", "vertex_n", 100, 0, 100);
    h0c_met_pt = subdir.make<TH1F>("met_pt", "met_pt", 1000, 0, 1000);
    h0c_met_phi = subdir.make<TH1F>("met_phi", "met_phi", 100, -pi, pi);
    h0c_leptons_n = subdir.make<TH1F>("leptons_n", "leptons_n", 10, 0, 10);
    h0c_jets_n = subdir.make<TH1F>("jets_n", "jets_n", 10, 0, 10);
    h0c_jets_ht = subdir.make<TH1F>("jets_ht", "jets_ht", 1000, 0, 1000);
    h0c_bjets_n = subdir.make<TH1F>("bjets_n", "bjets_n", 10, 0, 10);
    h0c_bjets_ht = subdir.make<TH1F>("bjets_ht", "bjets_ht", 1000, 0, 1000);

    subdir = dir.mkdir("step1");
    h1_vertex_n = subdir.make<TH1F>("vertex_n", "vertex_n", 100, 0, 100);
    h1_met_pt = subdir.make<TH1F>("met_pt", "met_pt", 1000, 0, 1000);
    h1_met_phi = subdir.make<TH1F>("met_phi", "met_phi", 100, -pi, pi);
    h1_leptons_n = subdir.make<TH1F>("leptons_n", "leptons_n", 10, 0, 10);

    h1_lepton1_pt  = subdir.make<TH1F>("lepton1_pt", "lepton1_pt", 1000, 0, 1000);
    h1_lepton1_eta = subdir.make<TH1F>("lepton1_eta", "lepton1_eta", 100, -maxeta, maxeta);
    h1_lepton1_phi = subdir.make<TH1F>("lepton1_phi", "lepton1_phi", 100, -pi, pi);
    h1_lepton1_q   = subdir.make<TH1F>("lepton1_q", "lepton1_q", 3, -1.5, 1.5);
       
    h1_lepton2_pt  = subdir.make<TH1F>("lepton2_pt", "lepton2_pt", 1000, 0, 1000);
    h1_lepton2_eta = subdir.make<TH1F>("lepton2_eta", "lepton2_eta", 100, -maxeta, maxeta);
    h1_lepton2_phi = subdir.make<TH1F>("lepton2_phi", "lepton2_phi", 100, -pi, pi);
    h1_lepton2_q   = subdir.make<TH1F>("lepton2_q", "lepton2_q", 3, -1.5, 1.5);
       
    h1_z_m   = subdir.make<TH1F>("z_m", "z_m", 1000, 0, 1000);
    h1_z_pt  = subdir.make<TH1F>("z_pt", "z_pt", 1000, 0, 1000);
    h1_z_eta = subdir.make<TH1F>("z_eta", "z_eta", 100, -maxeta, maxeta);
    h1_z_phi = subdir.make<TH1F>("z_phi", "z_phi", 100, -pi, pi);

    h1_jets_n = subdir.make<TH1F>("jets_n", "jets_n", 10, 0, 10);
    h1_jets_ht = subdir.make<TH1F>("jets_ht", "jets_ht", 1000, 0, 1000);

    h1_jet1_m   = subdir.make<TH1F>("jet1_m", "jet1_m", 500, 0, 500);
    h1_jet1_pt  = subdir.make<TH1F>("jet1_pt", "jet1_pt", 1000, 0, 1000);
    h1_jet1_eta = subdir.make<TH1F>("jet1_eta", "jet1_eta", 100, -maxeta, maxeta);
    h1_jet1_phi = subdir.make<TH1F>("jet1_phi", "jet1_phi", 100, -pi, pi);
    h1_jet1_btag = subdir.make<TH1F>("jet1_btag", "jet1_btag", 100, 0, 1);

    h1_jet2_m   = subdir.make<TH1F>("jet2_m", "jet2_m", 500, 0, 500);
    h1_jet2_pt  = subdir.make<TH1F>("jet2_pt", "jet2_pt", 1000, 0, 1000);
    h1_jet2_eta = subdir.make<TH1F>("jet2_eta", "jet2_eta", 100, -maxeta, maxeta);
    h1_jet2_phi = subdir.make<TH1F>("jet2_phi", "jet2_phi", 100, -pi, pi);
    h1_jet2_btag = subdir.make<TH1F>("jet2_btag", "jet2_btag", 100, 0, 1);

    h1_jet3_m   = subdir.make<TH1F>("jet3_m", "jet3_m", 500, 0, 500);
    h1_jet3_pt  = subdir.make<TH1F>("jet3_pt", "jet3_pt", 1000, 0, 1000);
    h1_jet3_eta = subdir.make<TH1F>("jet3_eta", "jet3_eta", 100, -maxeta, maxeta);
    h1_jet3_phi = subdir.make<TH1F>("jet3_phi", "jet3_phi", 100, -pi, pi);
    h1_jet3_btag = subdir.make<TH1F>("jet3_btag", "jet3_btag", 100, 0, 1);

    h1_jet4_m   = subdir.make<TH1F>("jet4_m", "jet4_m", 500, 0, 500);
    h1_jet4_pt  = subdir.make<TH1F>("jet4_pt", "jet4_pt", 1000, 0, 1000);
    h1_jet4_eta = subdir.make<TH1F>("jet4_eta", "jet4_eta", 100, -maxeta, maxeta);
    h1_jet4_phi = subdir.make<TH1F>("jet4_phi", "jet4_phi", 100, -pi, pi);
    h1_jet4_btag = subdir.make<TH1F>("jet4_btag", "jet4_btag", 100, 0, 1);

    h1_bjets_n = subdir.make<TH1F>("bjets_n", "bjets_n", 10, 0, 10);
    h1_bjets_ht = subdir.make<TH1F>("bjets_ht", "bjets_ht", 1000, 0, 1000);

    h1_event_st = subdir.make<TH1F>("event_st", "event_st", 1000, 0, 1000);

    subdir = dir.mkdir("step2");
    h2_vertex_n = subdir.make<TH1F>("vertex_n", "vertex_n", 100, 0, 100);
    h2_met_pt = subdir.make<TH1F>("met_pt", "met_pt", 1000, 0, 1000);
    h2_met_phi = subdir.make<TH1F>("met_phi", "met_phi", 100, -pi, pi);
    h2_leptons_n = subdir.make<TH1F>("leptons_n", "leptons_n", 10, 0, 10);

    h2_lepton1_pt  = subdir.make<TH1F>("lepton1_pt", "lepton1_pt", 1000, 0, 1000);
    h2_lepton1_eta = subdir.make<TH1F>("lepton1_eta", "lepton1_eta", 100, -maxeta, maxeta);
    h2_lepton1_phi = subdir.make<TH1F>("lepton1_phi", "lepton1_phi", 100, -pi, pi);
    h2_lepton1_q   = subdir.make<TH1F>("lepton1_q", "lepton1_q", 3, -1.5, 1.5);

    h2_lepton2_pt  = subdir.make<TH1F>("lepton2_pt", "lepton2_pt", 1000, 0, 1000);
    h2_lepton2_eta = subdir.make<TH1F>("lepton2_eta", "lepton2_eta", 100, -maxeta, maxeta);
    h2_lepton2_phi = subdir.make<TH1F>("lepton2_phi", "lepton2_phi", 100, -pi, pi);
    h2_lepton2_q   = subdir.make<TH1F>("lepton2_q", "lepton2_q", 3, -1.5, 1.5);

    h2_z_m   = subdir.make<TH1F>("z_m", "z_m", 1000, 0, 1000);
    h2_z_pt  = subdir.make<TH1F>("z_pt", "z_pt", 1000, 0, 1000);
    h2_z_eta = subdir.make<TH1F>("z_eta", "z_eta", 100, -maxeta, maxeta);
    h2_z_phi = subdir.make<TH1F>("z_phi", "z_phi", 100, -pi, pi);

    h2_jets_n = subdir.make<TH1F>("jets_n", "jets_n", 10, 0, 10);
    h2_jets_ht = subdir.make<TH1F>("jets_ht", "jets_ht", 1000, 0, 1000);

    h2_jet1_m   = subdir.make<TH1F>("jet1_m", "jet1_m", 500, 0, 500);
    h2_jet1_pt  = subdir.make<TH1F>("jet1_pt", "jet1_pt", 1000, 0, 1000);
    h2_jet1_eta = subdir.make<TH1F>("jet1_eta", "jet1_eta", 100, -maxeta, maxeta);
    h2_jet1_phi = subdir.make<TH1F>("jet1_phi", "jet1_phi", 100, -pi, pi);
    h2_jet1_btag = subdir.make<TH1F>("jet1_btag", "jet1_btag", 100, 0, 1);

    h2_jet2_m   = subdir.make<TH1F>("jet2_m", "jet2_m", 500, 0, 500);
    h2_jet2_pt  = subdir.make<TH1F>("jet2_pt", "jet2_pt", 1000, 0, 1000);
    h2_jet2_eta = subdir.make<TH1F>("jet2_eta", "jet2_eta", 100, -maxeta, maxeta);
    h2_jet2_phi = subdir.make<TH1F>("jet2_phi", "jet2_phi", 100, -pi, pi);
    h2_jet2_btag = subdir.make<TH1F>("jet2_btag", "jet2_btag", 100, 0, 1);

    h2_jet3_m   = subdir.make<TH1F>("jet3_m", "jet3_m", 500, 0, 500);
    h2_jet3_pt  = subdir.make<TH1F>("jet3_pt", "jet3_pt", 1000, 0, 1000);
    h2_jet3_eta = subdir.make<TH1F>("jet3_eta", "jet3_eta", 100, -maxeta, maxeta);
    h2_jet3_phi = subdir.make<TH1F>("jet3_phi", "jet3_phi", 100, -pi, pi);
    h2_jet3_btag = subdir.make<TH1F>("jet3_btag", "jet3_btag", 100, 0, 1);

    h2_jet4_m   = subdir.make<TH1F>("jet4_m", "jet4_m", 500, 0, 500);
    h2_jet4_pt  = subdir.make<TH1F>("jet4_pt", "jet4_pt", 1000, 0, 1000);
    h2_jet4_eta = subdir.make<TH1F>("jet4_eta", "jet4_eta", 100, -maxeta, maxeta);
    h2_jet4_phi = subdir.make<TH1F>("jet4_phi", "jet4_phi", 100, -pi, pi);
    h2_jet4_btag = subdir.make<TH1F>("jet4_btag", "jet4_btag", 100, 0, 1);

    h2_bjets_n = subdir.make<TH1F>("bjets_n", "bjets_n", 10, 0, 10);
    h2_bjets_ht = subdir.make<TH1F>("bjets_ht", "bjets_ht", 1000, 0, 1000);

    h2_event_st = subdir.make<TH1F>("event_st", "event_st", 1000, 0, 1000);

    subdir = dir.mkdir("step3");
    h3_vertex_n = subdir.make<TH1F>("vertex_n", "vertex_n", 100, 0, 100);
    h3_met_pt = subdir.make<TH1F>("met_pt", "met_pt", 1000, 0, 1000);
    h3_met_phi = subdir.make<TH1F>("met_phi", "met_phi", 100, -pi, pi);

    h3_z_m   = subdir.make<TH1F>("z_m", "z_m", 1000, 0, 1000);
    h3_z_pt  = subdir.make<TH1F>("z_pt", "z_pt", 1000, 0, 1000);
    h3_z_eta = subdir.make<TH1F>("z_eta", "z_eta", 100, -maxeta, maxeta);
    h3_z_phi = subdir.make<TH1F>("z_phi", "z_phi", 100, -pi, pi);

    h3_jets_n = subdir.make<TH1F>("jets_n", "jets_n", 10, 0, 10);
    h3_jets_ht = subdir.make<TH1F>("jets_ht", "jets_ht", 1000, 0, 1000);

    h3_jet1_m   = subdir.make<TH1F>("jet1_m", "jet1_m", 500, 0, 500);
    h3_jet1_pt  = subdir.make<TH1F>("jet1_pt", "jet1_pt", 1000, 0, 1000);
    h3_jet1_eta = subdir.make<TH1F>("jet1_eta", "jet1_eta", 100, -maxeta, maxeta);
    h3_jet1_phi = subdir.make<TH1F>("jet1_phi", "jet1_phi", 100, -pi, pi);
    h3_jet1_btag = subdir.make<TH1F>("jet1_btag", "jet1_btag", 100, 0, 1);

    h3_jet2_m   = subdir.make<TH1F>("jet2_m", "jet2_m", 500, 0, 500);
    h3_jet2_pt  = subdir.make<TH1F>("jet2_pt", "jet2_pt", 1000, 0, 1000);
    h3_jet2_eta = subdir.make<TH1F>("jet2_eta", "jet2_eta", 100, -maxeta, maxeta);
    h3_jet2_phi = subdir.make<TH1F>("jet2_phi", "jet2_phi", 100, -pi, pi);
    h3_jet2_btag = subdir.make<TH1F>("jet2_btag", "jet2_btag", 100, 0, 1);

    h3_jet3_m   = subdir.make<TH1F>("jet3_m", "jet3_m", 500, 0, 500);
    h3_jet3_pt  = subdir.make<TH1F>("jet3_pt", "jet3_pt", 1000, 0, 1000);
    h3_jet3_eta = subdir.make<TH1F>("jet3_eta", "jet3_eta", 100, -maxeta, maxeta);
    h3_jet3_phi = subdir.make<TH1F>("jet3_phi", "jet3_phi", 100, -pi, pi);
    h3_jet3_btag = subdir.make<TH1F>("jet3_btag", "jet3_btag", 100, 0, 1);

    h3_jet4_m   = subdir.make<TH1F>("jet4_m", "jet4_m", 500, 0, 500);
    h3_jet4_pt  = subdir.make<TH1F>("jet4_pt", "jet4_pt", 1000, 0, 1000);
    h3_jet4_eta = subdir.make<TH1F>("jet4_eta", "jet4_eta", 100, -maxeta, maxeta);
    h3_jet4_phi = subdir.make<TH1F>("jet4_phi", "jet4_phi", 100, -pi, pi);
    h3_jet4_btag = subdir.make<TH1F>("jet4_btag", "jet4_btag", 100, 0, 1);

    h3_bjets_n = subdir.make<TH1F>("bjets_n", "bjets_n", 10, 0, 10);
    h3_bjets_ht = subdir.make<TH1F>("bjets_ht", "bjets_ht", 1000, 0, 1000);

    h3_event_st = subdir.make<TH1F>("event_st", "event_st", 1000, 0, 1000);

    subdir = dir.mkdir("step4");
    h4_vertex_n = subdir.make<TH1F>("vertex_n", "vertex_n", 100, 0, 100);
    h4_met_pt = subdir.make<TH1F>("met_pt", "met_pt", 1000, 0, 1000);
    h4_met_phi = subdir.make<TH1F>("met_phi", "met_phi", 100, -pi, pi);

    h4_z_m   = subdir.make<TH1F>("z_m", "z_m", 1000, 0, 1000);
    h4_z_pt  = subdir.make<TH1F>("z_pt", "z_pt", 1000, 0, 1000);
    h4_z_eta = subdir.make<TH1F>("z_eta", "z_eta", 100, -maxeta, maxeta);
    h4_z_phi = subdir.make<TH1F>("z_phi", "z_phi", 100, -pi, pi);

    h4_jets_n = subdir.make<TH1F>("jets_n", "jets_n", 10, 0, 10);
    h4_jets_ht = subdir.make<TH1F>("jets_ht", "jets_ht", 1000, 0, 1000);

    h4_jet1_m   = subdir.make<TH1F>("jet1_m", "jet1_m", 500, 0, 500);
    h4_jet1_pt  = subdir.make<TH1F>("jet1_pt", "jet1_pt", 1000, 0, 1000);
    h4_jet1_eta = subdir.make<TH1F>("jet1_eta", "jet1_eta", 100, -maxeta, maxeta);
    h4_jet1_phi = subdir.make<TH1F>("jet1_phi", "jet1_phi", 100, -pi, pi);
    h4_jet1_btag = subdir.make<TH1F>("jet1_btag", "jet1_btag", 100, 0, 1);

    h4_jet2_m   = subdir.make<TH1F>("jet2_m", "jet2_m", 500, 0, 500);
    h4_jet2_pt  = subdir.make<TH1F>("jet2_pt", "jet2_pt", 1000, 0, 1000);
    h4_jet2_eta = subdir.make<TH1F>("jet2_eta", "jet2_eta", 100, -maxeta, maxeta);
    h4_jet2_phi = subdir.make<TH1F>("jet2_phi", "jet2_phi", 100, -pi, pi);
    h4_jet2_btag = subdir.make<TH1F>("jet2_btag", "jet2_btag", 100, 0, 1);

    h4_jet3_m   = subdir.make<TH1F>("jet3_m", "jet3_m", 500, 0, 500);
    h4_jet3_pt  = subdir.make<TH1F>("jet3_pt", "jet3_pt", 1000, 0, 1000);
    h4_jet3_eta = subdir.make<TH1F>("jet3_eta", "jet3_eta", 100, -maxeta, maxeta);
    h4_jet3_phi = subdir.make<TH1F>("jet3_phi", "jet3_phi", 100, -pi, pi);
    h4_jet3_btag = subdir.make<TH1F>("jet3_btag", "jet3_btag", 100, 0, 1);

    h4_jet4_m   = subdir.make<TH1F>("jet4_m", "jet4_m", 500, 0, 500);
    h4_jet4_pt  = subdir.make<TH1F>("jet4_pt", "jet4_pt", 1000, 0, 1000);
    h4_jet4_eta = subdir.make<TH1F>("jet4_eta", "jet4_eta", 100, -maxeta, maxeta);
    h4_jet4_phi = subdir.make<TH1F>("jet4_phi", "jet4_phi", 100, -pi, pi);
    h4_jet4_btag = subdir.make<TH1F>("jet4_btag", "jet4_btag", 100, 0, 1);

    h4_bjets_n = subdir.make<TH1F>("bjets_n", "bjets_n", 10, 0, 10);
    h4_bjets_ht = subdir.make<TH1F>("bjets_ht", "bjets_ht", 1000, 0, 1000);

    h4_event_st = subdir.make<TH1F>("event_st", "event_st", 1000, 0, 1000);

    subdir = dir.mkdir("step5");
    h5_vertex_n = subdir.make<TH1F>("vertex_n", "vertex_n", 100, 0, 100);
    h5_met_pt = subdir.make<TH1F>("met_pt", "met_pt", 1000, 0, 1000);
    h5_met_phi = subdir.make<TH1F>("met_phi", "met_phi", 100, -pi, pi);

    h5_z_m   = subdir.make<TH1F>("z_m", "z_m", 1000, 0, 1000);
    h5_z_pt  = subdir.make<TH1F>("z_pt", "z_pt", 1000, 0, 1000);
    h5_z_eta = subdir.make<TH1F>("z_eta", "z_eta", 100, -maxeta, maxeta);
    h5_z_phi = subdir.make<TH1F>("z_phi", "z_phi", 100, -pi, pi);

    h5_jets_n = subdir.make<TH1F>("jets_n", "jets_n", 10, 0, 10);
    h5_jets_ht = subdir.make<TH1F>("jets_ht", "jets_ht", 1000, 0, 1000);

    h5_jet1_m   = subdir.make<TH1F>("jet1_m", "jet1_m", 500, 0, 500);
    h5_jet1_pt  = subdir.make<TH1F>("jet1_pt", "jet1_pt", 1000, 0, 1000);
    h5_jet1_eta = subdir.make<TH1F>("jet1_eta", "jet1_eta", 100, -maxeta, maxeta);
    h5_jet1_phi = subdir.make<TH1F>("jet1_phi", "jet1_phi", 100, -pi, pi);
    h5_jet1_btag = subdir.make<TH1F>("jet1_btag", "jet1_btag", 100, 0, 1);

    h5_jet2_m   = subdir.make<TH1F>("jet2_m", "jet2_m", 500, 0, 500);
    h5_jet2_pt  = subdir.make<TH1F>("jet2_pt", "jet2_pt", 1000, 0, 1000);
    h5_jet2_eta = subdir.make<TH1F>("jet2_eta", "jet2_eta", 100, -maxeta, maxeta);
    h5_jet2_phi = subdir.make<TH1F>("jet2_phi", "jet2_phi", 100, -pi, pi);
    h5_jet2_btag = subdir.make<TH1F>("jet2_btag", "jet2_btag", 100, 0, 1);

    h5_jet3_m   = subdir.make<TH1F>("jet3_m", "jet3_m", 500, 0, 500);
    h5_jet3_pt  = subdir.make<TH1F>("jet3_pt", "jet3_pt", 1000, 0, 1000);
    h5_jet3_eta = subdir.make<TH1F>("jet3_eta", "jet3_eta", 100, -maxeta, maxeta);
    h5_jet3_phi = subdir.make<TH1F>("jet3_phi", "jet3_phi", 100, -pi, pi);
    h5_jet3_btag = subdir.make<TH1F>("jet3_btag", "jet3_btag", 100, 0, 1);

    h5_jet4_m   = subdir.make<TH1F>("jet4_m", "jet4_m", 500, 0, 500);
    h5_jet4_pt  = subdir.make<TH1F>("jet4_pt", "jet4_pt", 1000, 0, 1000);
    h5_jet4_eta = subdir.make<TH1F>("jet4_eta", "jet4_eta", 100, -maxeta, maxeta);
    h5_jet4_phi = subdir.make<TH1F>("jet4_phi", "jet4_phi", 100, -pi, pi);
    h5_jet4_btag = subdir.make<TH1F>("jet4_btag", "jet4_btag", 100, 0, 1);

    h5_bjets_n = subdir.make<TH1F>("bjets_n", "bjets_n", 10, 0, 10);
    h5_bjets_ht = subdir.make<TH1F>("bjets_ht", "bjets_ht", 1000, 0, 1000);

    h5_event_st = subdir.make<TH1F>("event_st", "event_st", 1000, 0, 1000);

  };
};

class TTLLEventSelector : public edm::one::EDFilter<edm::one::SharedResources>
{
public:
  TTLLEventSelector(const edm::ParameterSet& pset);
  bool filter(edm::Event& event, const edm::EventSetup&) override;

private:
  edm::EDGetTokenT<float> genWeightToken_, pileupWeightToken_;

  edm::EDGetTokenT<cat::MuonCollection> muonToken_;
  edm::EDGetTokenT<cat::ElectronCollection> electronToken_;
  edm::EDGetTokenT<cat::JetCollection> jetToken_;
  edm::EDGetTokenT<cat::METCollection> metToken_;

  edm::EDGetTokenT<int> recoFilterToken_;
  edm::EDGetTokenT<int> trigElElToken_, trigMuMuToken_, trigMuElToken_;
  edm::EDGetTokenT<int> nVertexToken_;

private:
  TH1F* h_weight, * h_pileupWeight, * h_genWeight;
  ControlPlots h_ee, h_mm, h_em;

private:
  bool isGoodMuon(const cat::Muon& mu)
  {
    if ( mu.pt() <= 20 or std::abs(mu.eta()) >= 2.4 ) return false;
    if ( mu.relIso(0.4) >= 0.12 ) return false;
    if ( !mu.isTightMuon() ) return false;
    return true;
  }
  bool isGoodElectron(const cat::Electron& el)
  {
    if ( el.pt() <= 20 or std::abs(el.eta()) >= 2.4 ) return false;
    //if ( el.relIso(0.3) >= 0.11 ) return false;
    if ( !el.electronID(elIdName_) ) return false;
    if ( !el.isPF() or !el.passConversionVeto() ) return false;
    const double scEta = std::abs(el.scEta());
    if ( scEta >= 1.4442 and scEta <= 1.566 ) return false;
    return true;
  }
  bool isBjet(const cat::Jet& jet)
  {
    if ( jet.bDiscriminator(bTagName_) >= bTagMin_ ) return true;
    return false;
  }

private:
  typedef reco::Candidate::LorentzVector LorentzVector;
  enum CHANNEL {
    CH_NONE=-1, CH_MUMU, CH_ELEL, CH_MUEL
  };

  const bool isMC_;
  const int filterCutStepBefore_;
  const std::string elIdName_;
  const std::string bTagName_;
  const double bTagMin_;
};

}

using namespace cat;

TTLLEventSelector::TTLLEventSelector(const edm::ParameterSet& pset):
  isMC_(pset.getParameter<bool>("isMC")),
  filterCutStepBefore_(pset.getParameter<int>("filterCutStepBefore")),
  elIdName_(pset.getParameter<string>("electronID")),
  bTagName_(pset.getParameter<string>("bTagName")),
  bTagMin_(pset.getParameter<double>("bTagMin"))
{
  muonToken_ = consumes<cat::MuonCollection>(pset.getParameter<edm::InputTag>("muons"));
  electronToken_ = consumes<cat::ElectronCollection>(pset.getParameter<edm::InputTag>("electrons"));
  jetToken_ = consumes<cat::JetCollection>(pset.getParameter<edm::InputTag>("jets"));
  metToken_ = consumes<cat::METCollection>(pset.getParameter<edm::InputTag>("mets"));
  nVertexToken_ = consumes<int>(pset.getParameter<edm::InputTag>("nVertex"));

  recoFilterToken_ = consumes<int>(pset.getParameter<edm::InputTag>("filterRECO"));
  trigElElToken_ = consumes<int>(pset.getParameter<edm::InputTag>("trigELEL"));
  trigMuMuToken_ = consumes<int>(pset.getParameter<edm::InputTag>("trigMUMU"));
  trigMuElToken_ = consumes<int>(pset.getParameter<edm::InputTag>("trigMUEL"));

  if ( isMC_ )
  {
    genWeightToken_ = consumes<float>(pset.getParameter<edm::InputTag>("genWeight"));
    pileupWeightToken_ = consumes<float>(pset.getParameter<edm::InputTag>("pileupWeight"));
  }

  // Fill histograms, etc
  usesResource("TFileService");
  edm::Service<TFileService> fs;

  auto doverall = fs->mkdir("overall", "overall");
  h_weight = doverall.make<TH1F>("weight", "weight", 200, -10, 10);
  if ( isMC_ )
  {
    h_genWeight = doverall.make<TH1F>("genWeight", "genWeight", 200, -10, 10);
    h_pileupWeight = doverall.make<TH1F>("pileupWeight", "pileupWeight", 200, -10, 10);
  }

  h_ee.book(fs->mkdir("ee"));
  h_mm.book(fs->mkdir("mm"));
  h_em.book(fs->mkdir("em"));

  produces<int>("channel");
  produces<float>("met");
  produces<float>("metphi");
  produces<std::vector<reco::CandidatePtr> >("leptons");
  produces<std::vector<reco::CandidatePtr> >("jets");
}

bool TTLLEventSelector::filter(edm::Event& event, const edm::EventSetup&)
{
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
  const double met_pt = metP4.pt();
  const double met_phi = metP4.phi();

  std::auto_ptr<std::vector<reco::CandidatePtr> > out_leptons(new std::vector<reco::CandidatePtr>());
  std::auto_ptr<std::vector<reco::CandidatePtr> > out_jets(new std::vector<reco::CandidatePtr>());

  // Compute event weight - from generator, pileup, etc
  double weight = 1.0;
  if ( isMC_ )
  {
    edm::Handle<float> fHandle;

    event.getByToken(genWeightToken_, fHandle);
    const float genWeight = *fHandle;

    event.getByToken(pileupWeightToken_, fHandle);
    const float pileupWeight = *fHandle;

    h_genWeight->Fill(genWeight);
    h_pileupWeight->Fill(pileupWeight);
    weight *= genWeight*pileupWeight;
  }
  h_weight->Fill(weight);

  // Get event filters and triggers
  edm::Handle<int> trigHandle;
  event.getByToken(recoFilterToken_, trigHandle);
  const int isRECOFilterOK = *trigHandle;

  event.getByToken(trigElElToken_, trigHandle);
  const int isTrigElEl = *trigHandle;
  event.getByToken(trigMuMuToken_, trigHandle);
  const int isTrigMuMu = *trigHandle;
  event.getByToken(trigMuElToken_, trigHandle);
  const int isTrigMuEl = *trigHandle;

  // Select good leptons
  double leptons_st = 0;
  for ( int i=0, n=muonHandle->size(); i<n; ++i )
  {
    auto& p = muonHandle->at(i);
    if ( p.pt() < 20 or abs(p.eta()) >= 2.4 ) continue;

    reco::CandidatePtr muonPtr = reco::CandidatePtr(muonHandle, i);

    if ( isGoodMuon(p) )
    {
      leptons_st += p.pt();
      out_leptons->push_back(muonPtr);
    }
  }
  for ( int i=0, n=electronHandle->size(); i<n; ++i )
  {
    auto& p = electronHandle->at(i);
    if ( p.pt() < 20 or std::abs(p.eta()) > 2.4 ) continue;

    reco::CandidatePtr electronPtr = reco::CandidatePtr(electronHandle, i);

    if ( isGoodElectron(p) )
    {
      leptons_st += p.pt();
      out_leptons->push_back(electronPtr);
    }
  }
  auto GtByPtPtr = [](reco::CandidatePtr a, reco::CandidatePtr b){return a->pt() > b->pt();};
  const int leptons_n = out_leptons->size();
  reco::CandidatePtr lepton1, lepton2;
  CHANNEL channel = CH_NONE;
  if ( leptons_n >= 2 ) {
    // Partial sort to select leading 2 leptons
    std::nth_element(out_leptons->begin(), out_leptons->begin()+2, out_leptons->end(), GtByPtPtr);

    // Set lepton1 and 2
    lepton1 = out_leptons->at(0); lepton2 = out_leptons->at(1);

    const int pdgId1 = std::abs(lepton1->pdgId());
    const int pdgId2 = std::abs(lepton2->pdgId());
    // Determine channel
    switch ( pdgId1+pdgId2 )
    {
      case 11+11: { channel = CH_ELEL; break; }
      case 13+13: { channel = CH_MUMU; break; }
      case 11+13: {
        channel = CH_MUEL;
        // Put electron front for emu channel
        if ( pdgId1 == 13 and pdgId2 == 11 ) std::swap(lepton1, lepton2);
      }
    }
  }

  // Select good jets
  int bjets_n = 0;
  double jets_ht = 0, bjets_ht = 0;
  for ( int i=0, n=jetHandle->size(); i<n; ++i )
  {
    auto& jet = jetHandle->at(i);

    if ( jet.pt() < 30 or std::abs(jet.eta()) > 2.5 ) continue;
    if ( leptons_n >= 1 and deltaR(jet.p4(), out_leptons->at(0)->p4()) < 0.5 ) continue;
    if ( leptons_n >= 2 and deltaR(jet.p4(), out_leptons->at(1)->p4()) < 0.5 ) continue;

    jets_ht += jet.pt();
    out_jets->push_back(reco::CandidatePtr(jetHandle, i));
    if ( isBjet(jet) )
    {
      bjets_ht += jet.pt();
      ++bjets_n;
    }
  }
  const int jets_n = out_jets->size();
  std::sort(out_jets->begin(), out_jets->end(), GtByPtPtr);

  // Check cut steps and fill histograms
  h_ee.hCutstep->Fill(-2, weight);
  h_ee.hCutstepNoweight->Fill(-2);
  h_ee.h0a_vertex_n->Fill(weight);

  h_mm.hCutstep->Fill(-2, weight);
  h_mm.hCutstepNoweight->Fill(-2);
  h_mm.h0a_vertex_n->Fill(weight);

  h_mm.hCutstep->Fill(-2, weight);
  h_mm.hCutstepNoweight->Fill(-2);
  h_mm.h0a_vertex_n->Fill(weight);

  // ElEl channel Cutstep 0b with trigger requirements
  int cutstep_ee = -2;
  if ( isTrigElEl )
  {
    ++cutstep_ee;
    h_ee.hCutstep->Fill(-1, weight);
    h_ee.hCutstepNoweight->Fill(-1);
    h_ee.h0b_vertex_n->Fill(weight);
    h_ee.h0b_met_pt->Fill(met_pt, weight);
    h_ee.h0b_met_phi->Fill(met_phi, weight);
    h_ee.h0b_leptons_n->Fill(leptons_n, weight);
    h_ee.h0b_jets_n->Fill(jets_n, weight);
    h_ee.h0b_bjets_n->Fill(bjets_n, weight);
    h_ee.h0b_bjets_ht->Fill(bjets_ht, weight);
    h_ee.h0b_jets_ht->Fill(jets_ht, weight);

    // Cutstep 0c with reco filters
    if ( isRECOFilterOK )
    {
      ++cutstep_ee;
      h_ee.hCutstep->Fill(0., weight);
      h_ee.hCutstepNoweight->Fill(0.);
      h_ee.h0c_vertex_n->Fill(weight);
      h_ee.h0c_met_pt->Fill(met_pt, weight);
      h_ee.h0c_met_phi->Fill(met_phi, weight);
      h_ee.h0c_leptons_n->Fill(leptons_n, weight);
      h_ee.h0c_jets_n->Fill(jets_n, weight);
      h_ee.h0c_bjets_n->Fill(bjets_n, weight);
      h_ee.h0c_bjets_ht->Fill(bjets_ht, weight);
      h_ee.h0c_jets_ht->Fill(jets_ht, weight);
    }
  }
  // MuMu channel Cutstep 0b with trigger requirements
  int cutstep_mm = -2;
  if ( isTrigMuMu )
  {
    ++cutstep_mm;
    h_mm.hCutstep->Fill(-1, weight);
    h_mm.hCutstepNoweight->Fill(-1);
    h_mm.h0b_vertex_n->Fill(weight);
    h_mm.h0b_met_pt->Fill(met_pt, weight);
    h_mm.h0b_met_phi->Fill(met_phi, weight);
    h_mm.h0b_leptons_n->Fill(leptons_n, weight);
    h_mm.h0b_jets_n->Fill(jets_n, weight);
    h_mm.h0b_bjets_n->Fill(bjets_n, weight);
    h_mm.h0b_bjets_ht->Fill(bjets_ht, weight);
    h_mm.h0b_jets_ht->Fill(jets_ht, weight);

    // Cutstep 0c with reco filters
    if ( isRECOFilterOK )
    {
      ++cutstep_mm;
      h_mm.hCutstep->Fill(0., weight);
      h_mm.hCutstepNoweight->Fill(0.);
      h_mm.h0c_vertex_n->Fill(weight);
      h_mm.h0c_met_pt->Fill(met_pt, weight);
      h_mm.h0c_met_phi->Fill(met_phi, weight);
      h_mm.h0c_leptons_n->Fill(leptons_n, weight);
      h_mm.h0c_jets_n->Fill(jets_n, weight);
      h_mm.h0c_bjets_n->Fill(bjets_n, weight);
      h_mm.h0c_bjets_ht->Fill(bjets_ht, weight);
      h_mm.h0c_jets_ht->Fill(jets_ht, weight);
    }
  }
  // MuEl channel Cutstep 0b with trigger requirements
  int cutstep_em = -2;
  if ( isTrigMuEl )
  {
    ++cutstep_em;
    h_em.hCutstep->Fill(-1, weight);
    h_em.hCutstepNoweight->Fill(-1);
    h_em.h0b_vertex_n->Fill(weight);
    h_em.h0b_met_pt->Fill(met_pt, weight);
    h_em.h0b_met_phi->Fill(met_phi, weight);
    h_em.h0b_leptons_n->Fill(leptons_n, weight);
    h_em.h0b_jets_n->Fill(jets_n, weight);
    h_em.h0b_bjets_n->Fill(bjets_n, weight);
    h_em.h0b_bjets_ht->Fill(bjets_ht, weight);
    h_em.h0b_jets_ht->Fill(jets_ht, weight);

    // Cutstep 0c with reco filters
    if ( isRECOFilterOK )
    {
      ++cutstep_em;
      h_em.hCutstep->Fill(0., weight);
      h_em.hCutstepNoweight->Fill(0.);
      h_em.h0c_vertex_n->Fill(weight);
      h_em.h0c_met_pt->Fill(met_pt, weight);
      h_em.h0c_met_phi->Fill(met_phi, weight);
      h_em.h0c_leptons_n->Fill(leptons_n, weight);
      h_em.h0c_jets_n->Fill(jets_n, weight);
      h_em.h0c_bjets_n->Fill(bjets_n, weight);
      h_em.h0c_bjets_ht->Fill(bjets_ht, weight);
      h_em.h0c_jets_ht->Fill(jets_ht, weight);
    }
  }

  int cutstep = 1;
  if ( leptons_n < 2 or lepton1.isNull() or lepton2.isNull() or
       (lepton1->p4()+lepton2->p4()).mass() < 20 or
       lepton1->charge()+lepton2->charge() != 0 )
  {
    channel = CH_NONE;
    cutstep = std::max(cutstep_ee, std::max(cutstep_mm, cutstep_em)); // reset the cut step
  }
  else
  {
    const auto zP4 = lepton1->p4()+lepton2->p4();
    const auto z_m = zP4.mass();
    switch ( channel )
    {
      ////////////////////////////////////////
      // First scan cut flow for ee channel //
      ////////////////////////////////////////
      case CH_ELEL: {
        cutstep = cutstep_ee+1;

        h_ee.hCutstep->Fill(cutstep, weight);
        h_ee.hCutstepNoweight->Fill(cutstep);
        h_ee.h1_vertex_n->Fill(weight);
        h_ee.h1_met_pt->Fill(met_pt, weight);
        h_ee.h1_met_phi->Fill(met_phi, weight);
        h_ee.h1_leptons_n->Fill(leptons_n, weight);
        h_ee.h1_lepton1_pt->Fill(lepton1->pt(), weight);
        h_ee.h1_lepton1_eta->Fill(lepton1->eta(), weight);
        h_ee.h1_lepton1_phi->Fill(lepton1->phi(), weight);
        h_ee.h1_lepton1_q->Fill(lepton1->charge(), weight);
        h_ee.h1_lepton2_pt->Fill(lepton2->pt(), weight);
        h_ee.h1_lepton2_eta->Fill(lepton2->eta(), weight);
        h_ee.h1_lepton2_phi->Fill(lepton2->phi(), weight);
        h_ee.h1_lepton2_q->Fill(lepton2->charge(), weight);
        h_ee.h1_z_m->Fill(z_m, weight);
        h_ee.h1_z_pt->Fill(zP4.pt(), weight);
        h_ee.h1_z_eta->Fill(zP4.eta(), weight);
        h_ee.h1_z_phi->Fill(zP4.phi(), weight);
        h_ee.h1_jets_n->Fill(jets_n, weight);
        h_ee.h1_jets_ht->Fill(jets_ht, weight);
        if ( jets_n >= 1 )
        {
          const auto& jet = dynamic_cast<const cat::Jet&>(*out_jets->at(0));
          h_ee.h1_jet1_m->Fill(jet.mass(), weight);
          h_ee.h1_jet1_pt->Fill(jet.pt(), weight);
          h_ee.h1_jet1_eta->Fill(jet.eta(), weight);
          h_ee.h1_jet1_phi->Fill(jet.phi(), weight);
          h_ee.h1_jet1_btag->Fill(jet.bDiscriminator(bTagName_), weight);
        }
        if ( jets_n >= 2 )
        {
          const auto& jet = dynamic_cast<const cat::Jet&>(*out_jets->at(1));
          h_ee.h1_jet2_m->Fill(jet.mass(), weight);
          h_ee.h1_jet2_pt->Fill(jet.pt(), weight);
          h_ee.h1_jet2_eta->Fill(jet.eta(), weight);
          h_ee.h1_jet2_phi->Fill(jet.phi(), weight);
          h_ee.h1_jet2_btag->Fill(jet.bDiscriminator(bTagName_), weight);
        }
        if ( jets_n >= 3 )
        {
          const auto& jet = dynamic_cast<const cat::Jet&>(*out_jets->at(2));
          h_ee.h1_jet3_m->Fill(jet.mass(), weight);
          h_ee.h1_jet3_pt->Fill(jet.pt(), weight);
          h_ee.h1_jet3_eta->Fill(jet.eta(), weight);
          h_ee.h1_jet3_phi->Fill(jet.phi(), weight);
          h_ee.h1_jet3_btag->Fill(jet.bDiscriminator(bTagName_), weight);
        }
        if ( jets_n >= 4 )
        {
          const auto& jet = dynamic_cast<const cat::Jet&>(*out_jets->at(3));
          h_ee.h1_jet4_m->Fill(jet.mass(), weight);
          h_ee.h1_jet4_pt->Fill(jet.pt(), weight);
          h_ee.h1_jet4_eta->Fill(jet.eta(), weight);
          h_ee.h1_jet4_phi->Fill(jet.phi(), weight);
          h_ee.h1_jet4_btag->Fill(jet.bDiscriminator(bTagName_), weight);
        }
        h_ee.h1_bjets_n->Fill(bjets_n, weight);
        h_ee.h1_bjets_ht->Fill(bjets_ht, weight);

        h_ee.h1_event_st->Fill(leptons_st+jets_ht+met_pt, weight);

        // Step2 Z mass veto
        if ( 76 <= z_m and z_m <= 106 ) break;
        ++cutstep;
        h_ee.hCutstep->Fill(cutstep, weight);
        h_ee.hCutstepNoweight->Fill(cutstep);
        h_ee.h2_vertex_n->Fill(weight);
        h_ee.h2_met_pt->Fill(met_pt, weight);
        h_ee.h2_met_phi->Fill(met_phi, weight);
        h_ee.h2_leptons_n->Fill(leptons_n, weight);
        h_ee.h2_lepton1_pt->Fill(lepton1->pt(), weight);
        h_ee.h2_lepton1_eta->Fill(lepton1->eta(), weight);
        h_ee.h2_lepton1_phi->Fill(lepton1->phi(), weight);
        h_ee.h2_lepton1_q->Fill(lepton1->charge(), weight);
        h_ee.h2_lepton2_pt->Fill(lepton2->pt(), weight);
        h_ee.h2_lepton2_eta->Fill(lepton2->eta(), weight);
        h_ee.h2_lepton2_phi->Fill(lepton2->phi(), weight);
        h_ee.h2_lepton2_q->Fill(lepton2->charge(), weight);
        h_ee.h2_z_m->Fill(zP4.mass(), weight);
        h_ee.h2_z_pt->Fill(zP4.pt(), weight);
        h_ee.h2_z_eta->Fill(zP4.eta(), weight);
        h_ee.h2_z_phi->Fill(zP4.phi(), weight);
        h_ee.h2_jets_n->Fill(jets_n, weight);
        h_ee.h2_jets_ht->Fill(jets_ht, weight);
        if ( jets_n >= 1 )
        {
          const auto& jet = dynamic_cast<const cat::Jet&>(*out_jets->at(0));
          h_ee.h2_jet1_m->Fill(jet.mass(), weight);
          h_ee.h2_jet1_pt->Fill(jet.pt(), weight);
          h_ee.h2_jet1_eta->Fill(jet.eta(), weight);
          h_ee.h2_jet1_phi->Fill(jet.phi(), weight);
          h_ee.h2_jet1_btag->Fill(jet.bDiscriminator(bTagName_), weight);
        }
        if ( jets_n >= 2 )
        {
          const auto& jet = dynamic_cast<const cat::Jet&>(*out_jets->at(1));
          h_ee.h2_jet2_m->Fill(jet.mass(), weight);
          h_ee.h2_jet2_pt->Fill(jet.pt(), weight);
          h_ee.h2_jet2_eta->Fill(jet.eta(), weight);
          h_ee.h2_jet2_phi->Fill(jet.phi(), weight);
          h_ee.h2_jet2_btag->Fill(jet.bDiscriminator(bTagName_), weight);
        }
        if ( jets_n >= 3 )
        {
          const auto& jet = dynamic_cast<const cat::Jet&>(*out_jets->at(2));
          h_ee.h2_jet3_m->Fill(jet.mass(), weight);
          h_ee.h2_jet3_pt->Fill(jet.pt(), weight);
          h_ee.h2_jet3_eta->Fill(jet.eta(), weight);
          h_ee.h2_jet3_phi->Fill(jet.phi(), weight);
          h_ee.h2_jet3_btag->Fill(jet.bDiscriminator(bTagName_), weight);
        }
        if ( jets_n >= 4 )
        {
          const auto& jet = dynamic_cast<const cat::Jet&>(*out_jets->at(3));
          h_ee.h2_jet4_m->Fill(jet.mass(), weight);
          h_ee.h2_jet4_pt->Fill(jet.pt(), weight);
          h_ee.h2_jet4_eta->Fill(jet.eta(), weight);
          h_ee.h2_jet4_phi->Fill(jet.phi(), weight);
          h_ee.h2_jet4_btag->Fill(jet.bDiscriminator(bTagName_), weight);
        }
        h_ee.h2_bjets_n->Fill(bjets_n, weight);
        h_ee.h2_bjets_ht->Fill(bjets_ht, weight);

        h_ee.h2_event_st->Fill(leptons_st+jets_ht+met_pt, weight);

        // Step3 Minimal jet multiplicity
        if ( jets_n < 2 ) break;
        ++cutstep;
        h_ee.hCutstep->Fill(cutstep, weight);
        h_ee.hCutstepNoweight->Fill(cutstep);
        h_ee.h3_vertex_n->Fill(weight);
        h_ee.h3_met_pt->Fill(met_pt, weight);
        h_ee.h3_met_phi->Fill(met_phi, weight);
        h_ee.h3_z_m->Fill(zP4.mass(), weight);
        h_ee.h3_z_pt->Fill(zP4.pt(), weight);
        h_ee.h3_z_eta->Fill(zP4.eta(), weight);
        h_ee.h3_z_phi->Fill(zP4.phi(), weight);
        h_ee.h3_jets_n->Fill(jets_n, weight);
        h_ee.h3_jets_ht->Fill(jets_ht, weight);
        const auto& jet1 = dynamic_cast<const cat::Jet&>(*out_jets->at(0));
        h_ee.h3_jet1_m->Fill(jet1.mass(), weight);
        h_ee.h3_jet1_pt->Fill(jet1.pt(), weight);
        h_ee.h3_jet1_eta->Fill(jet1.eta(), weight);
        h_ee.h3_jet1_phi->Fill(jet1.phi(), weight);
        h_ee.h3_jet1_btag->Fill(jet1.bDiscriminator(bTagName_), weight);
        const auto& jet2 = dynamic_cast<const cat::Jet&>(*out_jets->at(1));
        h_ee.h3_jet2_m->Fill(jet2.mass(), weight);
        h_ee.h3_jet2_pt->Fill(jet2.pt(), weight);
        h_ee.h3_jet2_eta->Fill(jet2.eta(), weight);
        h_ee.h3_jet2_phi->Fill(jet2.phi(), weight);
        h_ee.h3_jet2_btag->Fill(jet2.bDiscriminator(bTagName_), weight);
        if ( jets_n >= 3 )
        {
          const auto& jet = dynamic_cast<const cat::Jet&>(*out_jets->at(2));
          h_ee.h3_jet3_m->Fill(jet.mass(), weight);
          h_ee.h3_jet3_pt->Fill(jet.pt(), weight);
          h_ee.h3_jet3_eta->Fill(jet.eta(), weight);
          h_ee.h3_jet3_phi->Fill(jet.phi(), weight);
          h_ee.h3_jet3_btag->Fill(jet.bDiscriminator(bTagName_), weight);
        }
        if ( jets_n >= 4 )
        {
          const auto& jet = dynamic_cast<const cat::Jet&>(*out_jets->at(3));
          h_ee.h3_jet4_m->Fill(jet.mass(), weight);
          h_ee.h3_jet4_pt->Fill(jet.pt(), weight);
          h_ee.h3_jet4_eta->Fill(jet.eta(), weight);
          h_ee.h3_jet4_phi->Fill(jet.phi(), weight);
          h_ee.h3_jet4_btag->Fill(jet.bDiscriminator(bTagName_), weight);
        }
        h_ee.h3_bjets_n->Fill(bjets_n, weight);
        h_ee.h3_bjets_ht->Fill(bjets_ht, weight);

        h_ee.h3_event_st->Fill(leptons_st+jets_ht+met_pt, weight);

        // Step4 Missing transverse momentum
        if ( met_pt < 40 ) break;
        ++cutstep;
        h_ee.hCutstep->Fill(cutstep, weight);
        h_ee.hCutstepNoweight->Fill(cutstep);
        h_ee.h4_vertex_n->Fill(weight);
        h_ee.h4_met_pt->Fill(met_pt, weight);
        h_ee.h4_met_phi->Fill(met_phi, weight);
        h_ee.h4_z_m->Fill(zP4.mass(), weight);
        h_ee.h4_z_pt->Fill(zP4.pt(), weight);
        h_ee.h4_z_eta->Fill(zP4.eta(), weight);
        h_ee.h4_z_phi->Fill(zP4.phi(), weight);
        h_ee.h4_jets_n->Fill(jets_n, weight);
        h_ee.h4_jets_ht->Fill(jets_ht, weight);
        h_ee.h4_jet1_m->Fill(jet1.mass(), weight);
        h_ee.h4_jet1_pt->Fill(jet1.pt(), weight);
        h_ee.h4_jet1_eta->Fill(jet1.eta(), weight);
        h_ee.h4_jet1_phi->Fill(jet1.phi(), weight);
        h_ee.h4_jet1_btag->Fill(jet1.bDiscriminator(bTagName_), weight);
        h_ee.h4_jet2_m->Fill(jet2.mass(), weight);
        h_ee.h4_jet2_pt->Fill(jet2.pt(), weight);
        h_ee.h4_jet2_eta->Fill(jet2.eta(), weight);
        h_ee.h4_jet2_phi->Fill(jet2.phi(), weight);
        h_ee.h4_jet2_btag->Fill(jet2.bDiscriminator(bTagName_), weight);
        if ( jets_n >= 3 )
        {
          const auto& jet = dynamic_cast<const cat::Jet&>(*out_jets->at(2));
          h_ee.h4_jet3_m->Fill(jet.mass(), weight);
          h_ee.h4_jet3_pt->Fill(jet.pt(), weight);
          h_ee.h4_jet3_eta->Fill(jet.eta(), weight);
          h_ee.h4_jet3_phi->Fill(jet.phi(), weight);
          h_ee.h4_jet3_btag->Fill(jet.bDiscriminator(bTagName_), weight);
        }
        if ( jets_n >= 4 )
        {
          const auto& jet = dynamic_cast<const cat::Jet&>(*out_jets->at(3));
          h_ee.h4_jet4_m->Fill(jet.mass(), weight);
          h_ee.h4_jet4_pt->Fill(jet.pt(), weight);
          h_ee.h4_jet4_eta->Fill(jet.eta(), weight);
          h_ee.h4_jet4_phi->Fill(jet.phi(), weight);
          h_ee.h4_jet4_btag->Fill(jet.bDiscriminator(bTagName_), weight);
        }
        h_ee.h4_bjets_n->Fill(bjets_n, weight);
        h_ee.h4_bjets_ht->Fill(bjets_ht, weight);

        h_ee.h4_event_st->Fill(leptons_st+jets_ht+met_pt, weight);

        // Step5 one b jet
        if ( bjets_n < 1 ) break;
        ++cutstep;
        h_ee.hCutstep->Fill(cutstep, weight);
        h_ee.hCutstepNoweight->Fill(cutstep);
        h_ee.h5_vertex_n->Fill(weight);
        h_ee.h5_met_pt->Fill(met_pt, weight);
        h_ee.h5_met_phi->Fill(met_phi, weight);
        h_ee.h5_z_m->Fill(zP4.mass(), weight);
        h_ee.h5_z_pt->Fill(zP4.pt(), weight);
        h_ee.h5_z_eta->Fill(zP4.eta(), weight);
        h_ee.h5_z_phi->Fill(zP4.phi(), weight);
        h_ee.h5_jets_n->Fill(jets_n, weight);
        h_ee.h5_jets_ht->Fill(jets_ht, weight);
        h_ee.h5_jet1_m->Fill(jet1.mass(), weight);
        h_ee.h5_jet1_pt->Fill(jet1.pt(), weight);
        h_ee.h5_jet1_eta->Fill(jet1.eta(), weight);
        h_ee.h5_jet1_phi->Fill(jet1.phi(), weight);
        h_ee.h5_jet1_btag->Fill(jet1.bDiscriminator(bTagName_), weight);
        h_ee.h5_jet2_m->Fill(jet2.mass(), weight);
        h_ee.h5_jet2_pt->Fill(jet2.pt(), weight);
        h_ee.h5_jet2_eta->Fill(jet2.eta(), weight);
        h_ee.h5_jet2_phi->Fill(jet2.phi(), weight);
        h_ee.h5_jet2_btag->Fill(jet2.bDiscriminator(bTagName_), weight);
        if ( jets_n >= 3 )
        {
          const auto& jet = dynamic_cast<const cat::Jet&>(*out_jets->at(2));
          h_ee.h5_jet3_m->Fill(jet.mass(), weight);
          h_ee.h5_jet3_pt->Fill(jet.pt(), weight);
          h_ee.h5_jet3_eta->Fill(jet.eta(), weight);
          h_ee.h5_jet3_phi->Fill(jet.phi(), weight);
          h_ee.h5_jet3_btag->Fill(jet.bDiscriminator(bTagName_), weight);
        }
        if ( jets_n >= 4 )
        {
          const auto& jet = dynamic_cast<const cat::Jet&>(*out_jets->at(3));
          h_ee.h5_jet4_m->Fill(jet.mass(), weight);
          h_ee.h5_jet4_pt->Fill(jet.pt(), weight);
          h_ee.h5_jet4_eta->Fill(jet.eta(), weight);
          h_ee.h5_jet4_phi->Fill(jet.phi(), weight);
          h_ee.h5_jet4_btag->Fill(jet.bDiscriminator(bTagName_), weight);
        }
        h_ee.h5_bjets_n->Fill(bjets_n, weight);
        h_ee.h5_bjets_ht->Fill(bjets_ht, weight);

        h_ee.h5_event_st->Fill(leptons_st+jets_ht+met_pt, weight);
        break;
      }
      /////////////////////////////////////////
      // Next scan cut flow for mumu channel //
      /////////////////////////////////////////
      case CH_MUMU: {
        cutstep = cutstep_mm+1;

        h_mm.hCutstep->Fill(cutstep, weight);
        h_mm.hCutstepNoweight->Fill(cutstep);
        h_mm.h1_vertex_n->Fill(weight);
        h_mm.h1_met_pt->Fill(met_pt, weight);
        h_mm.h1_met_phi->Fill(met_phi, weight);
        h_mm.h1_leptons_n->Fill(leptons_n, weight);
        h_mm.h1_lepton1_pt->Fill(lepton1->pt(), weight);
        h_mm.h1_lepton1_eta->Fill(lepton1->eta(), weight);
        h_mm.h1_lepton1_phi->Fill(lepton1->phi(), weight);
        h_mm.h1_lepton1_q->Fill(lepton1->charge(), weight);
        h_mm.h1_lepton2_pt->Fill(lepton2->pt(), weight);
        h_mm.h1_lepton2_eta->Fill(lepton2->eta(), weight);
        h_mm.h1_lepton2_phi->Fill(lepton2->phi(), weight);
        h_mm.h1_lepton2_q->Fill(lepton2->charge(), weight);
        h_mm.h1_z_m->Fill(z_m, weight);
        h_mm.h1_z_pt->Fill(zP4.pt(), weight);
        h_mm.h1_z_eta->Fill(zP4.eta(), weight);
        h_mm.h1_z_phi->Fill(zP4.phi(), weight);
        h_mm.h1_jets_n->Fill(jets_n, weight);
        h_mm.h1_jets_ht->Fill(jets_ht, weight);
        if ( jets_n >= 1 )
        {
          const auto& jet = dynamic_cast<const cat::Jet&>(*out_jets->at(0));
          h_mm.h1_jet1_m->Fill(jet.mass(), weight);
          h_mm.h1_jet1_pt->Fill(jet.pt(), weight);
          h_mm.h1_jet1_eta->Fill(jet.eta(), weight);
          h_mm.h1_jet1_phi->Fill(jet.phi(), weight);
          h_mm.h1_jet1_btag->Fill(jet.bDiscriminator(bTagName_), weight);
        }
        if ( jets_n >= 2 )
        {
          const auto& jet = dynamic_cast<const cat::Jet&>(*out_jets->at(1));
          h_mm.h1_jet2_m->Fill(jet.mass(), weight);
          h_mm.h1_jet2_pt->Fill(jet.pt(), weight);
          h_mm.h1_jet2_eta->Fill(jet.eta(), weight);
          h_mm.h1_jet2_phi->Fill(jet.phi(), weight);
          h_mm.h1_jet2_btag->Fill(jet.bDiscriminator(bTagName_), weight);
        }
        if ( jets_n >= 3 )
        {
          const auto& jet = dynamic_cast<const cat::Jet&>(*out_jets->at(2));
          h_mm.h1_jet3_m->Fill(jet.mass(), weight);
          h_mm.h1_jet3_pt->Fill(jet.pt(), weight);
          h_mm.h1_jet3_eta->Fill(jet.eta(), weight);
          h_mm.h1_jet3_phi->Fill(jet.phi(), weight);
          h_mm.h1_jet3_btag->Fill(jet.bDiscriminator(bTagName_), weight);
        }
        if ( jets_n >= 4 )
        {
          const auto& jet = dynamic_cast<const cat::Jet&>(*out_jets->at(3));
          h_mm.h1_jet4_m->Fill(jet.mass(), weight);
          h_mm.h1_jet4_pt->Fill(jet.pt(), weight);
          h_mm.h1_jet4_eta->Fill(jet.eta(), weight);
          h_mm.h1_jet4_phi->Fill(jet.phi(), weight);
          h_mm.h1_jet4_btag->Fill(jet.bDiscriminator(bTagName_), weight);
        }
        h_mm.h1_bjets_n->Fill(bjets_n, weight);
        h_mm.h1_bjets_ht->Fill(bjets_ht, weight);

        h_mm.h1_event_st->Fill(leptons_st+jets_ht+met_pt, weight);

        // Step2 Z mass veto
        if ( 76 <= z_m and z_m <= 106 ) break;
        ++cutstep;
        h_mm.hCutstep->Fill(cutstep, weight);
        h_mm.hCutstepNoweight->Fill(cutstep);
        h_mm.h2_vertex_n->Fill(weight);
        h_mm.h2_met_pt->Fill(met_pt, weight);
        h_mm.h2_met_phi->Fill(met_phi, weight);
        h_mm.h2_leptons_n->Fill(leptons_n, weight);
        h_mm.h2_lepton1_pt->Fill(lepton1->pt(), weight);
        h_mm.h2_lepton1_eta->Fill(lepton1->eta(), weight);
        h_mm.h2_lepton1_phi->Fill(lepton1->phi(), weight);
        h_mm.h2_lepton1_q->Fill(lepton1->charge(), weight);
        h_mm.h2_lepton2_pt->Fill(lepton2->pt(), weight);
        h_mm.h2_lepton2_eta->Fill(lepton2->eta(), weight);
        h_mm.h2_lepton2_phi->Fill(lepton2->phi(), weight);
        h_mm.h2_lepton2_q->Fill(lepton2->charge(), weight);
        h_mm.h2_z_m->Fill(zP4.mass(), weight);
        h_mm.h2_z_pt->Fill(zP4.pt(), weight);
        h_mm.h2_z_eta->Fill(zP4.eta(), weight);
        h_mm.h2_z_phi->Fill(zP4.phi(), weight);
        h_mm.h2_jets_n->Fill(jets_n, weight);
        h_mm.h2_jets_ht->Fill(jets_ht, weight);
        if ( jets_n >= 1 )
        {
          const auto& jet = dynamic_cast<const cat::Jet&>(*out_jets->at(0));
          h_mm.h2_jet1_m->Fill(jet.mass(), weight);
          h_mm.h2_jet1_pt->Fill(jet.pt(), weight);
          h_mm.h2_jet1_eta->Fill(jet.eta(), weight);
          h_mm.h2_jet1_phi->Fill(jet.phi(), weight);
          h_mm.h2_jet1_btag->Fill(jet.bDiscriminator(bTagName_), weight);
        }
        if ( jets_n >= 2 )
        {
          const auto& jet = dynamic_cast<const cat::Jet&>(*out_jets->at(1));
          h_mm.h2_jet2_m->Fill(jet.mass(), weight);
          h_mm.h2_jet2_pt->Fill(jet.pt(), weight);
          h_mm.h2_jet2_eta->Fill(jet.eta(), weight);
          h_mm.h2_jet2_phi->Fill(jet.phi(), weight);
          h_mm.h2_jet2_btag->Fill(jet.bDiscriminator(bTagName_), weight);
        }
        if ( jets_n >= 3 )
        {
          const auto& jet = dynamic_cast<const cat::Jet&>(*out_jets->at(2));
          h_mm.h2_jet3_m->Fill(jet.mass(), weight);
          h_mm.h2_jet3_pt->Fill(jet.pt(), weight);
          h_mm.h2_jet3_eta->Fill(jet.eta(), weight);
          h_mm.h2_jet3_phi->Fill(jet.phi(), weight);
          h_mm.h2_jet3_btag->Fill(jet.bDiscriminator(bTagName_), weight);
        }
        if ( jets_n >= 4 )
        {
          const auto& jet = dynamic_cast<const cat::Jet&>(*out_jets->at(3));
          h_mm.h2_jet4_m->Fill(jet.mass(), weight);
          h_mm.h2_jet4_pt->Fill(jet.pt(), weight);
          h_mm.h2_jet4_eta->Fill(jet.eta(), weight);
          h_mm.h2_jet4_phi->Fill(jet.phi(), weight);
          h_mm.h2_jet4_btag->Fill(jet.bDiscriminator(bTagName_), weight);
        }
        h_mm.h2_bjets_n->Fill(bjets_n, weight);
        h_mm.h2_bjets_ht->Fill(bjets_ht, weight);

        h_mm.h2_event_st->Fill(leptons_st+jets_ht+met_pt, weight);

        // Step3 Minimal jet multiplicity
        if ( jets_n < 2 ) break;
        ++cutstep;
        h_mm.hCutstep->Fill(cutstep, weight);
        h_mm.hCutstepNoweight->Fill(cutstep);
        h_mm.h3_vertex_n->Fill(weight);
        h_mm.h3_met_pt->Fill(met_pt, weight);
        h_mm.h3_met_phi->Fill(met_phi, weight);
        h_mm.h3_z_m->Fill(zP4.mass(), weight);
        h_mm.h3_z_pt->Fill(zP4.pt(), weight);
        h_mm.h3_z_eta->Fill(zP4.eta(), weight);
        h_mm.h3_z_phi->Fill(zP4.phi(), weight);
        h_mm.h3_jets_n->Fill(jets_n, weight);
        h_mm.h3_jets_ht->Fill(jets_ht, weight);
        const auto& jet1 = dynamic_cast<const cat::Jet&>(*out_jets->at(0));
        h_mm.h3_jet1_m->Fill(jet1.mass(), weight);
        h_mm.h3_jet1_pt->Fill(jet1.pt(), weight);
        h_mm.h3_jet1_eta->Fill(jet1.eta(), weight);
        h_mm.h3_jet1_phi->Fill(jet1.phi(), weight);
        h_mm.h3_jet1_btag->Fill(jet1.bDiscriminator(bTagName_), weight);
        const auto& jet2 = dynamic_cast<const cat::Jet&>(*out_jets->at(1));
        h_mm.h3_jet2_m->Fill(jet2.mass(), weight);
        h_mm.h3_jet2_pt->Fill(jet2.pt(), weight);
        h_mm.h3_jet2_eta->Fill(jet2.eta(), weight);
        h_mm.h3_jet2_phi->Fill(jet2.phi(), weight);
        h_mm.h3_jet2_btag->Fill(jet2.bDiscriminator(bTagName_), weight);
        if ( jets_n >= 3 )
        {
          const auto& jet = dynamic_cast<const cat::Jet&>(*out_jets->at(2));
          h_mm.h3_jet3_m->Fill(jet.mass(), weight);
          h_mm.h3_jet3_pt->Fill(jet.pt(), weight);
          h_mm.h3_jet3_eta->Fill(jet.eta(), weight);
          h_mm.h3_jet3_phi->Fill(jet.phi(), weight);
          h_mm.h3_jet3_btag->Fill(jet.bDiscriminator(bTagName_), weight);
        }
        if ( jets_n >= 4 )
        {
          const auto& jet = dynamic_cast<const cat::Jet&>(*out_jets->at(3));
          h_mm.h3_jet4_m->Fill(jet.mass(), weight);
          h_mm.h3_jet4_pt->Fill(jet.pt(), weight);
          h_mm.h3_jet4_eta->Fill(jet.eta(), weight);
          h_mm.h3_jet4_phi->Fill(jet.phi(), weight);
          h_mm.h3_jet4_btag->Fill(jet.bDiscriminator(bTagName_), weight);
        }
        h_mm.h3_bjets_n->Fill(bjets_n, weight);
        h_mm.h3_bjets_ht->Fill(bjets_ht, weight);

        h_mm.h3_event_st->Fill(leptons_st+jets_ht+met_pt, weight);

        // Step4 Missing transverse momentum
        if ( met_pt < 40 ) break;
        ++cutstep;
        h_mm.hCutstep->Fill(cutstep, weight);
        h_mm.hCutstepNoweight->Fill(cutstep);
        h_mm.h4_vertex_n->Fill(weight);
        h_mm.h4_met_pt->Fill(met_pt, weight);
        h_mm.h4_met_phi->Fill(met_phi, weight);
        h_mm.h4_z_m->Fill(zP4.mass(), weight);
        h_mm.h4_z_pt->Fill(zP4.pt(), weight);
        h_mm.h4_z_eta->Fill(zP4.eta(), weight);
        h_mm.h4_z_phi->Fill(zP4.phi(), weight);
        h_mm.h4_jets_n->Fill(jets_n, weight);
        h_mm.h4_jets_ht->Fill(jets_ht, weight);
        h_mm.h4_jet1_m->Fill(jet1.mass(), weight);
        h_mm.h4_jet1_pt->Fill(jet1.pt(), weight);
        h_mm.h4_jet1_eta->Fill(jet1.eta(), weight);
        h_mm.h4_jet1_phi->Fill(jet1.phi(), weight);
        h_mm.h4_jet1_btag->Fill(jet1.bDiscriminator(bTagName_), weight);
        h_mm.h4_jet2_m->Fill(jet2.mass(), weight);
        h_mm.h4_jet2_pt->Fill(jet2.pt(), weight);
        h_mm.h4_jet2_eta->Fill(jet2.eta(), weight);
        h_mm.h4_jet2_phi->Fill(jet2.phi(), weight);
        h_mm.h4_jet2_btag->Fill(jet2.bDiscriminator(bTagName_), weight);
        if ( jets_n >= 3 )
        {
          const auto& jet = dynamic_cast<const cat::Jet&>(*out_jets->at(2));
          h_mm.h4_jet3_m->Fill(jet.mass(), weight);
          h_mm.h4_jet3_pt->Fill(jet.pt(), weight);
          h_mm.h4_jet3_eta->Fill(jet.eta(), weight);
          h_mm.h4_jet3_phi->Fill(jet.phi(), weight);
          h_mm.h4_jet3_btag->Fill(jet.bDiscriminator(bTagName_), weight);
        }
        if ( jets_n >= 4 )
        {
          const auto& jet = dynamic_cast<const cat::Jet&>(*out_jets->at(3));
          h_mm.h4_jet4_m->Fill(jet.mass(), weight);
          h_mm.h4_jet4_pt->Fill(jet.pt(), weight);
          h_mm.h4_jet4_eta->Fill(jet.eta(), weight);
          h_mm.h4_jet4_phi->Fill(jet.phi(), weight);
          h_mm.h4_jet4_btag->Fill(jet.bDiscriminator(bTagName_), weight);
        }
        h_mm.h4_bjets_n->Fill(bjets_n, weight);
        h_mm.h4_bjets_ht->Fill(bjets_ht, weight);

        h_mm.h4_event_st->Fill(leptons_st+jets_ht+met_pt, weight);

        // Step5 one b jet
        if ( bjets_n < 1 ) break;
        ++cutstep;
        h_mm.hCutstep->Fill(cutstep, weight);
        h_mm.hCutstepNoweight->Fill(cutstep);
        h_mm.h5_vertex_n->Fill(weight);
        h_mm.h5_met_pt->Fill(met_pt, weight);
        h_mm.h5_met_phi->Fill(met_phi, weight);
        h_mm.h5_z_m->Fill(zP4.mass(), weight);
        h_mm.h5_z_pt->Fill(zP4.pt(), weight);
        h_mm.h5_z_eta->Fill(zP4.eta(), weight);
        h_mm.h5_z_phi->Fill(zP4.phi(), weight);
        h_mm.h5_jets_n->Fill(jets_n, weight);
        h_mm.h5_jets_ht->Fill(jets_ht, weight);
        h_mm.h5_jet1_m->Fill(jet1.mass(), weight);
        h_mm.h5_jet1_pt->Fill(jet1.pt(), weight);
        h_mm.h5_jet1_eta->Fill(jet1.eta(), weight);
        h_mm.h5_jet1_phi->Fill(jet1.phi(), weight);
        h_mm.h5_jet1_btag->Fill(jet1.bDiscriminator(bTagName_), weight);
        h_mm.h5_jet2_m->Fill(jet2.mass(), weight);
        h_mm.h5_jet2_pt->Fill(jet2.pt(), weight);
        h_mm.h5_jet2_eta->Fill(jet2.eta(), weight);
        h_mm.h5_jet2_phi->Fill(jet2.phi(), weight);
        h_mm.h5_jet2_btag->Fill(jet2.bDiscriminator(bTagName_), weight);
        if ( jets_n >= 3 )
        {
          const auto& jet = dynamic_cast<const cat::Jet&>(*out_jets->at(2));
          h_mm.h5_jet3_m->Fill(jet.mass(), weight);
          h_mm.h5_jet3_pt->Fill(jet.pt(), weight);
          h_mm.h5_jet3_eta->Fill(jet.eta(), weight);
          h_mm.h5_jet3_phi->Fill(jet.phi(), weight);
          h_mm.h5_jet3_btag->Fill(jet.bDiscriminator(bTagName_), weight);
        }
        if ( jets_n >= 4 )
        {
          const auto& jet = dynamic_cast<const cat::Jet&>(*out_jets->at(3));
          h_mm.h5_jet4_m->Fill(jet.mass(), weight);
          h_mm.h5_jet4_pt->Fill(jet.pt(), weight);
          h_mm.h5_jet4_eta->Fill(jet.eta(), weight);
          h_mm.h5_jet4_phi->Fill(jet.phi(), weight);
          h_mm.h5_jet4_btag->Fill(jet.bDiscriminator(bTagName_), weight);
        }
        h_mm.h5_bjets_n->Fill(bjets_n, weight);
        h_mm.h5_bjets_ht->Fill(bjets_ht, weight);

        h_mm.h5_event_st->Fill(leptons_st+jets_ht+met_pt, weight);
        break;
      }
      ///////////////////////////////////////////
      // Finally scan cut flow for emu channel //
      ///////////////////////////////////////////
      case CH_MUEL: {
        cutstep = cutstep_em+1;

        h_em.hCutstep->Fill(cutstep, weight);
        h_em.hCutstepNoweight->Fill(cutstep);
        h_em.h1_vertex_n->Fill(weight);
        h_em.h1_met_pt->Fill(met_pt, weight);
        h_em.h1_met_phi->Fill(met_phi, weight);
        h_em.h1_leptons_n->Fill(leptons_n, weight);
        h_em.h1_lepton1_pt->Fill(lepton1->pt(), weight);
        h_em.h1_lepton1_eta->Fill(lepton1->eta(), weight);
        h_em.h1_lepton1_phi->Fill(lepton1->phi(), weight);
        h_em.h1_lepton1_q->Fill(lepton1->charge(), weight);
        h_em.h1_lepton2_pt->Fill(lepton2->pt(), weight);
        h_em.h1_lepton2_eta->Fill(lepton2->eta(), weight);
        h_em.h1_lepton2_phi->Fill(lepton2->phi(), weight);
        h_em.h1_lepton2_q->Fill(lepton2->charge(), weight);
        h_em.h1_z_m->Fill(z_m, weight);
        h_em.h1_z_pt->Fill(zP4.pt(), weight);
        h_em.h1_z_eta->Fill(zP4.eta(), weight);
        h_em.h1_z_phi->Fill(zP4.phi(), weight);
        h_em.h1_jets_n->Fill(jets_n, weight);
        h_em.h1_jets_ht->Fill(jets_ht, weight);
        if ( jets_n >= 1 )
        {
          const auto& jet = dynamic_cast<const cat::Jet&>(*out_jets->at(0));
          h_em.h1_jet1_m->Fill(jet.mass(), weight);
          h_em.h1_jet1_pt->Fill(jet.pt(), weight);
          h_em.h1_jet1_eta->Fill(jet.eta(), weight);
          h_em.h1_jet1_phi->Fill(jet.phi(), weight);
          h_em.h1_jet1_btag->Fill(jet.bDiscriminator(bTagName_), weight);
        }
        if ( jets_n >= 2 )
        {
          const auto& jet = dynamic_cast<const cat::Jet&>(*out_jets->at(1));
          h_em.h1_jet2_m->Fill(jet.mass(), weight);
          h_em.h1_jet2_pt->Fill(jet.pt(), weight);
          h_em.h1_jet2_eta->Fill(jet.eta(), weight);
          h_em.h1_jet2_phi->Fill(jet.phi(), weight);
          h_em.h1_jet2_btag->Fill(jet.bDiscriminator(bTagName_), weight);
        }
        if ( jets_n >= 3 )
        {
          const auto& jet = dynamic_cast<const cat::Jet&>(*out_jets->at(2));
          h_em.h1_jet3_m->Fill(jet.mass(), weight);
          h_em.h1_jet3_pt->Fill(jet.pt(), weight);
          h_em.h1_jet3_eta->Fill(jet.eta(), weight);
          h_em.h1_jet3_phi->Fill(jet.phi(), weight);
          h_em.h1_jet3_btag->Fill(jet.bDiscriminator(bTagName_), weight);
        }
        if ( jets_n >= 4 )
        {
          const auto& jet = dynamic_cast<const cat::Jet&>(*out_jets->at(3));
          h_em.h1_jet4_m->Fill(jet.mass(), weight);
          h_em.h1_jet4_pt->Fill(jet.pt(), weight);
          h_em.h1_jet4_eta->Fill(jet.eta(), weight);
          h_em.h1_jet4_phi->Fill(jet.phi(), weight);
          h_em.h1_jet4_btag->Fill(jet.bDiscriminator(bTagName_), weight);
        }
        h_em.h1_bjets_n->Fill(bjets_n, weight);
        h_em.h1_bjets_ht->Fill(bjets_ht, weight);

        h_em.h1_event_st->Fill(leptons_st+jets_ht+met_pt, weight);

        // Step2 Z mass veto - not applied for emu channel
        //if ( false and 76 <= z_m and z_m <= 106 ) break;
        ++cutstep;
        h_em.hCutstep->Fill(cutstep, weight);
        h_em.hCutstepNoweight->Fill(cutstep);
        h_em.h2_vertex_n->Fill(weight);
        h_em.h2_met_pt->Fill(met_pt, weight);
        h_em.h2_met_phi->Fill(met_phi, weight);
        h_em.h2_leptons_n->Fill(leptons_n, weight);
        h_em.h2_lepton1_pt->Fill(lepton1->pt(), weight);
        h_em.h2_lepton1_eta->Fill(lepton1->eta(), weight);
        h_em.h2_lepton1_phi->Fill(lepton1->phi(), weight);
        h_em.h2_lepton1_q->Fill(lepton1->charge(), weight);
        h_em.h2_lepton2_pt->Fill(lepton2->pt(), weight);
        h_em.h2_lepton2_eta->Fill(lepton2->eta(), weight);
        h_em.h2_lepton2_phi->Fill(lepton2->phi(), weight);
        h_em.h2_lepton2_q->Fill(lepton2->charge(), weight);
        h_em.h2_z_m->Fill(zP4.mass(), weight);
        h_em.h2_z_pt->Fill(zP4.pt(), weight);
        h_em.h2_z_eta->Fill(zP4.eta(), weight);
        h_em.h2_z_phi->Fill(zP4.phi(), weight);
        h_em.h2_jets_n->Fill(jets_n, weight);
        h_em.h2_jets_ht->Fill(jets_ht, weight);
        if ( jets_n >= 1 )
        {
          const auto& jet = dynamic_cast<const cat::Jet&>(*out_jets->at(0));
          h_em.h2_jet1_m->Fill(jet.mass(), weight);
          h_em.h2_jet1_pt->Fill(jet.pt(), weight);
          h_em.h2_jet1_eta->Fill(jet.eta(), weight);
          h_em.h2_jet1_phi->Fill(jet.phi(), weight);
          h_em.h2_jet1_btag->Fill(jet.bDiscriminator(bTagName_), weight);
        }
        if ( jets_n >= 2 )
        {
          const auto& jet = dynamic_cast<const cat::Jet&>(*out_jets->at(1));
          h_em.h2_jet2_m->Fill(jet.mass(), weight);
          h_em.h2_jet2_pt->Fill(jet.pt(), weight);
          h_em.h2_jet2_eta->Fill(jet.eta(), weight);
          h_em.h2_jet2_phi->Fill(jet.phi(), weight);
          h_em.h2_jet2_btag->Fill(jet.bDiscriminator(bTagName_), weight);
        }
        if ( jets_n >= 3 )
        {
          const auto& jet = dynamic_cast<const cat::Jet&>(*out_jets->at(2));
          h_em.h2_jet3_m->Fill(jet.mass(), weight);
          h_em.h2_jet3_pt->Fill(jet.pt(), weight);
          h_em.h2_jet3_eta->Fill(jet.eta(), weight);
          h_em.h2_jet3_phi->Fill(jet.phi(), weight);
          h_em.h2_jet3_btag->Fill(jet.bDiscriminator(bTagName_), weight);
        }
        if ( jets_n >= 4 )
        {
          const auto& jet = dynamic_cast<const cat::Jet&>(*out_jets->at(3));
          h_em.h2_jet4_m->Fill(jet.mass(), weight);
          h_em.h2_jet4_pt->Fill(jet.pt(), weight);
          h_em.h2_jet4_eta->Fill(jet.eta(), weight);
          h_em.h2_jet4_phi->Fill(jet.phi(), weight);
          h_em.h2_jet4_btag->Fill(jet.bDiscriminator(bTagName_), weight);
        }
        h_em.h2_bjets_n->Fill(bjets_n, weight);
        h_em.h2_bjets_ht->Fill(bjets_ht, weight);

        h_em.h2_event_st->Fill(leptons_st+jets_ht+met_pt, weight);

        // Step3 Minimal jet multiplicity
        if ( jets_n < 2 ) break;
        ++cutstep;
        h_em.hCutstep->Fill(cutstep, weight);
        h_em.hCutstepNoweight->Fill(cutstep);
        h_em.h3_vertex_n->Fill(weight);
        h_em.h3_met_pt->Fill(met_pt, weight);
        h_em.h3_met_phi->Fill(met_phi, weight);
        h_em.h3_z_m->Fill(zP4.mass(), weight);
        h_em.h3_z_pt->Fill(zP4.pt(), weight);
        h_em.h3_z_eta->Fill(zP4.eta(), weight);
        h_em.h3_z_phi->Fill(zP4.phi(), weight);
        h_em.h3_jets_n->Fill(jets_n, weight);
        h_em.h3_jets_ht->Fill(jets_ht, weight);
        const auto& jet1 = dynamic_cast<const cat::Jet&>(*out_jets->at(0));
        h_em.h3_jet1_m->Fill(jet1.mass(), weight);
        h_em.h3_jet1_pt->Fill(jet1.pt(), weight);
        h_em.h3_jet1_eta->Fill(jet1.eta(), weight);
        h_em.h3_jet1_phi->Fill(jet1.phi(), weight);
        h_em.h3_jet1_btag->Fill(jet1.bDiscriminator(bTagName_), weight);
        const auto& jet2 = dynamic_cast<const cat::Jet&>(*out_jets->at(1));
        h_em.h3_jet2_m->Fill(jet2.mass(), weight);
        h_em.h3_jet2_pt->Fill(jet2.pt(), weight);
        h_em.h3_jet2_eta->Fill(jet2.eta(), weight);
        h_em.h3_jet2_phi->Fill(jet2.phi(), weight);
        h_em.h3_jet2_btag->Fill(jet2.bDiscriminator(bTagName_), weight);
        if ( jets_n >= 3 )
        {
          const auto& jet = dynamic_cast<const cat::Jet&>(*out_jets->at(2));
          h_em.h3_jet3_m->Fill(jet.mass(), weight);
          h_em.h3_jet3_pt->Fill(jet.pt(), weight);
          h_em.h3_jet3_eta->Fill(jet.eta(), weight);
          h_em.h3_jet3_phi->Fill(jet.phi(), weight);
          h_em.h3_jet3_btag->Fill(jet.bDiscriminator(bTagName_), weight);
        }
        if ( jets_n >= 4 )
        {
          const auto& jet = dynamic_cast<const cat::Jet&>(*out_jets->at(3));
          h_em.h3_jet4_m->Fill(jet.mass(), weight);
          h_em.h3_jet4_pt->Fill(jet.pt(), weight);
          h_em.h3_jet4_eta->Fill(jet.eta(), weight);
          h_em.h3_jet4_phi->Fill(jet.phi(), weight);
          h_em.h3_jet4_btag->Fill(jet.bDiscriminator(bTagName_), weight);
        }
        h_em.h3_bjets_n->Fill(bjets_n, weight);
        h_em.h3_bjets_ht->Fill(bjets_ht, weight);

        h_em.h3_event_st->Fill(leptons_st+jets_ht+met_pt, weight);

        // Step4 Missing transverse momentum
        // Not applied for emu channel
        //if ( false and met_pt < 40 ) break;
        ++cutstep;
        h_em.hCutstep->Fill(cutstep, weight);
        h_em.hCutstepNoweight->Fill(cutstep);
        h_em.h4_vertex_n->Fill(weight);
        h_em.h4_met_pt->Fill(met_pt, weight);
        h_em.h4_met_phi->Fill(met_phi, weight);
        h_em.h4_z_m->Fill(zP4.mass(), weight);
        h_em.h4_z_pt->Fill(zP4.pt(), weight);
        h_em.h4_z_eta->Fill(zP4.eta(), weight);
        h_em.h4_z_phi->Fill(zP4.phi(), weight);
        h_em.h4_jets_n->Fill(jets_n, weight);
        h_em.h4_jets_ht->Fill(jets_ht, weight);
        h_em.h4_jet1_m->Fill(jet1.mass(), weight);
        h_em.h4_jet1_pt->Fill(jet1.pt(), weight);
        h_em.h4_jet1_eta->Fill(jet1.eta(), weight);
        h_em.h4_jet1_phi->Fill(jet1.phi(), weight);
        h_em.h4_jet1_btag->Fill(jet1.bDiscriminator(bTagName_), weight);
        h_em.h4_jet2_m->Fill(jet2.mass(), weight);
        h_em.h4_jet2_pt->Fill(jet2.pt(), weight);
        h_em.h4_jet2_eta->Fill(jet2.eta(), weight);
        h_em.h4_jet2_phi->Fill(jet2.phi(), weight);
        h_em.h4_jet2_btag->Fill(jet2.bDiscriminator(bTagName_), weight);
        if ( jets_n >= 3 )
        {
          const auto& jet = dynamic_cast<const cat::Jet&>(*out_jets->at(2));
          h_em.h4_jet3_m->Fill(jet.mass(), weight);
          h_em.h4_jet3_pt->Fill(jet.pt(), weight);
          h_em.h4_jet3_eta->Fill(jet.eta(), weight);
          h_em.h4_jet3_phi->Fill(jet.phi(), weight);
          h_em.h4_jet3_btag->Fill(jet.bDiscriminator(bTagName_), weight);
        }
        if ( jets_n >= 4 )
        {
          const auto& jet = dynamic_cast<const cat::Jet&>(*out_jets->at(3));
          h_em.h4_jet4_m->Fill(jet.mass(), weight);
          h_em.h4_jet4_pt->Fill(jet.pt(), weight);
          h_em.h4_jet4_eta->Fill(jet.eta(), weight);
          h_em.h4_jet4_phi->Fill(jet.phi(), weight);
          h_em.h4_jet4_btag->Fill(jet.bDiscriminator(bTagName_), weight);
        }
        h_em.h4_bjets_n->Fill(bjets_n, weight);
        h_em.h4_bjets_ht->Fill(bjets_ht, weight);

        h_em.h4_event_st->Fill(leptons_st+jets_ht+met_pt, weight);

        // Step5 one b jet
        if ( bjets_n < 1 ) break;
        ++cutstep;
        h_em.hCutstep->Fill(cutstep, weight);
        h_em.hCutstepNoweight->Fill(cutstep);
        h_em.h5_vertex_n->Fill(weight);
        h_em.h5_met_pt->Fill(met_pt, weight);
        h_em.h5_met_phi->Fill(met_phi, weight);
        h_em.h5_z_m->Fill(zP4.mass(), weight);
        h_em.h5_z_pt->Fill(zP4.pt(), weight);
        h_em.h5_z_eta->Fill(zP4.eta(), weight);
        h_em.h5_z_phi->Fill(zP4.phi(), weight);
        h_em.h5_jets_n->Fill(jets_n, weight);
        h_em.h5_jets_ht->Fill(jets_ht, weight);
        h_em.h5_jet1_m->Fill(jet1.mass(), weight);
        h_em.h5_jet1_pt->Fill(jet1.pt(), weight);
        h_em.h5_jet1_eta->Fill(jet1.eta(), weight);
        h_em.h5_jet1_phi->Fill(jet1.phi(), weight);
        h_em.h5_jet1_btag->Fill(jet1.bDiscriminator(bTagName_), weight);
        h_em.h5_jet2_m->Fill(jet2.mass(), weight);
        h_em.h5_jet2_pt->Fill(jet2.pt(), weight);
        h_em.h5_jet2_eta->Fill(jet2.eta(), weight);
        h_em.h5_jet2_phi->Fill(jet2.phi(), weight);
        h_em.h5_jet2_btag->Fill(jet2.bDiscriminator(bTagName_), weight);
        if ( jets_n >= 3 )
        {
          const auto& jet = dynamic_cast<const cat::Jet&>(*out_jets->at(2));
          h_em.h5_jet3_m->Fill(jet.mass(), weight);
          h_em.h5_jet3_pt->Fill(jet.pt(), weight);
          h_em.h5_jet3_eta->Fill(jet.eta(), weight);
          h_em.h5_jet3_phi->Fill(jet.phi(), weight);
          h_em.h5_jet3_btag->Fill(jet.bDiscriminator(bTagName_), weight);
        }
        if ( jets_n >= 4 )
        {
          const auto& jet = dynamic_cast<const cat::Jet&>(*out_jets->at(3));
          h_em.h5_jet4_m->Fill(jet.mass(), weight);
          h_em.h5_jet4_pt->Fill(jet.pt(), weight);
          h_em.h5_jet4_eta->Fill(jet.eta(), weight);
          h_em.h5_jet4_phi->Fill(jet.phi(), weight);
          h_em.h5_jet4_btag->Fill(jet.bDiscriminator(bTagName_), weight);
        }
        h_em.h5_bjets_n->Fill(bjets_n, weight);
        h_em.h5_bjets_ht->Fill(bjets_ht, weight);

        h_em.h5_event_st->Fill(leptons_st+jets_ht+met_pt, weight);
        break;
      }
      case CH_NONE: break;
    }
  }

  event.put(std::auto_ptr<int>(new int(channel)), "channel");
  event.put(std::auto_ptr<float>(new float(metP4.pt())), "met");
  event.put(std::auto_ptr<float>(new float(metP4.phi())), "metphi");
  event.put(out_leptons, "leptons");
  event.put(out_jets, "jets");

  return cutstep >= filterCutStepBefore_;
}

#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(TTLLEventSelector);

