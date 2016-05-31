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

struct ControlPlots
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
  H1 h1_lepton2_pt, h1_lepton2_eta, h1_lepton2_phi, h1_lepton2_q;
  H1 h1_z_m, h1_z_pt, h1_z_eta, h1_z_phi;
  H1 h1_z_m_noveto;
  H1 h1_jets_n, h1_jets_pt, h1_jets_eta, h1_jets_ht;
  H1 h1_jet1_m, h1_jet1_pt, h1_jet1_eta, h1_jet1_phi, h1_jet1_btag;
  H1 h1_jet2_m, h1_jet2_pt, h1_jet2_eta, h1_jet2_phi, h1_jet2_btag;
  H1 h1_jet3_m, h1_jet3_pt, h1_jet3_eta, h1_jet3_phi, h1_jet3_btag;
  H1 h1_jet4_m, h1_jet4_pt, h1_jet4_eta, h1_jet4_phi, h1_jet4_btag;
  H1 h1_bjets_n;
  H1 h1_event_st;

  H1 h2_vertex_n;
  H1 h2_met_pt, h2_met_phi;
  H1 h2_leptons_n;
  H1 h2_lepton1_pt, h2_lepton1_eta, h2_lepton1_phi, h2_lepton1_q;
  H1 h2_lepton2_pt, h2_lepton2_eta, h2_lepton2_phi, h2_lepton2_q;
  H1 h2_z_m, h2_z_pt, h2_z_eta, h2_z_phi;
  H1 h2_z_m_noveto;
  H1 h2_jets_n, h2_jets_pt, h2_jets_eta, h2_jets_ht;
  H1 h2_jet1_m, h2_jet1_pt, h2_jet1_eta, h2_jet1_phi, h2_jet1_btag;
  H1 h2_jet2_m, h2_jet2_pt, h2_jet2_eta, h2_jet2_phi, h2_jet2_btag;
  H1 h2_jet3_m, h2_jet3_pt, h2_jet3_eta, h2_jet3_phi, h2_jet3_btag;
  H1 h2_jet4_m, h2_jet4_pt, h2_jet4_eta, h2_jet4_phi, h2_jet4_btag;
  H1 h2_bjets_n;
  H1 h2_event_st;

  H1 h3_vertex_n;
  H1 h3_met_pt, h3_met_phi;
  H1 h3_z_m, h3_z_pt, h3_z_eta, h3_z_phi;
  H1 h3_z_m_noveto;
  H1 h3_jets_n, h3_jets_pt, h3_jets_eta, h3_jets_ht;
  H1 h3_jet1_m, h3_jet1_pt, h3_jet1_eta, h3_jet1_phi, h3_jet1_btag;
  H1 h3_jet2_m, h3_jet2_pt, h3_jet2_eta, h3_jet2_phi, h3_jet2_btag;
  H1 h3_jet3_m, h3_jet3_pt, h3_jet3_eta, h3_jet3_phi, h3_jet3_btag;
  H1 h3_jet4_m, h3_jet4_pt, h3_jet4_eta, h3_jet4_phi, h3_jet4_btag;
  H1 h3_bjets_n;
  H1 h3_event_st;

  H1 h4_vertex_n;
  H1 h4_met_pt, h4_met_phi;
  H1 h4_z_m, h4_z_pt, h4_z_eta, h4_z_phi;
  H1 h4_z_m_noveto;
  H1 h4_jets_n, h4_jets_pt, h4_jets_eta, h4_jets_ht;
  H1 h4_jet1_m, h4_jet1_pt, h4_jet1_eta, h4_jet1_phi, h4_jet1_btag;
  H1 h4_jet2_m, h4_jet2_pt, h4_jet2_eta, h4_jet2_phi, h4_jet2_btag;
  H1 h4_jet3_m, h4_jet3_pt, h4_jet3_eta, h4_jet3_phi, h4_jet3_btag;
  H1 h4_jet4_m, h4_jet4_pt, h4_jet4_eta, h4_jet4_phi, h4_jet4_btag;
  H1 h4_bjets_n;
  H1 h4_event_st;

  H1 h5_vertex_n;
  H1 h5_met_pt, h5_met_phi;
  H1 h5_z_m, h5_z_pt, h5_z_eta, h5_z_phi;
  H1 h5_z_m_noveto;
  H1 h5_jets_n, h5_jets_pt, h5_jets_eta, h5_jets_ht;
  H1 h5_jet1_m, h5_jet1_pt, h5_jet1_eta, h5_jet1_phi, h5_jet1_btag;
  H1 h5_jet2_m, h5_jet2_pt, h5_jet2_eta, h5_jet2_phi, h5_jet2_btag;
  H1 h5_jet3_m, h5_jet3_pt, h5_jet3_eta, h5_jet3_phi, h5_jet3_btag;
  H1 h5_jet4_m, h5_jet4_pt, h5_jet4_eta, h5_jet4_phi, h5_jet4_btag;
  H1 h5_bjets_n;
  H1 h5_event_st;

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
    hCutstep->GetXaxis()->SetBinLabel(4, "S1 Dilepton");
    hCutstep->GetXaxis()->SetBinLabel(5, "S2 Z veto");
    hCutstep->GetXaxis()->SetBinLabel(6, "S3 nJet2");
    hCutstep->GetXaxis()->SetBinLabel(7, "S4 MET40");
    hCutstep->GetXaxis()->SetBinLabel(8, "S5 nBJet1");

    hCutstepNoweight->GetXaxis()->SetBinLabel(1, "S0a all event");
    hCutstepNoweight->GetXaxis()->SetBinLabel(2, "S0b Trigger");
    hCutstepNoweight->GetXaxis()->SetBinLabel(3, "S0c Event filter");
    hCutstepNoweight->GetXaxis()->SetBinLabel(4, "S1 Dilepton");
    hCutstepNoweight->GetXaxis()->SetBinLabel(5, "S2 Z veto");
    hCutstepNoweight->GetXaxis()->SetBinLabel(6, "S3 nJet2");
    hCutstepNoweight->GetXaxis()->SetBinLabel(7, "S4 MET40");
    hCutstepNoweight->GetXaxis()->SetBinLabel(8, "S5 nBJet1");

    h2Cutstep->GetXaxis()->SetBinLabel(1, "S0a all event");
    h2Cutstep->GetXaxis()->SetBinLabel(2, "S0b Trigger");
    h2Cutstep->GetXaxis()->SetBinLabel(3, "S0c Event filter");
    h2Cutstep->GetXaxis()->SetBinLabel(4, "S1 Dilepton");
    h2Cutstep->GetXaxis()->SetBinLabel(5, "S2 Z veto");
    h2Cutstep->GetXaxis()->SetBinLabel(6, "S3 nJet2");
    h2Cutstep->GetXaxis()->SetBinLabel(7, "S4 MET40");
    h2Cutstep->GetXaxis()->SetBinLabel(8, "S5 nBJet1");

    h2Cutstep->GetYaxis()->SetBinLabel(1, "S0a all event");
    h2Cutstep->GetYaxis()->SetBinLabel(2, "S0b Trigger");
    h2Cutstep->GetYaxis()->SetBinLabel(3, "S0c Event filter");
    h2Cutstep->GetYaxis()->SetBinLabel(4, "S1 Dilepton");
    h2Cutstep->GetYaxis()->SetBinLabel(5, "S2 Z veto");
    h2Cutstep->GetYaxis()->SetBinLabel(6, "S3 nJet2");
    h2Cutstep->GetYaxis()->SetBinLabel(7, "S4 MET40");
    h2Cutstep->GetYaxis()->SetBinLabel(8, "S5 nBJet1");

    h2CutstepNoweight->GetXaxis()->SetBinLabel(1, "S0a all event");
    h2CutstepNoweight->GetXaxis()->SetBinLabel(2, "S0b Trigger");
    h2CutstepNoweight->GetXaxis()->SetBinLabel(3, "S0c Event filter");
    h2CutstepNoweight->GetXaxis()->SetBinLabel(4, "S1 Dilepton");
    h2CutstepNoweight->GetXaxis()->SetBinLabel(5, "S2 Z veto");
    h2CutstepNoweight->GetXaxis()->SetBinLabel(6, "S3 nJet2");
    h2CutstepNoweight->GetXaxis()->SetBinLabel(7, "S4 MET40");
    h2CutstepNoweight->GetXaxis()->SetBinLabel(8, "S5 nBJet1");

    h2CutstepNoweight->GetYaxis()->SetBinLabel(1, "S0a all event");
    h2CutstepNoweight->GetYaxis()->SetBinLabel(2, "S0b Trigger");
    h2CutstepNoweight->GetYaxis()->SetBinLabel(3, "S0c Event filter");
    h2CutstepNoweight->GetYaxis()->SetBinLabel(4, "S1 Dilepton");
    h2CutstepNoweight->GetYaxis()->SetBinLabel(5, "S2 Z veto");
    h2CutstepNoweight->GetYaxis()->SetBinLabel(6, "S3 nJet2");
    h2CutstepNoweight->GetYaxis()->SetBinLabel(7, "S4 MET40");
    h2CutstepNoweight->GetYaxis()->SetBinLabel(8, "S5 nBJet1");

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

    h1_lepton2_pt  = subdir.make<TH1D>("lepton2_pt", "lepton2_pt", 1000, 0, 1000);
    h1_lepton2_eta = subdir.make<TH1D>("lepton2_eta", "lepton2_eta", 100, -maxeta, maxeta);
    h1_lepton2_phi = subdir.make<TH1D>("lepton2_phi", "lepton2_phi", 100, -pi, pi);
    h1_lepton2_q   = subdir.make<TH1D>("lepton2_q", "lepton2_q", 3, -1.5, 1.5);

    h1_z_m   = subdir.make<TH1D>("z_m", "z_m", 1000, 0, 1000);
    h1_z_pt  = subdir.make<TH1D>("z_pt", "z_pt", 1000, 0, 1000);
    h1_z_eta = subdir.make<TH1D>("z_eta", "z_eta", 100, -maxeta, maxeta);
    h1_z_phi = subdir.make<TH1D>("z_phi", "z_phi", 100, -pi, pi);
    h1_z_m_noveto = subdir.make<TH1D>("z_m_noveto", "z_m_noveto", 1000, 0, 1000);

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

    h1_bjets_n = subdir.make<TH1D>("bjets_n", "bjets_n", 10, 0, 10);

    h1_event_st = subdir.make<TH1D>("event_st", "event_st", 1000, 0, 1000);

    subdir = dir.mkdir("step2");
    h2_vertex_n = subdir.make<TH1D>("vertex_n", "vertex_n", 100, 0, 100);
    h2_met_pt = subdir.make<TH1D>("met_pt", "met_pt", 1000, 0, 1000);
    h2_met_phi = subdir.make<TH1D>("met_phi", "met_phi", 100, -pi, pi);
    h2_leptons_n = subdir.make<TH1D>("leptons_n", "leptons_n", 10, 0, 10);

    h2_lepton1_pt  = subdir.make<TH1D>("lepton1_pt", "lepton1_pt", 1000, 0, 1000);
    h2_lepton1_eta = subdir.make<TH1D>("lepton1_eta", "lepton1_eta", 100, -maxeta, maxeta);
    h2_lepton1_phi = subdir.make<TH1D>("lepton1_phi", "lepton1_phi", 100, -pi, pi);
    h2_lepton1_q   = subdir.make<TH1D>("lepton1_q", "lepton1_q", 3, -1.5, 1.5);

    h2_lepton2_pt  = subdir.make<TH1D>("lepton2_pt", "lepton2_pt", 1000, 0, 1000);
    h2_lepton2_eta = subdir.make<TH1D>("lepton2_eta", "lepton2_eta", 100, -maxeta, maxeta);
    h2_lepton2_phi = subdir.make<TH1D>("lepton2_phi", "lepton2_phi", 100, -pi, pi);
    h2_lepton2_q   = subdir.make<TH1D>("lepton2_q", "lepton2_q", 3, -1.5, 1.5);

    h2_z_m   = subdir.make<TH1D>("z_m", "z_m", 1000, 0, 1000);
    h2_z_pt  = subdir.make<TH1D>("z_pt", "z_pt", 1000, 0, 1000);
    h2_z_eta = subdir.make<TH1D>("z_eta", "z_eta", 100, -maxeta, maxeta);
    h2_z_phi = subdir.make<TH1D>("z_phi", "z_phi", 100, -pi, pi);
    h2_z_m_noveto = subdir.make<TH1D>("z_m_noveto", "z_m_noveto", 1000, 0, 1000);

    h2_jets_n = subdir.make<TH1D>("jets_n", "jets_n", 10, 0, 10);
    h2_jets_pt  = subdir.make<TH1D>("jets_pt", "jets_pt", 1000, 0, 1000);
    h2_jets_eta = subdir.make<TH1D>("jets_eta", "jets_eta", 100, -maxeta, maxeta);
    h2_jets_ht = subdir.make<TH1D>("jets_ht", "jets_ht", 1000, 0, 1000);

    h2_jet1_m   = subdir.make<TH1D>("jet1_m", "jet1_m", 500, 0, 500);
    h2_jet1_pt  = subdir.make<TH1D>("jet1_pt", "jet1_pt", 1000, 0, 1000);
    h2_jet1_eta = subdir.make<TH1D>("jet1_eta", "jet1_eta", 100, -maxeta, maxeta);
    h2_jet1_phi = subdir.make<TH1D>("jet1_phi", "jet1_phi", 100, -pi, pi);
    h2_jet1_btag = subdir.make<TH1D>("jet1_btag", "jet1_btag", 100, 0, 1);

    h2_jet2_m   = subdir.make<TH1D>("jet2_m", "jet2_m", 500, 0, 500);
    h2_jet2_pt  = subdir.make<TH1D>("jet2_pt", "jet2_pt", 1000, 0, 1000);
    h2_jet2_eta = subdir.make<TH1D>("jet2_eta", "jet2_eta", 100, -maxeta, maxeta);
    h2_jet2_phi = subdir.make<TH1D>("jet2_phi", "jet2_phi", 100, -pi, pi);
    h2_jet2_btag = subdir.make<TH1D>("jet2_btag", "jet2_btag", 100, 0, 1);

    h2_jet3_m   = subdir.make<TH1D>("jet3_m", "jet3_m", 500, 0, 500);
    h2_jet3_pt  = subdir.make<TH1D>("jet3_pt", "jet3_pt", 1000, 0, 1000);
    h2_jet3_eta = subdir.make<TH1D>("jet3_eta", "jet3_eta", 100, -maxeta, maxeta);
    h2_jet3_phi = subdir.make<TH1D>("jet3_phi", "jet3_phi", 100, -pi, pi);
    h2_jet3_btag = subdir.make<TH1D>("jet3_btag", "jet3_btag", 100, 0, 1);

    h2_jet4_m   = subdir.make<TH1D>("jet4_m", "jet4_m", 500, 0, 500);
    h2_jet4_pt  = subdir.make<TH1D>("jet4_pt", "jet4_pt", 1000, 0, 1000);
    h2_jet4_eta = subdir.make<TH1D>("jet4_eta", "jet4_eta", 100, -maxeta, maxeta);
    h2_jet4_phi = subdir.make<TH1D>("jet4_phi", "jet4_phi", 100, -pi, pi);
    h2_jet4_btag = subdir.make<TH1D>("jet4_btag", "jet4_btag", 100, 0, 1);

    h2_bjets_n = subdir.make<TH1D>("bjets_n", "bjets_n", 10, 0, 10);

    h2_event_st = subdir.make<TH1D>("event_st", "event_st", 1000, 0, 1000);

    subdir = dir.mkdir("step3");
    h3_vertex_n = subdir.make<TH1D>("vertex_n", "vertex_n", 100, 0, 100);
    h3_met_pt = subdir.make<TH1D>("met_pt", "met_pt", 1000, 0, 1000);
    h3_met_phi = subdir.make<TH1D>("met_phi", "met_phi", 100, -pi, pi);

    h3_z_m   = subdir.make<TH1D>("z_m", "z_m", 1000, 0, 1000);
    h3_z_pt  = subdir.make<TH1D>("z_pt", "z_pt", 1000, 0, 1000);
    h3_z_eta = subdir.make<TH1D>("z_eta", "z_eta", 100, -maxeta, maxeta);
    h3_z_phi = subdir.make<TH1D>("z_phi", "z_phi", 100, -pi, pi);
    h3_z_m_noveto = subdir.make<TH1D>("z_m_noveto", "z_m_noveto", 1000, 0, 1000);

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

    h3_bjets_n = subdir.make<TH1D>("bjets_n", "bjets_n", 10, 0, 10);

    h3_event_st = subdir.make<TH1D>("event_st", "event_st", 1000, 0, 1000);

    subdir = dir.mkdir("step4");
    h4_vertex_n = subdir.make<TH1D>("vertex_n", "vertex_n", 100, 0, 100);
    h4_met_pt = subdir.make<TH1D>("met_pt", "met_pt", 1000, 0, 1000);
    h4_met_phi = subdir.make<TH1D>("met_phi", "met_phi", 100, -pi, pi);

    h4_z_m   = subdir.make<TH1D>("z_m", "z_m", 1000, 0, 1000);
    h4_z_pt  = subdir.make<TH1D>("z_pt", "z_pt", 1000, 0, 1000);
    h4_z_eta = subdir.make<TH1D>("z_eta", "z_eta", 100, -maxeta, maxeta);
    h4_z_phi = subdir.make<TH1D>("z_phi", "z_phi", 100, -pi, pi);
    h4_z_m_noveto = subdir.make<TH1D>("z_m_noveto", "z_m_noveto", 1000, 0, 1000);

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

    h4_bjets_n = subdir.make<TH1D>("bjets_n", "bjets_n", 10, 0, 10);

    h4_event_st = subdir.make<TH1D>("event_st", "event_st", 1000, 0, 1000);

    subdir = dir.mkdir("step5");
    h5_vertex_n = subdir.make<TH1D>("vertex_n", "vertex_n", 100, 0, 100);
    h5_met_pt = subdir.make<TH1D>("met_pt", "met_pt", 1000, 0, 1000);
    h5_met_phi = subdir.make<TH1D>("met_phi", "met_phi", 100, -pi, pi);

    h5_z_m   = subdir.make<TH1D>("z_m", "z_m", 1000, 0, 1000);
    h5_z_pt  = subdir.make<TH1D>("z_pt", "z_pt", 1000, 0, 1000);
    h5_z_eta = subdir.make<TH1D>("z_eta", "z_eta", 100, -maxeta, maxeta);
    h5_z_phi = subdir.make<TH1D>("z_phi", "z_phi", 100, -pi, pi);
    h5_z_m_noveto = subdir.make<TH1D>("z_m_noveto", "z_m_noveto", 1000, 0, 1000);

    h5_jets_n = subdir.make<TH1D>("jets_n", "jets_n", 10, 0, 10);
    h5_jets_pt  = subdir.make<TH1D>("jets_pt", "jets_pt", 1000, 0, 1000);
    h5_jets_eta = subdir.make<TH1D>("jets_eta", "jets_eta", 100, -maxeta, maxeta);
    h5_jets_ht = subdir.make<TH1D>("jets_ht", "jets_ht", 1000, 0, 1000);

    h5_jet1_m   = subdir.make<TH1D>("jet1_m", "jet1_m", 500, 0, 500);
    h5_jet1_pt  = subdir.make<TH1D>("jet1_pt", "jet1_pt", 1000, 0, 1000);
    h5_jet1_eta = subdir.make<TH1D>("jet1_eta", "jet1_eta", 100, -maxeta, maxeta);
    h5_jet1_phi = subdir.make<TH1D>("jet1_phi", "jet1_phi", 100, -pi, pi);
    h5_jet1_btag = subdir.make<TH1D>("jet1_btag", "jet1_btag", 100, 0, 1);

    h5_jet2_m   = subdir.make<TH1D>("jet2_m", "jet2_m", 500, 0, 500);
    h5_jet2_pt  = subdir.make<TH1D>("jet2_pt", "jet2_pt", 1000, 0, 1000);
    h5_jet2_eta = subdir.make<TH1D>("jet2_eta", "jet2_eta", 100, -maxeta, maxeta);
    h5_jet2_phi = subdir.make<TH1D>("jet2_phi", "jet2_phi", 100, -pi, pi);
    h5_jet2_btag = subdir.make<TH1D>("jet2_btag", "jet2_btag", 100, 0, 1);

    h5_jet3_m   = subdir.make<TH1D>("jet3_m", "jet3_m", 500, 0, 500);
    h5_jet3_pt  = subdir.make<TH1D>("jet3_pt", "jet3_pt", 1000, 0, 1000);
    h5_jet3_eta = subdir.make<TH1D>("jet3_eta", "jet3_eta", 100, -maxeta, maxeta);
    h5_jet3_phi = subdir.make<TH1D>("jet3_phi", "jet3_phi", 100, -pi, pi);
    h5_jet3_btag = subdir.make<TH1D>("jet3_btag", "jet3_btag", 100, 0, 1);

    h5_jet4_m   = subdir.make<TH1D>("jet4_m", "jet4_m", 500, 0, 500);
    h5_jet4_pt  = subdir.make<TH1D>("jet4_pt", "jet4_pt", 1000, 0, 1000);
    h5_jet4_eta = subdir.make<TH1D>("jet4_eta", "jet4_eta", 100, -maxeta, maxeta);
    h5_jet4_phi = subdir.make<TH1D>("jet4_phi", "jet4_phi", 100, -pi, pi);
    h5_jet4_btag = subdir.make<TH1D>("jet4_btag", "jet4_btag", 100, 0, 1);

    h5_bjets_n = subdir.make<TH1D>("bjets_n", "bjets_n", 10, 0, 10);

    h5_event_st = subdir.make<TH1D>("event_st", "event_st", 1000, 0, 1000);

  };
};
const int ControlPlots::nMaxCutstep = 8; // 5+3

class TTLLEventSelector : public edm::one::EDFilter<edm::one::SharedResources>
{
public:
  TTLLEventSelector(const edm::ParameterSet& pset);
  bool filter(edm::Event& event, const edm::EventSetup&) override;
  ~TTLLEventSelector();

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
  edm::EDGetTokenT<int> trigElElToken_, trigMuMuToken_, trigMuElToken_;
  edm::EDGetTokenT<int> nVertexToken_;

  std::vector<edm::EDGetTokenT<float> > extWeightTokensF_;
  std::vector<edm::EDGetTokenT<double> > extWeightTokensD_;

private:
  TH1D* h_weight, * h_pileupWeight, * h_genWeight;
  ControlPlots h_ee, h_mm, h_em;

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
    if ( std::abs(mu.eta()) > 2.4 ) return false;
    if ( shiftedMuonPt(mu) < 20 ) return false;

    if ( mu.relIso(0.4) > 0.15 ) return false;
    if ( !mu.isTightMuon() ) return false;
    return true;
  }
  bool isGoodElectron(const cat::Electron& el)
  {
    if ( std::abs(el.eta()) > 2.4 ) return false;
    if ( shiftedElectronPt(el) < 20 ) return false;

    if ( isMVAElectronSel_ and !el.isTrigMVAValid() ) return false;

    //if ( el.relIso(0.3) >= 0.11 ) return false;
    if ( !el.electronID(elIdName_) ) return false;
    //if ( !el.isPF() or !el.passConversionVeto() ) return false;
    const double scEta = std::abs(el.scEta());
    if ( scEta > 1.4442 and scEta < 1.566 ) return false;
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

  bool isMC_;
  bool isIgnoreTrig_; // Accept event even if it does not pass HLT. Needed for synchronization
  const int applyFilterAt_;

  // ID variables
  bool isMVAElectronSel_;
  std::string bTagName_;
  std::string elIdName_;
  enum class BTagWP { CSVL, CSVM, CSVT } bTagWP_;

};

}

using namespace cat;

TTLLEventSelector::TTLLEventSelector(const edm::ParameterSet& pset):
  isMC_(pset.getParameter<bool>("isMC")),
  applyFilterAt_(pset.getParameter<int>("applyFilterAt"))
{
  const auto muonSet = pset.getParameter<edm::ParameterSet>("muon");
  muonToken_ = consumes<cat::MuonCollection>(muonSet.getParameter<edm::InputTag>("src"));
  muonScale_ = muonSet.getParameter<int>("scaleDirection");
  if ( isMC_ )
  {
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
  electronScale_ = electronSet.getParameter<int>("scaleDirection");
  if ( isMC_ )
  {
    const auto electronSFSet = electronSet.getParameter<edm::ParameterSet>("efficiencySF");
    // FIXME : for electrons, eta bins are NOT folded - always double check this with cfg
    electronSF_.set(electronSFSet.getParameter<vdouble>("pt_bins"),
                    electronSFSet.getParameter<vdouble>("abseta_bins"),
                    electronSFSet.getParameter<vdouble>("values"),
                    electronSFSet.getParameter<vdouble>("errors"));
    electronSFShift_ = electronSet.getParameter<int>("efficiencySFDirection");
  }
  isMVAElectronSel_ = false;
  if ( elIdName_.substr(0,3) == "mva" ) isMVAElectronSel_ = true;

  const auto jetSet = pset.getParameter<edm::ParameterSet>("jet");
  jetToken_ = consumes<cat::JetCollection>(jetSet.getParameter<edm::InputTag>("src"));
  jetScale_ = jetSet.getParameter<int>("scaleDirection");
  jetResol_ = jetSet.getParameter<int>("resolDirection");
  bTagName_ = jetSet.getParameter<string>("bTagName");
  const auto bTagWPStr = jetSet.getParameter<string>("bTagWP");
  if      ( bTagWPStr == "CSVL" ) bTagWP_ = BTagWP::CSVL;
  else if ( bTagWPStr == "CSVM" ) bTagWP_ = BTagWP::CSVM;
  else if ( bTagWPStr == "CSVT" ) bTagWP_ = BTagWP::CSVT;
  else edm::LogError("TTLLEventSelector") << "Wrong bTagWP parameter " << bTagWPStr;
  isSkipJER_ = jetSet.getParameter<bool>("skipJER");

  const auto metSet = pset.getParameter<edm::ParameterSet>("met");
  metToken_ = consumes<cat::METCollection>(metSet.getParameter<edm::InputTag>("src"));

  const auto vertexSet = pset.getParameter<edm::ParameterSet>("vertex");
  nVertexToken_ = consumes<int>(vertexSet.getParameter<edm::InputTag>("nVertex"));
  pileupWeightToken_ = consumes<float>(vertexSet.getParameter<edm::InputTag>("pileupWeight"));

  const auto filterSet = pset.getParameter<edm::ParameterSet>("filters");
  recoFilterToken_ = consumes<int>(filterSet.getParameter<edm::InputTag>("filterRECO"));
  trigElElToken_ = consumes<int>(filterSet.getParameter<edm::InputTag>("trigELEL"));
  trigMuMuToken_ = consumes<int>(filterSet.getParameter<edm::InputTag>("trigMUMU"));
  trigMuElToken_ = consumes<int>(filterSet.getParameter<edm::InputTag>("trigMUEL"));
  isIgnoreTrig_ = filterSet.getParameter<bool>("ignoreTrig");

  if ( isMC_ )
  {
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
  if ( isMC_ )
  {
    h_genWeight = doverall.make<TH1D>("genWeight", "genWeight", 200, -10, 10);
    h_pileupWeight = doverall.make<TH1D>("pileupWeight", "pileupWeight", 200, -10, 10);
  }

  h_ee.book(fs->mkdir("ee"));
  h_mm.book(fs->mkdir("mm"));
  h_em.book(fs->mkdir("em"));

  produces<int>("channel");
  produces<float>("weight");
  produces<float>("met");
  produces<float>("metphi");
  produces<std::vector<cat::Lepton> >("leptons");
  produces<std::vector<cat::Jet> >("jets");
}

bool TTLLEventSelector::filter(edm::Event& event, const edm::EventSetup&)
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
  if ( isMC_ )
  {
    float genWeight = 1.;
    edm::Handle<float> fHandle;
    edm::Handle<vfloat> vfHandle;

    if ( genWeightIndex_ < 0 ) {
      event.getByToken(genWeightToken_, vfHandle);
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

  event.getByToken(trigElElToken_, trigHandle);
  const int isTrigElEl = *trigHandle;
  event.getByToken(trigMuMuToken_, trigHandle);
  const int isTrigMuMu = *trigHandle;
  event.getByToken(trigMuElToken_, trigHandle);
  const int isTrigMuEl = *trigHandle;

  // Select good leptons
  double leptons_st = 0;
  cat::MuonCollection selMuons;
  for ( int i=0, n=muonHandle->size(); i<n; ++i )
  {
    auto& p = muonHandle->at(i);
    if ( !isGoodMuon(p) ) continue;
    const double pt = shiftedMuonPt(p);
    const double scale = pt/p.pt();

    cat::Muon lep(p);
    lep.setP4(p.p4()*scale);
    selMuons.push_back(lep);

    leptons_st += pt;
    metDpx += lep.px()-p.px();
    metDpy += lep.py()-p.py();
  }
  cat::ElectronCollection selElectrons;
  for ( int i=0, n=electronHandle->size(); i<n; ++i )
  {
    auto& p = electronHandle->at(i);
    if ( !isGoodElectron(p) ) continue;
    const double pt = shiftedElectronPt(p);
    const double scale = pt/p.pt();

    cat::Electron lep(p);
    lep.setP4(p.p4()*scale);
    selElectrons.push_back(lep);

    leptons_st += pt;
    metDpx += lep.px()-p.px();
    metDpy += lep.py()-p.py();
  }
  std::vector<const cat::Lepton*> selLeptons;
  for ( auto& x : selMuons ) selLeptons.push_back(&x);
  for ( auto& x : selElectrons ) selLeptons.push_back(&x);
  std::sort(selLeptons.begin(), selLeptons.end(),
            [&](const cat::Lepton* a, const cat::Lepton* b){return a->pt() > b->pt();});
  // Copy selLeptons to out_leptons
  for ( auto x : selLeptons ) out_leptons->push_back(*x);
  const int leptons_n = selLeptons.size();
  const cat::Lepton* lepton1 = 0, * lepton2 = 0;
  int channel = CH_NOLL;
  if ( leptons_n >= 2 ) {
    // Set lepton1 and 2
    lepton1 = selLeptons.at(0);
    lepton2 = selLeptons.at(1);

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
    // Apply lepton SF
    if ( channel == CH_ELEL )
    {
      const auto e1 = dynamic_cast<const cat::Electron*>(lepton1);
      const auto e2 = dynamic_cast<const cat::Electron*>(lepton2);
      const double w1 = electronSF_(lepton1->pt(), std::abs(e1->scEta()), electronSFShift_);
      const double w2 = electronSF_(lepton2->pt(), std::abs(e2->scEta()), electronSFShift_);
      weight *= w1*w2;
    }
    else if ( channel == CH_MUMU )
    {
      const double w1 = muonSF_(lepton1->pt(), std::abs(lepton1->eta()), muonSFShift_);
      const double w2 = muonSF_(lepton2->pt(), std::abs(lepton2->eta()), muonSFShift_);
      weight *= w1*w2;
    }
    else if ( channel == CH_MUEL )
    {
      const auto e1 = dynamic_cast<const cat::Electron*>(lepton1);
      const double w1 = electronSF_(lepton1->pt(), std::abs(e1->scEta()), electronSFShift_);
      const double w2 = muonSF_(lepton2->pt(), std::abs(lepton2->eta()), muonSFShift_);
      weight *= w1*w2;
    }
    else edm::LogError("TTLLEventSelector") << "Strange event with nLepton >=2 but not falling info ee,mumu,emu category";
  }
  const double z_m = leptons_n < 2 ? -1 : (lepton1->p4()+lepton2->p4()).mass();

  // Select good jets
  int bjets_n = 0;
  double jets_ht = 0;
  for ( int i=0, n=jetHandle->size(); i<n; ++i )
  {
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
    if ( leptons_n >= 2 and deltaR(jet.p4(), out_leptons->at(1).p4()) < 0.4 ) continue;

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

  h_ee.hCutstep->Fill(-2, weight);
  h_ee.hCutstepNoweight->Fill(-2);
  h_ee.h0a_vertex_n->Fill(nVertex, weight);

  h_mm.hCutstep->Fill(-2, weight);
  h_mm.hCutstepNoweight->Fill(-2);
  h_mm.h0a_vertex_n->Fill(nVertex, weight);

  h_em.hCutstep->Fill(-2, weight);
  h_em.hCutstepNoweight->Fill(-2);
  h_em.h0a_vertex_n->Fill(nVertex, weight);

  // ElEl channel Cutstep 0b with trigger requirements
  int cutstep_ee = -2;
  if ( isIgnoreTrig_ or isTrigElEl )
  {
    ++cutstep_ee;
    h_ee.hCutstep->Fill(-1, weight);
    h_ee.hCutstepNoweight->Fill(-1);
    h_ee.h0b_vertex_n->Fill(nVertex, weight);
    h_ee.h0b_met_pt->Fill(met_pt, weight);
    h_ee.h0b_met_phi->Fill(met_phi, weight);
    h_ee.h0b_leptons_n->Fill(leptons_n, weight);
    h_ee.h0b_jets_n->Fill(jets_n, weight);
    h_ee.h0b_bjets_n->Fill(bjets_n, weight);
    h_ee.h0b_jets_ht->Fill(jets_ht, weight);
    for ( auto jet : *out_jets ) {
      h_ee.h0b_jets_pt->Fill(jet.pt(), weight);
      h_ee.h0b_jets_eta->Fill(jet.eta(), weight);
    }

    // Cutstep 0c with reco filters
    if ( isMC_ or isRECOFilterOK )
    {
      ++cutstep_ee;
      h_ee.hCutstep->Fill(0., weight);
      h_ee.hCutstepNoweight->Fill(0.);
      h_ee.h0c_vertex_n->Fill(nVertex, weight);
      h_ee.h0c_met_pt->Fill(met_pt, weight);
      h_ee.h0c_met_phi->Fill(met_phi, weight);
      h_ee.h0c_leptons_n->Fill(leptons_n, weight);
      h_ee.h0c_jets_n->Fill(jets_n, weight);
      h_ee.h0c_bjets_n->Fill(bjets_n, weight);
      h_ee.h0c_jets_ht->Fill(jets_ht, weight);
      for ( auto jet : *out_jets ) {
        h_ee.h0c_jets_pt->Fill(jet.pt(), weight);
        h_ee.h0c_jets_eta->Fill(jet.eta(), weight);
      }
    }
  }
  // MuMu channel Cutstep 0b with trigger requirements
  int cutstep_mm = -2;
  if ( isIgnoreTrig_ or isTrigMuMu )
  {
    ++cutstep_mm;
    h_mm.hCutstep->Fill(-1, weight);
    h_mm.hCutstepNoweight->Fill(-1);
    h_mm.h0b_vertex_n->Fill(nVertex, weight);
    h_mm.h0b_met_pt->Fill(met_pt, weight);
    h_mm.h0b_met_phi->Fill(met_phi, weight);
    h_mm.h0b_leptons_n->Fill(leptons_n, weight);
    h_mm.h0b_jets_n->Fill(jets_n, weight);
    h_mm.h0b_bjets_n->Fill(bjets_n, weight);
    h_mm.h0b_jets_ht->Fill(jets_ht, weight);
    for ( auto jet : *out_jets ) {
      h_mm.h0b_jets_pt->Fill(jet.pt(), weight);
      h_mm.h0b_jets_eta->Fill(jet.eta(), weight);
    }

    // Cutstep 0c with reco filters
    if ( isMC_ or isRECOFilterOK )
    {
      ++cutstep_mm;
      h_mm.hCutstep->Fill(0., weight);
      h_mm.hCutstepNoweight->Fill(0.);
      h_mm.h0c_vertex_n->Fill(nVertex, weight);
      h_mm.h0c_met_pt->Fill(met_pt, weight);
      h_mm.h0c_met_phi->Fill(met_phi, weight);
      h_mm.h0c_leptons_n->Fill(leptons_n, weight);
      h_mm.h0c_jets_n->Fill(jets_n, weight);
      h_mm.h0c_bjets_n->Fill(bjets_n, weight);
      h_mm.h0c_jets_ht->Fill(jets_ht, weight);
      for ( auto jet : *out_jets ) {
        h_mm.h0c_jets_pt->Fill(jet.pt(), weight);
        h_mm.h0c_jets_eta->Fill(jet.eta(), weight);
      }
    }
  }
  // MuEl channel Cutstep 0b with trigger requirements
  int cutstep_em = -2;
  if ( isIgnoreTrig_ or isTrigMuEl )
  {
    ++cutstep_em;
    h_em.hCutstep->Fill(-1, weight);
    h_em.hCutstepNoweight->Fill(-1);
    h_em.h0b_vertex_n->Fill(nVertex, weight);
    h_em.h0b_met_pt->Fill(met_pt, weight);
    h_em.h0b_met_phi->Fill(met_phi, weight);
    h_em.h0b_leptons_n->Fill(leptons_n, weight);
    h_em.h0b_jets_n->Fill(jets_n, weight);
    h_em.h0b_bjets_n->Fill(bjets_n, weight);
    h_em.h0b_jets_ht->Fill(jets_ht, weight);
    for ( auto jet : *out_jets ) {
      h_em.h0b_jets_pt->Fill(jet.pt(), weight);
      h_em.h0b_jets_eta->Fill(jet.eta(), weight);
    }

    // Cutstep 0c with reco filters
    if ( isMC_ or isRECOFilterOK )
    {
      ++cutstep_em;
      h_em.hCutstep->Fill(0., weight);
      h_em.hCutstepNoweight->Fill(0.);
      h_em.h0c_vertex_n->Fill(nVertex, weight);
      h_em.h0c_met_pt->Fill(met_pt, weight);
      h_em.h0c_met_phi->Fill(met_phi, weight);
      h_em.h0c_leptons_n->Fill(leptons_n, weight);
      h_em.h0c_jets_n->Fill(jets_n, weight);
      h_em.h0c_bjets_n->Fill(bjets_n, weight);
      h_em.h0c_jets_ht->Fill(jets_ht, weight);
      for ( auto jet : *out_jets ) {
        h_em.h0b_jets_pt->Fill(jet.pt(), weight);
        h_em.h0b_jets_eta->Fill(jet.eta(), weight);
      }
    }
  }

  // Check each cut steps
  int cutstep = -1;
  // bitset for the cut steps, fill the results only for events that pass step0a,0b,0c
  std::bitset<ControlPlots::nMaxCutstep-2> cutstepBits(0);
  //for ( auto x : cutstepBits ) x = false;
  if ( (channel == CH_ELEL and cutstep_ee == 0) or
       (channel == CH_MUMU and cutstep_mm == 0) or
       (channel == CH_MUEL and cutstep_em == 0) )
  {
    // Step1 Dilepton
    if ( leptons_n >= 2 and z_m >= 20 and lepton1->charge()+lepton2->charge() == 0 )
    {
      cutstepBits[0] = true;
      // Step2 Z mass veto : Step1 have to be required by construction
      if ( channel == CH_MUEL or !(76 <= z_m and z_m <= 106) ) cutstepBits[1] = true;
    }
    // Step3 Minimal jet multiplicity
    if ( jets_n >= 2 ) cutstepBits[2] = true;
    // Step4 Missing transverse momentum
    if ( channel == CH_MUEL or met_pt >= 40 ) cutstepBits[3] = true;
    // Step5 one b jet
    if ( bjets_n >= 1 ) cutstepBits[4] = true;

    // Set the cut step of this event
    const int nstep = cutstepBits.size();
    for ( cutstep=0; cutstep<nstep; ++cutstep )
    {
      if ( !cutstepBits[cutstep] ) break;
    }

  }
  else
  {
    cutstep = std::max(cutstep_ee, std::max(cutstep_mm, cutstep_em)); // reset the cut step
  }

  // Cut step is ready. Now proceed to fill histograms
  switch (1) default: {
    int icutstep = 0;

    if ( cutstep <= 0 ) break;
    ++icutstep; // =1

    const auto lepton1P4 = shiftedLepPt(*lepton1)/lepton1->pt()*lepton1->p4();
    const auto lepton2P4 = shiftedLepPt(*lepton2)/lepton2->pt()*lepton2->p4();
    const auto zP4 = lepton1P4+lepton2P4;

    if ( channel == CH_ELEL )
    {
      h_ee.hCutstep->Fill(icutstep, weight);
      h_ee.hCutstepNoweight->Fill(icutstep);
      h_ee.h1_vertex_n->Fill(nVertex, weight);
      h_ee.h1_met_pt->Fill(met_pt, weight);
      h_ee.h1_met_phi->Fill(met_phi, weight);
      h_ee.h1_leptons_n->Fill(leptons_n, weight);
      h_ee.h1_lepton1_pt->Fill(lepton1P4.pt(), weight);
      h_ee.h1_lepton1_eta->Fill(lepton1->eta(), weight);
      h_ee.h1_lepton1_phi->Fill(lepton1->phi(), weight);
      h_ee.h1_lepton1_q->Fill(lepton1->charge(), weight);
      h_ee.h1_lepton2_pt->Fill(lepton2P4.pt(), weight);
      h_ee.h1_lepton2_eta->Fill(lepton2->eta(), weight);
      h_ee.h1_lepton2_phi->Fill(lepton2->phi(), weight);
      h_ee.h1_lepton2_q->Fill(lepton2->charge(), weight);
      h_ee.h1_z_m->Fill(z_m, weight);
      h_ee.h1_z_pt->Fill(zP4.pt(), weight);
      h_ee.h1_z_eta->Fill(zP4.eta(), weight);
      h_ee.h1_z_phi->Fill(zP4.phi(), weight);
      h_ee.h1_jets_n->Fill(jets_n, weight);
      h_ee.h1_jets_ht->Fill(jets_ht, weight);
      for ( auto jet : *out_jets ) {
        h_ee.h1_jets_pt->Fill(jet.pt(), weight);
        h_ee.h1_jets_eta->Fill(jet.eta(), weight);
      }
      if ( jets_n >= 1 )
      {
        const auto& jet = out_jets->at(0);
        h_ee.h1_jet1_m->Fill(jet.mass(), weight);
        h_ee.h1_jet1_pt->Fill(jet.pt(), weight);
        h_ee.h1_jet1_eta->Fill(jet.eta(), weight);
        h_ee.h1_jet1_phi->Fill(jet.phi(), weight);
        h_ee.h1_jet1_btag->Fill(jet.bDiscriminator(bTagName_), weight);
      }
      if ( jets_n >= 2 )
      {
        const auto& jet = out_jets->at(1);
        h_ee.h1_jet2_m->Fill(jet.mass(), weight);
        h_ee.h1_jet2_pt->Fill(jet.pt(), weight);
        h_ee.h1_jet2_eta->Fill(jet.eta(), weight);
        h_ee.h1_jet2_phi->Fill(jet.phi(), weight);
        h_ee.h1_jet2_btag->Fill(jet.bDiscriminator(bTagName_), weight);
      }
      if ( jets_n >= 3 )
      {
        const auto& jet = out_jets->at(2);
        h_ee.h1_jet3_m->Fill(jet.mass(), weight);
        h_ee.h1_jet3_pt->Fill(jet.pt(), weight);
        h_ee.h1_jet3_eta->Fill(jet.eta(), weight);
        h_ee.h1_jet3_phi->Fill(jet.phi(), weight);
        h_ee.h1_jet3_btag->Fill(jet.bDiscriminator(bTagName_), weight);
      }
      if ( jets_n >= 4 )
      {
        const auto& jet = out_jets->at(3);
        h_ee.h1_jet4_m->Fill(jet.mass(), weight);
        h_ee.h1_jet4_pt->Fill(jet.pt(), weight);
        h_ee.h1_jet4_eta->Fill(jet.eta(), weight);
        h_ee.h1_jet4_phi->Fill(jet.phi(), weight);
        h_ee.h1_jet4_btag->Fill(jet.bDiscriminator(bTagName_), weight);
      }
      h_ee.h1_bjets_n->Fill(bjets_n, weight);

      h_ee.h1_event_st->Fill(leptons_st+jets_ht+met_pt, weight);
    }
    else if ( channel == CH_MUMU )
    {
      h_mm.hCutstep->Fill(icutstep, weight);
      h_mm.hCutstepNoweight->Fill(icutstep);
      h_mm.h1_vertex_n->Fill(nVertex, weight);
      h_mm.h1_met_pt->Fill(met_pt, weight);
      h_mm.h1_met_phi->Fill(met_phi, weight);
      h_mm.h1_leptons_n->Fill(leptons_n, weight);
      h_mm.h1_lepton1_pt->Fill(lepton1P4.pt(), weight);
      h_mm.h1_lepton1_eta->Fill(lepton1->eta(), weight);
      h_mm.h1_lepton1_phi->Fill(lepton1->phi(), weight);
      h_mm.h1_lepton1_q->Fill(lepton1->charge(), weight);
      h_mm.h1_lepton2_pt->Fill(lepton2P4.pt(), weight);
      h_mm.h1_lepton2_eta->Fill(lepton2->eta(), weight);
      h_mm.h1_lepton2_phi->Fill(lepton2->phi(), weight);
      h_mm.h1_lepton2_q->Fill(lepton2->charge(), weight);
      h_mm.h1_z_m->Fill(z_m, weight);
      h_mm.h1_z_pt->Fill(zP4.pt(), weight);
      h_mm.h1_z_eta->Fill(zP4.eta(), weight);
      h_mm.h1_z_phi->Fill(zP4.phi(), weight);
      h_mm.h1_jets_n->Fill(jets_n, weight);
      h_mm.h1_jets_ht->Fill(jets_ht, weight);
      for ( auto jet : *out_jets ) {
        h_mm.h1_jets_pt->Fill(jet.pt(), weight);
        h_mm.h1_jets_eta->Fill(jet.eta(), weight);
      }
      if ( jets_n >= 1 )
      {
        const auto& jet = out_jets->at(0);
        h_mm.h1_jet1_m->Fill(jet.mass(), weight);
        h_mm.h1_jet1_pt->Fill(jet.pt(), weight);
        h_mm.h1_jet1_eta->Fill(jet.eta(), weight);
        h_mm.h1_jet1_phi->Fill(jet.phi(), weight);
        h_mm.h1_jet1_btag->Fill(jet.bDiscriminator(bTagName_), weight);
      }
      if ( jets_n >= 2 )
      {
        const auto& jet = out_jets->at(1);
        h_mm.h1_jet2_m->Fill(jet.mass(), weight);
        h_mm.h1_jet2_pt->Fill(jet.pt(), weight);
        h_mm.h1_jet2_eta->Fill(jet.eta(), weight);
        h_mm.h1_jet2_phi->Fill(jet.phi(), weight);
        h_mm.h1_jet2_btag->Fill(jet.bDiscriminator(bTagName_), weight);
      }
      if ( jets_n >= 3 )
      {
        const auto& jet = out_jets->at(2);
        h_mm.h1_jet3_m->Fill(jet.mass(), weight);
        h_mm.h1_jet3_pt->Fill(jet.pt(), weight);
        h_mm.h1_jet3_eta->Fill(jet.eta(), weight);
        h_mm.h1_jet3_phi->Fill(jet.phi(), weight);
        h_mm.h1_jet3_btag->Fill(jet.bDiscriminator(bTagName_), weight);
      }
      if ( jets_n >= 4 )
      {
        const auto& jet = out_jets->at(3);
        h_mm.h1_jet4_m->Fill(jet.mass(), weight);
        h_mm.h1_jet4_pt->Fill(jet.pt(), weight);
        h_mm.h1_jet4_eta->Fill(jet.eta(), weight);
        h_mm.h1_jet4_phi->Fill(jet.phi(), weight);
        h_mm.h1_jet4_btag->Fill(jet.bDiscriminator(bTagName_), weight);
      }
      h_mm.h1_bjets_n->Fill(bjets_n, weight);

      h_mm.h1_event_st->Fill(leptons_st+jets_ht+met_pt, weight);
    }
    else if ( channel == CH_MUEL )
    {
      h_em.hCutstep->Fill(icutstep, weight);
      h_em.hCutstepNoweight->Fill(icutstep);
      h_em.h1_vertex_n->Fill(nVertex, weight);
      h_em.h1_met_pt->Fill(met_pt, weight);
      h_em.h1_met_phi->Fill(met_phi, weight);
      h_em.h1_leptons_n->Fill(leptons_n, weight);
      h_em.h1_lepton1_pt->Fill(lepton1P4.pt(), weight);
      h_em.h1_lepton1_eta->Fill(lepton1->eta(), weight);
      h_em.h1_lepton1_phi->Fill(lepton1->phi(), weight);
      h_em.h1_lepton1_q->Fill(lepton1->charge(), weight);
      h_em.h1_lepton2_pt->Fill(lepton2P4.pt(), weight);
      h_em.h1_lepton2_eta->Fill(lepton2->eta(), weight);
      h_em.h1_lepton2_phi->Fill(lepton2->phi(), weight);
      h_em.h1_lepton2_q->Fill(lepton2->charge(), weight);
      h_em.h1_z_m->Fill(z_m, weight);
      h_em.h1_z_pt->Fill(zP4.pt(), weight);
      h_em.h1_z_eta->Fill(zP4.eta(), weight);
      h_em.h1_z_phi->Fill(zP4.phi(), weight);
      h_em.h1_jets_n->Fill(jets_n, weight);
      h_em.h1_jets_ht->Fill(jets_ht, weight);
      for ( auto jet : *out_jets ) {
        h_em.h1_jets_pt->Fill(jet.pt(), weight);
        h_em.h1_jets_eta->Fill(jet.eta(), weight);
      }
      if ( jets_n >= 1 )
      {
        const auto& jet = out_jets->at(0);
        h_em.h1_jet1_m->Fill(jet.mass(), weight);
        h_em.h1_jet1_pt->Fill(jet.pt(), weight);
        h_em.h1_jet1_eta->Fill(jet.eta(), weight);
        h_em.h1_jet1_phi->Fill(jet.phi(), weight);
        h_em.h1_jet1_btag->Fill(jet.bDiscriminator(bTagName_), weight);
      }
      if ( jets_n >= 2 )
      {
        const auto& jet = out_jets->at(1);
        h_em.h1_jet2_m->Fill(jet.mass(), weight);
        h_em.h1_jet2_pt->Fill(jet.pt(), weight);
        h_em.h1_jet2_eta->Fill(jet.eta(), weight);
        h_em.h1_jet2_phi->Fill(jet.phi(), weight);
        h_em.h1_jet2_btag->Fill(jet.bDiscriminator(bTagName_), weight);
      }
      if ( jets_n >= 3 )
      {
        const auto& jet = out_jets->at(2);
        h_em.h1_jet3_m->Fill(jet.mass(), weight);
        h_em.h1_jet3_pt->Fill(jet.pt(), weight);
        h_em.h1_jet3_eta->Fill(jet.eta(), weight);
        h_em.h1_jet3_phi->Fill(jet.phi(), weight);
        h_em.h1_jet3_btag->Fill(jet.bDiscriminator(bTagName_), weight);
      }
      if ( jets_n >= 4 )
      {
        const auto& jet = out_jets->at(3);
        h_em.h1_jet4_m->Fill(jet.mass(), weight);
        h_em.h1_jet4_pt->Fill(jet.pt(), weight);
        h_em.h1_jet4_eta->Fill(jet.eta(), weight);
        h_em.h1_jet4_phi->Fill(jet.phi(), weight);
        h_em.h1_jet4_btag->Fill(jet.bDiscriminator(bTagName_), weight);
      }
      h_em.h1_bjets_n->Fill(bjets_n, weight);

      h_em.h1_event_st->Fill(leptons_st+jets_ht+met_pt, weight);

    }

    // Finalize cut step 1 and start cut step2
    if ( cutstep <= 1 ) break;
    ++icutstep; // =2

    if ( channel == CH_ELEL )
    {
      h_ee.hCutstep->Fill(icutstep, weight);
      h_ee.hCutstepNoweight->Fill(icutstep);
      h_ee.h2_vertex_n->Fill(nVertex, weight);
      h_ee.h2_met_pt->Fill(met_pt, weight);
      h_ee.h2_met_phi->Fill(met_phi, weight);
      h_ee.h2_leptons_n->Fill(leptons_n, weight);
      h_ee.h2_lepton1_pt->Fill(lepton1P4.pt(), weight);
      h_ee.h2_lepton1_eta->Fill(lepton1->eta(), weight);
      h_ee.h2_lepton1_phi->Fill(lepton1->phi(), weight);
      h_ee.h2_lepton1_q->Fill(lepton1->charge(), weight);
      h_ee.h2_lepton2_pt->Fill(lepton2P4.pt(), weight);
      h_ee.h2_lepton2_eta->Fill(lepton2->eta(), weight);
      h_ee.h2_lepton2_phi->Fill(lepton2->phi(), weight);
      h_ee.h2_lepton2_q->Fill(lepton2->charge(), weight);
      h_ee.h2_z_m->Fill(zP4.mass(), weight);
      h_ee.h2_z_pt->Fill(zP4.pt(), weight);
      h_ee.h2_z_eta->Fill(zP4.eta(), weight);
      h_ee.h2_z_phi->Fill(zP4.phi(), weight);
      h_ee.h2_jets_n->Fill(jets_n, weight);
      h_ee.h2_jets_ht->Fill(jets_ht, weight);
      for ( auto jet : *out_jets ) {
        h_ee.h2_jets_pt->Fill(jet.pt(), weight);
        h_ee.h2_jets_eta->Fill(jet.eta(), weight);
      }
      if ( jets_n >= 1 )
      {
        const auto& jet = out_jets->at(0);
        h_ee.h2_jet1_m->Fill(jet.mass(), weight);
        h_ee.h2_jet1_pt->Fill(jet.pt(), weight);
        h_ee.h2_jet1_eta->Fill(jet.eta(), weight);
        h_ee.h2_jet1_phi->Fill(jet.phi(), weight);
        h_ee.h2_jet1_btag->Fill(jet.bDiscriminator(bTagName_), weight);
      }
      if ( jets_n >= 2 )
      {
        const auto& jet = out_jets->at(1);
        h_ee.h2_jet2_m->Fill(jet.mass(), weight);
        h_ee.h2_jet2_pt->Fill(jet.pt(), weight);
        h_ee.h2_jet2_eta->Fill(jet.eta(), weight);
        h_ee.h2_jet2_phi->Fill(jet.phi(), weight);
        h_ee.h2_jet2_btag->Fill(jet.bDiscriminator(bTagName_), weight);
      }
      if ( jets_n >= 3 )
      {
        const auto& jet = out_jets->at(2);
        h_ee.h2_jet3_m->Fill(jet.mass(), weight);
        h_ee.h2_jet3_pt->Fill(jet.pt(), weight);
        h_ee.h2_jet3_eta->Fill(jet.eta(), weight);
        h_ee.h2_jet3_phi->Fill(jet.phi(), weight);
        h_ee.h2_jet3_btag->Fill(jet.bDiscriminator(bTagName_), weight);
      }
      if ( jets_n >= 4 )
      {
        const auto& jet = out_jets->at(3);
        h_ee.h2_jet4_m->Fill(jet.mass(), weight);
        h_ee.h2_jet4_pt->Fill(jet.pt(), weight);
        h_ee.h2_jet4_eta->Fill(jet.eta(), weight);
        h_ee.h2_jet4_phi->Fill(jet.phi(), weight);
        h_ee.h2_jet4_btag->Fill(jet.bDiscriminator(bTagName_), weight);
      }
      h_ee.h2_bjets_n->Fill(bjets_n, weight);

      h_ee.h2_event_st->Fill(leptons_st+jets_ht+met_pt, weight);
    }
    else if ( channel == CH_MUMU )
    {
      h_mm.hCutstep->Fill(icutstep, weight);
      h_mm.hCutstepNoweight->Fill(icutstep);
      h_mm.h2_vertex_n->Fill(nVertex, weight);
      h_mm.h2_met_pt->Fill(met_pt, weight);
      h_mm.h2_met_phi->Fill(met_phi, weight);
      h_mm.h2_leptons_n->Fill(leptons_n, weight);
      h_mm.h2_lepton1_pt->Fill(lepton1P4.pt(), weight);
      h_mm.h2_lepton1_eta->Fill(lepton1->eta(), weight);
      h_mm.h2_lepton1_phi->Fill(lepton1->phi(), weight);
      h_mm.h2_lepton1_q->Fill(lepton1->charge(), weight);
      h_mm.h2_lepton2_pt->Fill(lepton2P4.pt(), weight);
      h_mm.h2_lepton2_eta->Fill(lepton2->eta(), weight);
      h_mm.h2_lepton2_phi->Fill(lepton2->phi(), weight);
      h_mm.h2_lepton2_q->Fill(lepton2->charge(), weight);
      h_mm.h2_z_m->Fill(zP4.mass(), weight);
      h_mm.h2_z_pt->Fill(zP4.pt(), weight);
      h_mm.h2_z_eta->Fill(zP4.eta(), weight);
      h_mm.h2_z_phi->Fill(zP4.phi(), weight);
      h_mm.h2_jets_n->Fill(jets_n, weight);
      h_mm.h2_jets_ht->Fill(jets_ht, weight);
      for ( auto jet : *out_jets ) {
        h_mm.h2_jets_pt->Fill(jet.pt(), weight);
        h_mm.h2_jets_eta->Fill(jet.eta(), weight);
      }
      if ( jets_n >= 1 )
      {
        const auto& jet = out_jets->at(0);
        h_mm.h2_jet1_m->Fill(jet.mass(), weight);
        h_mm.h2_jet1_pt->Fill(jet.pt(), weight);
        h_mm.h2_jet1_eta->Fill(jet.eta(), weight);
        h_mm.h2_jet1_phi->Fill(jet.phi(), weight);
        h_mm.h2_jet1_btag->Fill(jet.bDiscriminator(bTagName_), weight);
      }
      if ( jets_n >= 2 )
      {
        const auto& jet = out_jets->at(1);
        h_mm.h2_jet2_m->Fill(jet.mass(), weight);
        h_mm.h2_jet2_pt->Fill(jet.pt(), weight);
        h_mm.h2_jet2_eta->Fill(jet.eta(), weight);
        h_mm.h2_jet2_phi->Fill(jet.phi(), weight);
        h_mm.h2_jet2_btag->Fill(jet.bDiscriminator(bTagName_), weight);
      }
      if ( jets_n >= 3 )
      {
        const auto& jet = out_jets->at(2);
        h_mm.h2_jet3_m->Fill(jet.mass(), weight);
        h_mm.h2_jet3_pt->Fill(jet.pt(), weight);
        h_mm.h2_jet3_eta->Fill(jet.eta(), weight);
        h_mm.h2_jet3_phi->Fill(jet.phi(), weight);
        h_mm.h2_jet3_btag->Fill(jet.bDiscriminator(bTagName_), weight);
      }
      if ( jets_n >= 4 )
      {
        const auto& jet = out_jets->at(3);
        h_mm.h2_jet4_m->Fill(jet.mass(), weight);
        h_mm.h2_jet4_pt->Fill(jet.pt(), weight);
        h_mm.h2_jet4_eta->Fill(jet.eta(), weight);
        h_mm.h2_jet4_phi->Fill(jet.phi(), weight);
        h_mm.h2_jet4_btag->Fill(jet.bDiscriminator(bTagName_), weight);
      }
      h_mm.h2_bjets_n->Fill(bjets_n, weight);

      h_mm.h2_event_st->Fill(leptons_st+jets_ht+met_pt, weight);
    }
    else if ( channel == CH_MUEL )
    {
      h_em.hCutstep->Fill(icutstep, weight);
      h_em.hCutstepNoweight->Fill(icutstep);
      h_em.h2_vertex_n->Fill(nVertex, weight);
      h_em.h2_met_pt->Fill(met_pt, weight);
      h_em.h2_met_phi->Fill(met_phi, weight);
      h_em.h2_leptons_n->Fill(leptons_n, weight);
      h_em.h2_lepton1_pt->Fill(lepton1P4.pt(), weight);
      h_em.h2_lepton1_eta->Fill(lepton1->eta(), weight);
      h_em.h2_lepton1_phi->Fill(lepton1->phi(), weight);
      h_em.h2_lepton1_q->Fill(lepton1->charge(), weight);
      h_em.h2_lepton2_pt->Fill(lepton2P4.pt(), weight);
      h_em.h2_lepton2_eta->Fill(lepton2->eta(), weight);
      h_em.h2_lepton2_phi->Fill(lepton2->phi(), weight);
      h_em.h2_lepton2_q->Fill(lepton2->charge(), weight);
      h_em.h2_z_m->Fill(zP4.mass(), weight);
      h_em.h2_z_pt->Fill(zP4.pt(), weight);
      h_em.h2_z_eta->Fill(zP4.eta(), weight);
      h_em.h2_z_phi->Fill(zP4.phi(), weight);
      h_em.h2_jets_n->Fill(jets_n, weight);
      h_em.h2_jets_ht->Fill(jets_ht, weight);
      for ( auto jet : *out_jets ) {
        h_em.h2_jets_pt->Fill(jet.pt(), weight);
        h_em.h2_jets_eta->Fill(jet.eta(), weight);
      }
      if ( jets_n >= 1 )
      {
        const auto& jet = out_jets->at(0);
        h_em.h2_jet1_m->Fill(jet.mass(), weight);
        h_em.h2_jet1_pt->Fill(jet.pt(), weight);
        h_em.h2_jet1_eta->Fill(jet.eta(), weight);
        h_em.h2_jet1_phi->Fill(jet.phi(), weight);
        h_em.h2_jet1_btag->Fill(jet.bDiscriminator(bTagName_), weight);
      }
      if ( jets_n >= 2 )
      {
        const auto& jet = out_jets->at(1);
        h_em.h2_jet2_m->Fill(jet.mass(), weight);
        h_em.h2_jet2_pt->Fill(jet.pt(), weight);
        h_em.h2_jet2_eta->Fill(jet.eta(), weight);
        h_em.h2_jet2_phi->Fill(jet.phi(), weight);
        h_em.h2_jet2_btag->Fill(jet.bDiscriminator(bTagName_), weight);
      }
      if ( jets_n >= 3 )
      {
        const auto& jet = out_jets->at(2);
        h_em.h2_jet3_m->Fill(jet.mass(), weight);
        h_em.h2_jet3_pt->Fill(jet.pt(), weight);
        h_em.h2_jet3_eta->Fill(jet.eta(), weight);
        h_em.h2_jet3_phi->Fill(jet.phi(), weight);
        h_em.h2_jet3_btag->Fill(jet.bDiscriminator(bTagName_), weight);
      }
      if ( jets_n >= 4 )
      {
        const auto& jet = out_jets->at(3);
        h_em.h2_jet4_m->Fill(jet.mass(), weight);
        h_em.h2_jet4_pt->Fill(jet.pt(), weight);
        h_em.h2_jet4_eta->Fill(jet.eta(), weight);
        h_em.h2_jet4_phi->Fill(jet.phi(), weight);
        h_em.h2_jet4_btag->Fill(jet.bDiscriminator(bTagName_), weight);
      }
      h_em.h2_bjets_n->Fill(bjets_n, weight);

      h_em.h2_event_st->Fill(leptons_st+jets_ht+met_pt, weight);

    }

    // Finalize cut step 2 and start cut step3
    if ( cutstep <= 2 ) break;
    ++icutstep; // =3

    const auto& jet1 = out_jets->at(0);
    const auto& jet2 = out_jets->at(1);

    if ( channel == CH_ELEL )
    {
      h_ee.hCutstep->Fill(icutstep, weight);
      h_ee.hCutstepNoweight->Fill(icutstep);
      h_ee.h3_vertex_n->Fill(nVertex, weight);
      h_ee.h3_met_pt->Fill(met_pt, weight);
      h_ee.h3_met_phi->Fill(met_phi, weight);
      h_ee.h3_z_m->Fill(zP4.mass(), weight);
      h_ee.h3_z_pt->Fill(zP4.pt(), weight);
      h_ee.h3_z_eta->Fill(zP4.eta(), weight);
      h_ee.h3_z_phi->Fill(zP4.phi(), weight);
      h_ee.h3_jets_n->Fill(jets_n, weight);
      h_ee.h3_jets_ht->Fill(jets_ht, weight);
      for ( auto jet : *out_jets ) {
        h_ee.h3_jets_pt->Fill(jet.pt(), weight);
        h_ee.h3_jets_eta->Fill(jet.eta(), weight);
      }
      h_ee.h3_jet1_m->Fill(jet1.mass(), weight);
      h_ee.h3_jet1_pt->Fill(jet1.pt(), weight);
      h_ee.h3_jet1_eta->Fill(jet1.eta(), weight);
      h_ee.h3_jet1_phi->Fill(jet1.phi(), weight);
      h_ee.h3_jet1_btag->Fill(jet1.bDiscriminator(bTagName_), weight);
      h_ee.h3_jet2_m->Fill(jet2.mass(), weight);
      h_ee.h3_jet2_pt->Fill(jet2.pt(), weight);
      h_ee.h3_jet2_eta->Fill(jet2.eta(), weight);
      h_ee.h3_jet2_phi->Fill(jet2.phi(), weight);
      h_ee.h3_jet2_btag->Fill(jet2.bDiscriminator(bTagName_), weight);
      if ( jets_n >= 3 )
      {
        const auto& jet = out_jets->at(2);
        h_ee.h3_jet3_m->Fill(jet.mass(), weight);
        h_ee.h3_jet3_pt->Fill(jet.pt(), weight);
        h_ee.h3_jet3_eta->Fill(jet.eta(), weight);
        h_ee.h3_jet3_phi->Fill(jet.phi(), weight);
        h_ee.h3_jet3_btag->Fill(jet.bDiscriminator(bTagName_), weight);
      }
      if ( jets_n >= 4 )
      {
        const auto& jet = out_jets->at(3);
        h_ee.h3_jet4_m->Fill(jet.mass(), weight);
        h_ee.h3_jet4_pt->Fill(jet.pt(), weight);
        h_ee.h3_jet4_eta->Fill(jet.eta(), weight);
        h_ee.h3_jet4_phi->Fill(jet.phi(), weight);
        h_ee.h3_jet4_btag->Fill(jet.bDiscriminator(bTagName_), weight);
      }
      h_ee.h3_bjets_n->Fill(bjets_n, weight);

      h_ee.h3_event_st->Fill(leptons_st+jets_ht+met_pt, weight);
    }
    else if ( channel == CH_MUMU )
    {
      h_mm.hCutstep->Fill(icutstep, weight);
      h_mm.hCutstepNoweight->Fill(icutstep);
      h_mm.h3_vertex_n->Fill(nVertex, weight);
      h_mm.h3_met_pt->Fill(met_pt, weight);
      h_mm.h3_met_phi->Fill(met_phi, weight);
      h_mm.h3_z_m->Fill(zP4.mass(), weight);
      h_mm.h3_z_pt->Fill(zP4.pt(), weight);
      h_mm.h3_z_eta->Fill(zP4.eta(), weight);
      h_mm.h3_z_phi->Fill(zP4.phi(), weight);
      h_mm.h3_jets_n->Fill(jets_n, weight);
      h_mm.h3_jets_ht->Fill(jets_ht, weight);
      for ( auto jet : *out_jets ) {
        h_mm.h3_jets_pt->Fill(jet.pt(), weight);
        h_mm.h3_jets_eta->Fill(jet.eta(), weight);
      }
      h_mm.h3_jet1_m->Fill(jet1.mass(), weight);
      h_mm.h3_jet1_pt->Fill(jet1.pt(), weight);
      h_mm.h3_jet1_eta->Fill(jet1.eta(), weight);
      h_mm.h3_jet1_phi->Fill(jet1.phi(), weight);
      h_mm.h3_jet1_btag->Fill(jet1.bDiscriminator(bTagName_), weight);
      h_mm.h3_jet2_m->Fill(jet2.mass(), weight);
      h_mm.h3_jet2_pt->Fill(jet2.pt(), weight);
      h_mm.h3_jet2_eta->Fill(jet2.eta(), weight);
      h_mm.h3_jet2_phi->Fill(jet2.phi(), weight);
      h_mm.h3_jet2_btag->Fill(jet2.bDiscriminator(bTagName_), weight);
      if ( jets_n >= 3 )
      {
        const auto& jet = out_jets->at(2);
        h_mm.h3_jet3_m->Fill(jet.mass(), weight);
        h_mm.h3_jet3_pt->Fill(jet.pt(), weight);
        h_mm.h3_jet3_eta->Fill(jet.eta(), weight);
        h_mm.h3_jet3_phi->Fill(jet.phi(), weight);
        h_mm.h3_jet3_btag->Fill(jet.bDiscriminator(bTagName_), weight);
      }
      if ( jets_n >= 4 )
      {
        const auto& jet = out_jets->at(3);
        h_mm.h3_jet4_m->Fill(jet.mass(), weight);
        h_mm.h3_jet4_pt->Fill(jet.pt(), weight);
        h_mm.h3_jet4_eta->Fill(jet.eta(), weight);
        h_mm.h3_jet4_phi->Fill(jet.phi(), weight);
        h_mm.h3_jet4_btag->Fill(jet.bDiscriminator(bTagName_), weight);
      }
      h_mm.h3_bjets_n->Fill(bjets_n, weight);

      h_mm.h3_event_st->Fill(leptons_st+jets_ht+met_pt, weight);

    }
    else if ( channel == CH_MUEL )
    {
      h_em.hCutstep->Fill(icutstep, weight);
      h_em.hCutstepNoweight->Fill(icutstep);
      h_em.h3_vertex_n->Fill(nVertex, weight);
      h_em.h3_met_pt->Fill(met_pt, weight);
      h_em.h3_met_phi->Fill(met_phi, weight);
      h_em.h3_z_m->Fill(zP4.mass(), weight);
      h_em.h3_z_pt->Fill(zP4.pt(), weight);
      h_em.h3_z_eta->Fill(zP4.eta(), weight);
      h_em.h3_z_phi->Fill(zP4.phi(), weight);
      h_em.h3_jets_n->Fill(jets_n, weight);
      h_em.h3_jets_ht->Fill(jets_ht, weight);
      for ( auto jet : *out_jets ) {
        h_em.h3_jets_pt->Fill(jet.pt(), weight);
        h_em.h3_jets_eta->Fill(jet.eta(), weight);
      }
      h_em.h3_jet1_m->Fill(jet1.mass(), weight);
      h_em.h3_jet1_pt->Fill(jet1.pt(), weight);
      h_em.h3_jet1_eta->Fill(jet1.eta(), weight);
      h_em.h3_jet1_phi->Fill(jet1.phi(), weight);
      h_em.h3_jet1_btag->Fill(jet1.bDiscriminator(bTagName_), weight);
      h_em.h3_jet2_m->Fill(jet2.mass(), weight);
      h_em.h3_jet2_pt->Fill(jet2.pt(), weight);
      h_em.h3_jet2_eta->Fill(jet2.eta(), weight);
      h_em.h3_jet2_phi->Fill(jet2.phi(), weight);
      h_em.h3_jet2_btag->Fill(jet2.bDiscriminator(bTagName_), weight);
      if ( jets_n >= 3 )
      {
        const auto& jet = out_jets->at(2);
        h_em.h3_jet3_m->Fill(jet.mass(), weight);
        h_em.h3_jet3_pt->Fill(jet.pt(), weight);
        h_em.h3_jet3_eta->Fill(jet.eta(), weight);
        h_em.h3_jet3_phi->Fill(jet.phi(), weight);
        h_em.h3_jet3_btag->Fill(jet.bDiscriminator(bTagName_), weight);
      }
      if ( jets_n >= 4 )
      {
        const auto& jet = out_jets->at(3);
        h_em.h3_jet4_m->Fill(jet.mass(), weight);
        h_em.h3_jet4_pt->Fill(jet.pt(), weight);
        h_em.h3_jet4_eta->Fill(jet.eta(), weight);
        h_em.h3_jet4_phi->Fill(jet.phi(), weight);
        h_em.h3_jet4_btag->Fill(jet.bDiscriminator(bTagName_), weight);
      }
      h_em.h3_bjets_n->Fill(bjets_n, weight);

      h_em.h3_event_st->Fill(leptons_st+jets_ht+met_pt, weight);

    }

    // Finalize cut step 3 and start cut step 4
    if ( cutstep <= 3 ) break;
    ++icutstep; // =4

    if ( channel == CH_ELEL )
    {
      h_ee.hCutstep->Fill(icutstep, weight);
      h_ee.hCutstepNoweight->Fill(icutstep);
      h_ee.h4_vertex_n->Fill(nVertex, weight);
      h_ee.h4_met_pt->Fill(met_pt, weight);
      h_ee.h4_met_phi->Fill(met_phi, weight);
      h_ee.h4_z_m->Fill(zP4.mass(), weight);
      h_ee.h4_z_pt->Fill(zP4.pt(), weight);
      h_ee.h4_z_eta->Fill(zP4.eta(), weight);
      h_ee.h4_z_phi->Fill(zP4.phi(), weight);
      h_ee.h4_jets_n->Fill(jets_n, weight);
      h_ee.h4_jets_ht->Fill(jets_ht, weight);
      for ( auto jet : *out_jets ) {
        h_ee.h4_jets_pt->Fill(jet.pt(), weight);
        h_ee.h4_jets_eta->Fill(jet.eta(), weight);
      }
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
        const auto& jet = out_jets->at(2);
        h_ee.h4_jet3_m->Fill(jet.mass(), weight);
        h_ee.h4_jet3_pt->Fill(jet.pt(), weight);
        h_ee.h4_jet3_eta->Fill(jet.eta(), weight);
        h_ee.h4_jet3_phi->Fill(jet.phi(), weight);
        h_ee.h4_jet3_btag->Fill(jet.bDiscriminator(bTagName_), weight);
      }
      if ( jets_n >= 4 )
      {
        const auto& jet = out_jets->at(3);
        h_ee.h4_jet4_m->Fill(jet.mass(), weight);
        h_ee.h4_jet4_pt->Fill(jet.pt(), weight);
        h_ee.h4_jet4_eta->Fill(jet.eta(), weight);
        h_ee.h4_jet4_phi->Fill(jet.phi(), weight);
        h_ee.h4_jet4_btag->Fill(jet.bDiscriminator(bTagName_), weight);
      }
      h_ee.h4_bjets_n->Fill(bjets_n, weight);

      h_ee.h4_event_st->Fill(leptons_st+jets_ht+met_pt, weight);
    }
    else if ( channel == CH_MUMU )
    {
      h_mm.hCutstep->Fill(icutstep, weight);
      h_mm.hCutstepNoweight->Fill(icutstep);
      h_mm.h4_vertex_n->Fill(nVertex, weight);
      h_mm.h4_met_pt->Fill(met_pt, weight);
      h_mm.h4_met_phi->Fill(met_phi, weight);
      h_mm.h4_z_m->Fill(zP4.mass(), weight);
      h_mm.h4_z_pt->Fill(zP4.pt(), weight);
      h_mm.h4_z_eta->Fill(zP4.eta(), weight);
      h_mm.h4_z_phi->Fill(zP4.phi(), weight);
      h_mm.h4_jets_n->Fill(jets_n, weight);
      h_mm.h4_jets_ht->Fill(jets_ht, weight);
      for ( auto jet : *out_jets ) {
        h_mm.h4_jets_pt->Fill(jet.pt(), weight);
        h_mm.h4_jets_eta->Fill(jet.eta(), weight);
      }
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
        const auto& jet = out_jets->at(2);
        h_mm.h4_jet3_m->Fill(jet.mass(), weight);
        h_mm.h4_jet3_pt->Fill(jet.pt(), weight);
        h_mm.h4_jet3_eta->Fill(jet.eta(), weight);
        h_mm.h4_jet3_phi->Fill(jet.phi(), weight);
        h_mm.h4_jet3_btag->Fill(jet.bDiscriminator(bTagName_), weight);
      }
      if ( jets_n >= 4 )
      {
        const auto& jet = out_jets->at(3);
        h_mm.h4_jet4_m->Fill(jet.mass(), weight);
        h_mm.h4_jet4_pt->Fill(jet.pt(), weight);
        h_mm.h4_jet4_eta->Fill(jet.eta(), weight);
        h_mm.h4_jet4_phi->Fill(jet.phi(), weight);
        h_mm.h4_jet4_btag->Fill(jet.bDiscriminator(bTagName_), weight);
      }
      h_mm.h4_bjets_n->Fill(bjets_n, weight);

      h_mm.h4_event_st->Fill(leptons_st+jets_ht+met_pt, weight);

    }
    else if ( channel == CH_MUEL )
    {
      h_em.hCutstep->Fill(icutstep, weight);
      h_em.hCutstepNoweight->Fill(icutstep);
      h_em.h4_vertex_n->Fill(nVertex, weight);
      h_em.h4_met_pt->Fill(met_pt, weight);
      h_em.h4_met_phi->Fill(met_phi, weight);
      h_em.h4_z_m->Fill(zP4.mass(), weight);
      h_em.h4_z_pt->Fill(zP4.pt(), weight);
      h_em.h4_z_eta->Fill(zP4.eta(), weight);
      h_em.h4_z_phi->Fill(zP4.phi(), weight);
      h_em.h4_jets_n->Fill(jets_n, weight);
      h_em.h4_jets_ht->Fill(jets_ht, weight);
      for ( auto jet : *out_jets ) {
        h_em.h4_jets_pt->Fill(jet.pt(), weight);
        h_em.h4_jets_eta->Fill(jet.eta(), weight);
      }
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
        const auto& jet = out_jets->at(2);
        h_em.h4_jet3_m->Fill(jet.mass(), weight);
        h_em.h4_jet3_pt->Fill(jet.pt(), weight);
        h_em.h4_jet3_eta->Fill(jet.eta(), weight);
        h_em.h4_jet3_phi->Fill(jet.phi(), weight);
        h_em.h4_jet3_btag->Fill(jet.bDiscriminator(bTagName_), weight);
      }
      if ( jets_n >= 4 )
      {
        const auto& jet = out_jets->at(3);
        h_em.h4_jet4_m->Fill(jet.mass(), weight);
        h_em.h4_jet4_pt->Fill(jet.pt(), weight);
        h_em.h4_jet4_eta->Fill(jet.eta(), weight);
        h_em.h4_jet4_phi->Fill(jet.phi(), weight);
        h_em.h4_jet4_btag->Fill(jet.bDiscriminator(bTagName_), weight);
      }
      h_em.h4_bjets_n->Fill(bjets_n, weight);

      h_em.h4_event_st->Fill(leptons_st+jets_ht+met_pt, weight);

    }

    // Finalize cut step 4 and start cut step 5
    if ( cutstep <= 4 ) break;
    ++icutstep; // =5

    if ( channel == CH_ELEL )
    {
      h_ee.hCutstep->Fill(icutstep, weight);
      h_ee.hCutstepNoweight->Fill(icutstep);
      h_ee.h5_vertex_n->Fill(nVertex, weight);
      h_ee.h5_met_pt->Fill(met_pt, weight);
      h_ee.h5_met_phi->Fill(met_phi, weight);
      h_ee.h5_z_m->Fill(zP4.mass(), weight);
      h_ee.h5_z_pt->Fill(zP4.pt(), weight);
      h_ee.h5_z_eta->Fill(zP4.eta(), weight);
      h_ee.h5_z_phi->Fill(zP4.phi(), weight);
      h_ee.h5_jets_n->Fill(jets_n, weight);
      h_ee.h5_jets_ht->Fill(jets_ht, weight);
      for ( auto jet : *out_jets ) {
        h_ee.h5_jets_pt->Fill(jet.pt(), weight);
        h_ee.h5_jets_eta->Fill(jet.eta(), weight);
      }
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
        const auto& jet = out_jets->at(2);
        h_ee.h5_jet3_m->Fill(jet.mass(), weight);
        h_ee.h5_jet3_pt->Fill(shiftedJetPt(jet), weight);
        h_ee.h5_jet3_eta->Fill(jet.eta(), weight);
        h_ee.h5_jet3_phi->Fill(jet.phi(), weight);
        h_ee.h5_jet3_btag->Fill(jet.bDiscriminator(bTagName_), weight);
      }
      if ( jets_n >= 4 )
      {
        const auto& jet = out_jets->at(3);
        h_ee.h5_jet4_m->Fill(jet.mass(), weight);
        h_ee.h5_jet4_pt->Fill(shiftedJetPt(jet), weight);
        h_ee.h5_jet4_eta->Fill(jet.eta(), weight);
        h_ee.h5_jet4_phi->Fill(jet.phi(), weight);
        h_ee.h5_jet4_btag->Fill(jet.bDiscriminator(bTagName_), weight);
      }
      h_ee.h5_bjets_n->Fill(bjets_n, weight);

      h_ee.h5_event_st->Fill(leptons_st+jets_ht+met_pt, weight);
    }
    else if ( channel == CH_MUMU )
    {
      h_mm.hCutstep->Fill(icutstep, weight);
      h_mm.hCutstepNoweight->Fill(icutstep);
      h_mm.h5_vertex_n->Fill(nVertex, weight);
      h_mm.h5_met_pt->Fill(met_pt, weight);
      h_mm.h5_met_phi->Fill(met_phi, weight);
      h_mm.h5_z_m->Fill(zP4.mass(), weight);
      h_mm.h5_z_pt->Fill(zP4.pt(), weight);
      h_mm.h5_z_eta->Fill(zP4.eta(), weight);
      h_mm.h5_z_phi->Fill(zP4.phi(), weight);
      h_mm.h5_jets_n->Fill(jets_n, weight);
      h_mm.h5_jets_ht->Fill(jets_ht, weight);
      for ( auto jet : *out_jets ) {
        h_mm.h5_jets_pt->Fill(jet.pt(), weight);
        h_mm.h5_jets_eta->Fill(jet.eta(), weight);
      }
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
        const auto& jet = out_jets->at(2);
        h_mm.h5_jet3_m->Fill(jet.mass(), weight);
        h_mm.h5_jet3_pt->Fill(shiftedJetPt(jet), weight);
        h_mm.h5_jet3_eta->Fill(jet.eta(), weight);
        h_mm.h5_jet3_phi->Fill(jet.phi(), weight);
        h_mm.h5_jet3_btag->Fill(jet.bDiscriminator(bTagName_), weight);
      }
      if ( jets_n >= 4 )
      {
        const auto& jet = out_jets->at(3);
        h_mm.h5_jet4_m->Fill(jet.mass(), weight);
        h_mm.h5_jet4_pt->Fill(shiftedJetPt(jet), weight);
        h_mm.h5_jet4_eta->Fill(jet.eta(), weight);
        h_mm.h5_jet4_phi->Fill(jet.phi(), weight);
        h_mm.h5_jet4_btag->Fill(jet.bDiscriminator(bTagName_), weight);
      }
      h_mm.h5_bjets_n->Fill(bjets_n, weight);

      h_mm.h5_event_st->Fill(leptons_st+jets_ht+met_pt, weight);

    }
    else if ( channel == CH_MUEL )
    {
      h_em.hCutstep->Fill(icutstep, weight);
      h_em.hCutstepNoweight->Fill(icutstep);
      h_em.h5_vertex_n->Fill(nVertex, weight);
      h_em.h5_met_pt->Fill(met_pt, weight);
      h_em.h5_met_phi->Fill(met_phi, weight);
      h_em.h5_z_m->Fill(zP4.mass(), weight);
      h_em.h5_z_pt->Fill(zP4.pt(), weight);
      h_em.h5_z_eta->Fill(zP4.eta(), weight);
      h_em.h5_z_phi->Fill(zP4.phi(), weight);
      h_em.h5_jets_n->Fill(jets_n, weight);
      h_em.h5_jets_ht->Fill(jets_ht, weight);
      for ( auto jet : *out_jets ) {
        h_em.h5_jets_pt->Fill(jet.pt(), weight);
        h_em.h5_jets_eta->Fill(jet.eta(), weight);
      }
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
        const auto& jet = out_jets->at(2);
        h_em.h5_jet3_m->Fill(jet.mass(), weight);
        h_em.h5_jet3_pt->Fill(shiftedJetPt(jet), weight);
        h_em.h5_jet3_eta->Fill(jet.eta(), weight);
        h_em.h5_jet3_phi->Fill(jet.phi(), weight);
        h_em.h5_jet3_btag->Fill(jet.bDiscriminator(bTagName_), weight);
      }
      if ( jets_n >= 4 )
      {
        const auto& jet = out_jets->at(3);
        h_em.h5_jet4_m->Fill(jet.mass(), weight);
        h_em.h5_jet4_pt->Fill(shiftedJetPt(jet), weight);
        h_em.h5_jet4_eta->Fill(jet.eta(), weight);
        h_em.h5_jet4_phi->Fill(jet.phi(), weight);
        h_em.h5_jet4_btag->Fill(jet.bDiscriminator(bTagName_), weight);
      }
      h_em.h5_bjets_n->Fill(bjets_n, weight);

      h_em.h5_event_st->Fill(leptons_st+jets_ht+met_pt, weight);
    }
  } // switch(1)

  // Cutsomized cutflow without z-veto cut to be used in DY estimation and other studies
  for ( int istep=1, nstep=cutstepBits.size(); istep<=nstep; ++istep )
  {
    if ( istep != 2 and !cutstepBits[istep-1] ) break;

    if ( channel == CH_ELEL )
    {
      switch ( istep )
      {
        case 5: h_ee.h5_z_m_noveto->Fill(z_m, weight);
        case 4: h_ee.h4_z_m_noveto->Fill(z_m, weight);
        case 3: h_ee.h3_z_m_noveto->Fill(z_m, weight);
        case 2: h_ee.h2_z_m_noveto->Fill(z_m, weight);
        case 1: h_ee.h1_z_m_noveto->Fill(z_m, weight);
      }
    }
    else if ( channel == CH_MUMU )
    {
      switch ( istep )
      {
        case 5: h_mm.h5_z_m_noveto->Fill(z_m, weight);
        case 4: h_mm.h4_z_m_noveto->Fill(z_m, weight);
        case 3: h_mm.h3_z_m_noveto->Fill(z_m, weight);
        case 2: h_mm.h2_z_m_noveto->Fill(z_m, weight);
        case 1: h_mm.h1_z_m_noveto->Fill(z_m, weight);
      }
    }
    else if ( channel == CH_MUEL )
    {
      switch ( istep )
      {
        case 5: h_em.h5_z_m_noveto->Fill(z_m, weight);
        case 4: h_em.h4_z_m_noveto->Fill(z_m, weight);
        case 3: h_em.h3_z_m_noveto->Fill(z_m, weight);
        case 2: h_em.h2_z_m_noveto->Fill(z_m, weight);
        case 1: h_em.h1_z_m_noveto->Fill(z_m, weight);
      }
    }
  }

  // Fill cut flow 2D plot
  for ( int istep=1, nstep=cutstepBits.size(); istep<=nstep; ++istep )
  {
    const bool res1 = cutstepBits[istep-1];

    // Fill diagonal terms
    if      ( channel == CH_ELEL ) {
      h_ee.h2Cutstep->Fill(istep, istep, res1*weight);
      h_ee.h2CutstepNoweight->Fill(istep, istep, res1);
    }
    else if ( channel == CH_MUMU ) {
      h_mm.h2Cutstep->Fill(istep, istep, res1*weight);
      h_mm.h2CutstepNoweight->Fill(istep, istep, res1);
    }
    else if ( channel == CH_MUEL ) {
      h_em.h2Cutstep->Fill(istep, istep, res1*weight);
      h_em.h2CutstepNoweight->Fill(istep, istep, res1);
    }

    // Fill correlations and anti-correlations
    for ( int jstep=1; jstep<istep; ++jstep )
    {
      const bool res2 = cutstepBits[jstep-1];
      const int result = res1 && res2;
      const int aresult = res1 && !res2;
      if      ( channel == CH_ELEL ) {
        h_ee.h2Cutstep->Fill(istep, jstep, result*weight);
        h_ee.h2CutstepNoweight->Fill(istep, jstep, result);
        h_ee.h2Cutstep->Fill(jstep, istep, aresult*weight);
        h_ee.h2CutstepNoweight->Fill(jstep, istep, aresult);
      }
      else if ( channel == CH_MUMU ) {
        h_mm.h2Cutstep->Fill(istep, jstep, result*weight);
        h_mm.h2CutstepNoweight->Fill(istep, jstep, result);
        h_mm.h2Cutstep->Fill(jstep, istep, aresult*weight);
        h_mm.h2CutstepNoweight->Fill(jstep, istep, aresult);
      }
      else if ( channel == CH_MUEL ) {
        h_em.h2Cutstep->Fill(istep, jstep, result*weight);
        h_em.h2CutstepNoweight->Fill(istep, jstep, result);
        h_em.h2Cutstep->Fill(jstep, istep, aresult*weight);
        h_em.h2CutstepNoweight->Fill(jstep, istep, aresult);
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

TTLLEventSelector::~TTLLEventSelector()
{
  if ( h_em.hCutstepNoweight )
  {
    cout << "---- cut flows without weight ----\n";
    cout << "Step\tee\tmumu\temu\n";
    const int n = h_em.hCutstepNoweight->GetNbinsX();
    for ( int i=1; i<=n; ++i )
    {
      const string name(h_ee.hCutstepNoweight->GetXaxis()->GetBinLabel(i));
      if ( name.empty() ) break;
      cout << name.substr(0, name.find(' '));
      cout << '\t' << h_ee.hCutstepNoweight->GetBinContent(i);
      cout << '\t' << h_mm.hCutstepNoweight->GetBinContent(i);
      cout << '\t' << h_em.hCutstepNoweight->GetBinContent(i) << endl;
    }
    cout << "-----------------------------------\n";
  }
}

#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(TTLLEventSelector);

