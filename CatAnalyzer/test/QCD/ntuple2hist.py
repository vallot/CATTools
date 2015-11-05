#!/usr/bin/env python
import ROOT, os, sys, copy
from array import array

ROOT.gROOT.SetBatch()
pi = ROOT.TMath.Pi()

in_type = sys.argv[1]
in_f = sys.argv[2]
#out_f = sys.argv[3]
out_f = in_f.replace("_ntuple.root", "_hist.root")

      

in_rf = ROOT.TFile(in_f)
mc = True
if not in_type == 'mc':
  mc = False


def hist_maker(name, title, bin_set, x_name, y_name, tr, br, w):
  if bin_set[2] == 2500:
    bin = []
    for x in xrange(15):
      bin.append(2500.0/2.0*float(x)/15.0)
    bin.append(1500)
    bin.append(2000)
    bin.append(2500)
    pt_bin = array('d', bin)
    hist = ROOT.TH1F(name, title, len(pt_bin)-1, pt_bin)
  else:
    hist = ROOT.TH1F(name, title, bin_set[0], bin_set[1], bin_set[2])
  hist.GetXaxis().SetTitle(x_name)
  hist.GetYaxis().SetTitle(y_name)
  hist.Sumw2()
  tr.Project(name, br, w)
  return hist

def hist2_maker(name,tr,br_reco,br_gen, w):
  hist2 = ROOT.TH2F(name, name, 18,0,pi, 18,0,pi)
  hist2.GetXaxis().SetTitle("reco_beta")
  hist2.GetYaxis().SetTitle("gen_beta")
  hist2.Sumw2()
  tr.Project(name, br_gen+":"+br_reco, w)
  return hist2
### analysis cuts for data

if mc:
    sys_e = ["nom", "jer_u", "jer_d", "jar", "pu_u", "pu_d"]
else:
    sys_e = ["nom", "jes_u", "jes_d"]    
hist_l = []
for sys in sys_e:
  if mc:
    e_w = "(1.0*pileupWeight)"
  else:
    e_w = "(1.0)"
  ## event selection
  r_mass_cut = "(raw_mass > 220)*"
  #r_mass_cut = "(%s_raw_mass > 250)*"%sys
  dr_cut = "(del_r>0.5)*(del_r<1.5)*"
  dr12_cut = "(abs(abs(del_r12)-%f)<1.0)*"%pi
  eta_cut = "(abs(jet1_eta)<2.5)*"
  min_pt_cut = "(jet1_pt>30)*(jet2_pt>30)*(jet3_pt>30)*"
  met_cut = "(metSig<0.3)*"
  #d_cut = r_mass_cut+dr_cut+dr12_cut+eta_cut+min_pt_cut+met_cut
  d_cut = r_mass_cut+dr_cut+dr12_cut+eta_cut+min_pt_cut+met_cut
  ## event classification

  l_eta = "(abs(jet2_eta)>0.0)*(abs(jet2_eta)<0.8)*"
  #m_eta = "(abs(%s_jet2_eta)>0.8)*(abs(%s_jet2_eta)<1.5)*"%(sys,sys)
  h_eta = "(abs(jet2_eta)>0.8)*(abs(jet2_eta)<2.5)*"

  h_pt = "(jet1_pt>740)*(jet1_pt<2500)*(hlt_450_pass == 1)"
  m_pt = "(jet1_pt>640)*(jet1_pt<2500)*(hlt_400_pass == 1)"
  l_pt = "(jet1_pt>510)*(jet1_pt<2500)*(hlt_320_pass == 1)"

  ## cc hist
  eta_bin = ["low_eta", "high_eta"]
  #eta_bin = ["low_eta", "medium_eta", "high_eta"]
  #eta_bin_cut = [l_eta, m_eta, h_eta]
  eta_bin_cut = [l_eta, h_eta]
  pt_bin = ["low_pt", "medium_pt", "high_pt"]
  pt_bin_cut = [l_pt, m_pt, h_pt]
  beta_l = ["beta", "del_eta", "del_phi", "del_r", "raw_mass","del_r12"] 
  beta_bin = [[18, 0, pi], [30, -3, 3], [30, -3, 3], [30, 0, 3], [50, 0, 5000], [30,0,3]]
  jet_l = ["pt", "eta", "phi"]
  jet_bin = [[[30,0,2500],[30,-3,3],[30,-pi,pi]],[[30,0,2500],[30,-3,3],[30,-pi,pi]],[[30,0,2500],[30,-3,3],[30,-pi,pi]]]
  ev_l = ["njet", "metSig", "nvtx"]  
  ev_bin = [[30, 0, 30],[100, 0, 1], [50, 0, 50]]
    


  for eta_i, eta_loop in enumerate(eta_bin):
    for pt_i, pt_loop in enumerate(pt_bin):
      cut = d_cut + eta_bin_cut[eta_i] + pt_bin_cut[pt_i]
      if sys.startswith("pu"):
        tr = in_rf.Get("cc/nom").Clone(eta_loop+pt_loop)
        if sys == "pu_u":
          e_w = "(1.0*pileupWeight_up)"
        else:
          e_w = "(1.0*pileupWeight_dn)"
      else:
        tr = in_rf.Get("cc/%s"%sys).Clone(eta_loop+pt_loop)
      tr.Draw(">>eventList", cut)
      el = ROOT.gDirectory.Get("eventList")
      tr.SetEventList(el)
 
      for i, beta_loop in enumerate(beta_l):
        name = eta_loop+"_"+pt_loop+"_"+sys+"_"+beta_loop
        title = name
        bin_set = beta_bin[i]
        x_name = beta_loop
        y_name = "count"
        br = beta_loop
        hist_l.append(copy.deepcopy(hist_maker(name, title, bin_set, x_name, y_name, tr, br, e_w)))
      for ji in xrange(3):
        for i, jet_loop in enumerate(jet_l):
          name = eta_loop+"_"+pt_loop+"_"+sys+"_jet%d_"%(ji+1)+jet_loop
          title = name
          bin_set = jet_bin[pt_i][i]
          x_name = jet_loop
          y_name = "count"
          br = "jet%d_"%(ji+1)+jet_loop
          hist_l.append(copy.deepcopy(hist_maker(name, title, bin_set, x_name, y_name, tr, br, e_w)))
      for i, ev_loop in enumerate(ev_l):
        name = eta_loop+"_"+pt_loop+"_"+sys+"_"+ev_loop
        title = name
        bin_set = ev_bin[i]
        x_name = ev_loop
        y_name = "count"
        br = ev_loop
        hist_l.append(copy.deepcopy(hist_maker(name, title, bin_set, x_name, y_name, tr, br, e_w)))

      if (sys == "nom" and mc):
        name = eta_loop+"_"+pt_loop+"_gen_beta"
        bin_set = beta_bin[0]
        x_name = "#beta GEN"
        y_name = "count"
        br = "gbeta"
        hist_l.append(copy.deepcopy(hist_maker(name, title, bin_set, x_name, y_name, tr, br, e_w)))

        name = eta_loop+"_"+pt_loop+"_resM_beta"
        hist_l.append(copy.deepcopy(hist2_maker(name, tr, "beta", "gbeta", e_w)))

      del tr
      del el


out_rf = ROOT.TFile(out_f, "RECREATE")
for x in hist_l:
  print x.GetName(), x.GetEntries()
  x.Write()
out_rf.Write()
out_rf.Close()
in_rf.Close()
