#!/usr/bin/env python
import subprocess
#from CATTools.CatAnalyzer.histoHelper import *

inputfile="catTuple.root"



import ROOT #, CATTools.CatAnalyzer.CMS_lumi, json, os,copy
ROOT.gROOT.SetBatch(True)

ROOT.gSystem.Load("libFWCoreFWLite.so");
ROOT.gSystem.Load("libDataFormatsFWLite.so");
ROOT.AutoLibraryLoader.enable();



from DataFormats.FWLite import Handle, Events
from matching import Matching
events = Events(inputfile)

genParticleHandle, genParticleLabel = Handle("vector<reco::GenParticle>"), "genParticles"

V1Handle, V1Label = Handle("vector<cat::SecVertex>"), "catSecVertexs"
V2Handle, V2Label = Handle("vector<cat::SecVertexV2>"), "catSecVertexsV2"
V3Handle, V3Label = Handle("vector<cat::SecVertexV2>"), "catSecVertexsV3"
V4Handle, V4Label = Handle("vector<cat::SecVertexV2>"), "catSecVertexsV4"

output=ROOT.TFile("output.root","RECREATE")
profile1 = ROOT.TProfile("eff_prod","Efficiency of each J/#psi reconstruction modules; Module type ; Efficiency",4,1,5)
profile2 = ROOT.TProfile("jpsi_averge","Average J/#psi mulplicity of each J/#psi reconstruction modules; Module type ; Average",4,1,5)

profile1.GetXaxis().SetBinLabel(1, "Lepton+Lepton")
profile1.GetXaxis().SetBinLabel(2, "Lepton+generalTracks")
profile1.GetXaxis().SetBinLabel(3, "miniAOD's SecVtx")
profile1.GetXaxis().SetBinLabel(4, "pfCandidates+pfCandidates")

profile2.GetXaxis().SetBinLabel(1, "Lepton+Lepton")
profile2.GetXaxis().SetBinLabel(2, "Lepton+generalTracks")
profile2.GetXaxis().SetBinLabel(3, "miniAOD's SecVtx")
profile2.GetXaxis().SetBinLabel(4, "pfCandidates+pfCandidates")
print events.size()
for i,event in enumerate(events) :
  event.getByLabel( genParticleLabel, genParticleHandle)
  event.getByLabel( V1Label, V1Handle)
  event.getByLabel( V2Label, V2Handle)
  event.getByLabel( V3Label, V3Handle)
  event.getByLabel( V4Label, V4Handle)
  jpsi_idx = 0
  genJpsis =[]
  for gen in genParticleHandle.product() :
    if ( abs(gen.pdgId()) == 443 and gen.numberOfDaughters()==2 and (  abs(gen.daughter(0).pdgId())==13 or abs(gen.daughter(0).pdgId()) == 1 ) ) :
      jpsi_idx = jpsi_idx +1 
      genJpsis.append( gen )
      print "genJpsi pt : %f , eta %f, phi %f"%(gen.pt(), gen.eta(), gen.phi())
  if jpsi_idx ==0 : continue
  V1 = V1Handle.product()
  V2 = V2Handle.product()
  V3 = V3Handle.product()
  V4 = V4Handle.product()
  mV1 = Matching( genJpsis, V1)
  mV2 = Matching( genJpsis, V2)
  mV3 = Matching( genJpsis, V3)
  mV4 = Matching( genJpsis, V4)

  print "genJpsi : %d // Size of J/psi v1  : %d  // j/psi v2 : %d  // j/psi v3 :  %d  // j/psi v4 :  %d"%(jpsi_idx, len(V1), len(V2),len(V3), len(V4))
  if ( len(mV1) > 0 ) : profile1.Fill(1,1)
  else : profile1.Fill(1,0)
  if ( len(mV2) > 0 ) : profile1.Fill(2,1)
  else : profile1.Fill(2,0)
  if ( len(mV3) > 0 ) : profile1.Fill(3,1)
  else : profile1.Fill(3,0)
  if ( len(mV4) > 0 ) : profile1.Fill(4,1)
  else : profile1.Fill(4,0)
  #print "genJpsi : %d // size of Matched J/psi  // j/psi v2 : %d  // j/psi v3 :  %d  // j/psi v4 :  %d"%(jpsi_idx, len(matchV2),len(matchV3), len(matchV4))
  profile2.Fill(1,V1.size())
  profile2.Fill(2,V2.size())
  profile2.Fill(3,V3.size())
  profile2.Fill(4,V4.size())

output.Write()
output.Close()
