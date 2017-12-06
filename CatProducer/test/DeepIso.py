#! /usr/bin/env python
import os
import glob
import math

import ROOT
from ROOT import TLorentzVector
# Import stuff from FWLite
import sys
from DataFormats.FWLite import Events, Handle

files = ["catTuple.root"]
events = Events (files)

# Create the output file. 
#f = ROOT.TFile("hist_QCD_1event.root", "recreate")
#f.cd()

f_txt = open("qcd_iso.txt",'w')

chIsoHist = ROOT.TH2F('chIsoHist', "Charged Hadron", 25,-0.5,0.5,25,-0.5,0.5) 
nhIsoHist = ROOT.TH2F('nhIsoHist', "Neutron Hadron", 25,-0.5,0.5,25,-0.5,0.5) 
phIsoHist = ROOT.TH2F('phIsoHist', "Photon",         25,-0.5,0.5,25,-0.5,0.5) 
IsoHist = ROOT.TH2F('IsoHist', "Iso",         25,-0.5,0.5,25,-0.5,0.5) 

muonHandle = Handle( "std::vector<cat::Muon>" )
muonLabel = ("catMuons","","CAT")

for event in events:
  event.getByLabel (muonLabel, muonHandle)
  muons = muonHandle.product()
  #print "size = ",  len( muons )
  found = False
  for imu in range( 0, len( muons ) ) :
   
    passMuon = False
    passMuon = muons[imu].pt() > 20 and abs(muons[imu].eta()) < 2.1 and muons[imu].isTightMuon()

    if not passMuon:
      continue

    chIsoCands = muons[imu].chIsoCandidates()
    nhIsoCands = muons[imu].nhIsoCandidates()
    phIsoCands = muons[imu].phIsoCandidates()
    isoChIsoSum = 0
    isoNhIsoSum = 0
    isoPhIsoSum = 0
 
    Ncands = len(chIsoCands) + len(nhIsoCands) + len(phIsoCands) 
    f_txt.write(str(Ncands))
    f_txt.write('\n')

    p4muon = TLorentzVector(muons[imu].px(), muons[imu].py(), muons[imu].pz(), muons[imu].energy())

    f_txt.write(' ')
    printp4muon = '{:>17} {:>17} {:>17} {:>17}\n'.format(p4muon.Px(), p4muon.Py(), p4muon.Pz(), p4muon.E())
    f_txt.write(printp4muon)

    for i in range(0, len(chIsoCands)):
      p4cand = TLorentzVector(chIsoCands[i].px(), chIsoCands[i].py(), chIsoCands[i].pz(), chIsoCands[i].energy())

      printp4cand = ' {:>17} {:>17} {:>17} {:>17}\n'.format(p4cand.Px(), p4cand.Py(), p4cand.Pz(), p4cand.E())
      f_txt.write(printp4cand)
 
      deta = p4cand.Eta()-p4muon.Eta()
      dphi = p4muon.DeltaPhi(p4cand)
      chIsoHist.Fill(deta, dphi, p4cand.Pt())
      IsoHist.Fill(deta, dphi, p4cand.Pt())
      #print "charged= ", deta, " " , dphi, " " ,p4cand.Pt() 
      isoChIsoSum += chIsoCands[i].pt()
    for i in range(0, len(nhIsoCands)):
      p4cand = TLorentzVector(nhIsoCands[i].px(), nhIsoCands[i].py(), nhIsoCands[i].pz(), nhIsoCands[i].energy())
 
      printp4cand = ' {:>17} {:>17} {:>17} {:>17}\n'.format(p4cand.Px(), p4cand.Py(), p4cand.Pz(), p4cand.E())
      f_txt.write(printp4cand)
 
      deta = p4cand.Eta()-p4muon.Eta()
      dphi = p4muon.DeltaPhi(p4cand)
      nhIsoHist.Fill(deta, dphi, p4cand.Pt())
      IsoHist.Fill(deta, dphi, p4cand.Pt())
      #print "neutral= ", deta, " " , dphi, " " ,p4cand.Pt()
      isoNhIsoSum += nhIsoCands[i].pt()
    for i in range(0, len(phIsoCands)):
      p4cand = TLorentzVector(phIsoCands[i].px(), phIsoCands[i].py(), phIsoCands[i].pz(), phIsoCands[i].energy())

      printp4cand = ' {:>17} {:>17} {:>17} {:>17}\n'.format(p4cand.Px(), p4cand.Py(), p4cand.Pz(), p4cand.E())
      f_txt.write(printp4cand)

      deta = p4cand.Eta()-p4muon.Eta()
      dphi = p4muon.DeltaPhi(p4cand)
      phIsoHist.Fill(deta, dphi, p4cand.Pt())
      IsoHist.Fill(deta, dphi, p4cand.Pt())
      #print "photon= ", deta, " " , dphi, " " ,p4cand.Pt()
      isoPhIsoSum += phIsoCands[i].pt()

    found = True
    #print isoChIsoSum , " " , isoNhIsoSum , " " , isoPhIsoSum 
    #print muons[imu].chargedHadronIso(), " " , muons[imu].neutralHadronIso(), " " , muons[imu].photonIso()

  #if found: 
  #  break

f_txt.close()

#f.cd()
#f.Write()
#f.Close()

