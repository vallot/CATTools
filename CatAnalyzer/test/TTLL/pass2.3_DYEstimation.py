#!/usr/bin/python

from ROOT import *
import math
datalumi = 2.11*1000

rootfileDir = "pass2/central"

stepList = [1,2,3,4,5]

fDY = TFile.Open(rootfileDir+"/Z__gamma_rightarrow_ll.root")
fEM = TFile.Open(rootfileDir+"/MuonEG.root")
fMM = TFile.Open(rootfileDir+"/DoubleMuon.root")
fEE = TFile.Open(rootfileDir+"/DoubleEG.root")

rd_ee_out_loose = fEE.Get("/ttll/ee/step1/z_m").Integral()
rd_ee_in_loose = fEE.Get("/ttll/ee/step1/z_m_noveto").Integral()-rd_ee_out_loose  

rd_mm_out_loose = fMM.Get("/ttll/mm/step1/z_m").Integral()
rd_mm_in_loose = fMM.Get("/ttll/mm/step1/z_m_noveto").Integral()-rd_mm_out_loose

kEE = math.sqrt(rd_ee_in_loose/rd_mm_in_loose)/2. 
kMM = math.sqrt(rd_mm_in_loose/rd_ee_in_loose)/2. 

print "kMM = ", kMM, "kEE = ", kEE

for i,step in enumerate(stepList):
  dy_ee_out = fDY.Get("/ttll/ee/step%s/z_m" % step).Integral()*datalumi
  dy_ee_in = fDY.Get("/ttll/ee/step%s/z_m_noveto" % step).Integral()*datalumi-dy_ee_out  

  dy_mm_out = fDY.Get("/ttll/mm/step%s/z_m" % step).Integral()*datalumi
  dy_mm_in = fDY.Get("/ttll/mm/step%s/z_m_noveto" % step).Integral()*datalumi-dy_mm_out 

  rd_ee_out = fEE.Get("/ttll/ee/step%s/z_m" % step).Integral()
  rd_ee_in = fEE.Get("/ttll/ee/step%s/z_m_noveto" % step).Integral()-rd_ee_out  

  rd_mm_out = fMM.Get("/ttll/mm/step%s/z_m" % step).Integral()
  rd_mm_in = fMM.Get("/ttll/mm/step%s/z_m_noveto" % step).Integral()-rd_mm_out 

  rd_em_out = fMM.Get("/ttll/em/step%s/z_m" % step).Integral()
  rd_em_in = fMM.Get("/ttll/em/step%s/z_m_noveto" % step).Integral()-rd_em_out 

  rDY_ee = dy_ee_out/dy_ee_in 
  rDY_mm = dy_mm_out/dy_mm_in 

  nOutEst_ee = rDY_ee*(rd_ee_in - rd_em_in*kEE)
  nOutEst_mm = rDY_mm*(rd_mm_in - rd_em_in*kMM)

  scale_ee = nOutEst_ee/dy_ee_out
  scale_mm = nOutEst_mm/dy_mm_out

  print "DY estimation for step", step, "ee =",scale_ee, "mm =",scale_mm

