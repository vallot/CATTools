#!/usr/bin/env python

from ROOT import *
from math import sqrt
import json

rootfileDir = "pass2/nominal"

hName = "eventsTTLL/%s/step%d/z_m_noveto"
fDY = TFile.Open(rootfileDir+"/Z__gamma_rightarrow_ll.root")
fEE = TFile(rootfileDir+"/DoubleEG.root")
fMM = TFile(rootfileDir+"/DoubleMuon.root")
fEM = TFile(rootfileDir+"/MuonEG.root")

dsIn = json.loads(open("pass2/dataset.json").read())
lumi = 1000*sum(x['lumi'] for x in dsIn['DoubleMuon']['subsamples'])

## Calculate efficiency factor kMM and kEE using the step1
hRD_ee = fEE.Get(hName % ("ee", 1))
hRD_mm = fMM.Get(hName % ("mm", 1))
hRD_em = fEM.Get(hName % ("em", 1))
hDY_ee = fDY.Get(hName % ("ee", 1))
hDY_mm = fDY.Get(hName % ("mm", 1))

bin1 = hRD_ee.FindBin(91-15+1e-5)
bin2 = hRD_ee.FindBin(91+15-1e-5)

nRD_ee_in = hRD_ee.Integral(bin1, bin2)
nRD_mm_in = hRD_mm.Integral(bin1, bin2)

kEE = sqrt(nRD_ee_in/nRD_mm_in)/2.
kMM = sqrt(nRD_mm_in/nRD_ee_in)/2.

print "kMM = ", kMM, "kEE = ", kEE

out = {
    "file":"Z__gamma_rightarrow_ll.root",
    "kMM":kMM, "kEE":kEE,
    "scale":{},
}

for step in range(1,6):
    hRD_ee = fEE.Get(hName % ("ee", step))
    hRD_mm = fMM.Get(hName % ("mm", step))
    hRD_em = fEM.Get(hName % ("em", step))
    hDY_ee = fDY.Get(hName % ("ee", step))
    hDY_mm = fDY.Get(hName % ("mm", step))

    nRD_ee_in = hRD_ee.Integral(bin1, bin2)
    nRD_mm_in = hRD_mm.Integral(bin1, bin2)
    nRD_em_in = hRD_em.Integral(bin1, bin2)
    nDY_ee_in = lumi*hDY_ee.Integral(bin1, bin2)
    nDY_mm_in = lumi*hDY_mm.Integral(bin1, bin2)

    nRD_ee_out = hRD_ee.Integral()-nRD_ee_in
    nRD_mm_out = hRD_mm.Integral()-nRD_mm_in
    nRD_em_out = hRD_em.Integral()-nRD_em_in
    nDY_ee_out = lumi*hDY_ee.Integral()-nDY_ee_in
    nDY_mm_out = lumi*hDY_mm.Integral()-nDY_mm_in

    r_ee = nDY_ee_out/nDY_ee_in
    r_mm = nDY_mm_out/nDY_mm_in

    nOutEst_ee = r_ee*(nRD_ee_in - nRD_em_in*kEE)
    nOutEst_mm = r_mm*(nRD_mm_in - nRD_em_in*kMM)

    scale_ee = nOutEst_ee/nDY_ee_out
    scale_mm = nOutEst_mm/nDY_mm_out

    print "DY estimation for step", step, "ee =",scale_ee, "mm =",scale_mm

    out["scale"]["eventsTTLL/ee/step%d" % step] = scale_ee
    out["scale"]["eventsTTLL/mm/step%d" % step] = scale_mm

f = open("pass2/scaler_DY.json", "w")
f.write(json.dumps(out, indent=4, sort_keys=True))
f.close()
