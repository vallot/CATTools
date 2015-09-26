#!/usr/bin/env python

from ROOT import *

f = TFile("histtest.root")

genLevels = [
    ("FullParton", kBlack),
    ("FiducialParton", kBlue),
    ("CommonParton", kRed),
    ("CommonParticle", kRed),
    ("Particle", kMagenta),
]

plots = [
    "SL_topPt", "SL_topPtTtbarSys", "SL_topY", "SL_ttbarDelPhi", "SL_topPtLead", "SL_topPtSubLead", "SL_ttbarPt", "SL_ttbarMass",
    "DL_topPt", "DL_topPtTtbarSys", "DL_topY", "DL_ttbarDelPhi", "DL_topPtLead", "DL_topPtSubLead", "DL_ttbarPt", "DL_ttbarMass",
]

objs = []

for p in plots:
    hists = []
    c = TCanvas("c%s" % p, p, 600, 500)
    c.SetRightMargin(0.225)

    for i, (dirName, color) in enumerate(genLevels):
        h = f.Get("ana/%s/%s" % (dirName, p))
        if h == None:
            print "Cannot find %s/%s" % (dirName, p)
            continue
        h = h.Clone()
        h.SetName(dirName)
        h.SetLineColor(color)
        h.SetMinimum(0)

        for i in range(1,h.GetNbinsX()+1):
            h.SetBinContent(i, h.GetBinContent(i)/h.GetBinWidth(i))

        if len(hists) == 0: h.Draw()
        else: h.Draw("sames")

        hists.append(h)
    objs.append([c, hists])

for c, hists in objs:
    c.Update()
    height = (1-c.GetTopMargin())*min(0.2, 1./len(hists))
    for i, h in enumerate(hists):
        st = h.FindObject("stats")
        if st == None: continue

        y = 1-i*height-c.GetTopMargin()-i*0.01
        st.SetY2NDC(y)
        st.SetY1NDC(y-height)
        st.SetLineColor(h.GetLineColor())
    c.Modified()
    c.Update()

    c.Print("%s.png" % c.GetName())
