#!/usr/bin/env python

from ROOT import *
from array import array

f = TFile("hist.root")
d = f.Get("ttbarNoTau")

genLevels = [
    ("FullParton", kBlack),
    ("FiducialParton", kBlue),
    ("CommonParton", kRed),
    ("CommonParticle", kMagenta),
    ("Particle", kGreen+1),
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
        h = d.Get("%s/%s" % (dirName, p))
        if h == None:
            print "Cannot find %s/%s" % (dirName, p)
            continue
        h = h.Clone()
        h.SetName(dirName)
        h.SetLineColor(color)
        h.SetMinimum(0)

        stats = array('d', [0,0,0,0])
        h.GetStats(stats)
        nEntry = h.GetEntries()
        for i in range(1,h.GetNbinsX()+1):
            h.SetBinContent(i, h.GetBinContent(i)/h.GetBinWidth(i))
        h.PutStats(stats)
        h.SetEntries(nEntry)

        if len(hists) == 0: h.Draw("hist")
        else: h.Draw("sameshist")

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

for p in plots:
    cCommon = TCanvas("cCommon%s" % p, p, 500, 500)

    hCommon = d.Get("CommonParticle/resp_%s" % p)
    hCommon.Draw("COLZ")

    cAll = TCanvas("cAll%s" % p, p, 500, 500)

    hAll = d.Get("Particle/resp_%s" % p)
    hAll.Draw("COLZ")

    objs.extend([cCommon, hCommon, cAll, hAll])

    cCommon.Print("%s.png" % cCommon.GetName())
    cAll.Print("%s.png" % cAll.GetName())

for p in plots:
    hParton = d.Get("FullParton/%s" % p)
    hPseudo = d.Get("Particle/%s" % p)

    vals = []
    for b in range(1,hParton.GetNbinsX()+1):
        nParton = hParton.GetBinContent(b)
        nPseudo = hPseudo.GetBinContent(b)
        #c = (nParton/hParton.Integral())/(nPseudo/hPseudo.Integral())
        c = (nParton)/(nPseudo)

        binX = hParton.GetXaxis().GetBinLowEdge(b)
        binW = hParton.GetXaxis().GetBinWidth(b)

        vals.append("%.3f" % c)

    print '"%s":[%s]' % (p, (",".join(vals)))
