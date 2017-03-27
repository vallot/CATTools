#!/usr/bin/env python
## merge.py - produce reduced set of histograms to go into final fit

import json
from ROOT import *
from array import array
from math import sqrt
import os

outdir = "hists"
if not os.path.exists(outdir): os.mkdir(outdir)
if not os.path.exists("pass2"): os.system("ln -s ../pass2")

module = "eventsTTLJ"
modes = ("el", "mu")
#step = "step5d" ## which cut step to start (step5d for nJet4)
step = "step6" ## which cut step to start (step6 for nBjet1)
hists = {
    'bjets_n':{
        'title':'bjets_n;b-jet multiplicity;Events',
        'bins':range(7),
    },
    'met_pt':{
        'title':'met_pt;MET (GeV);Events',
        'bins':[0,20,30,40,50,70,100,150,200,300,400,600,1000,],
    },
    'event_st':{
        'title':'event_st;#Sigma p_{T} (GeV);Events/10 GeV',
        'bins':range(0,1001,10),
    },
    'event_mT':{
        'title':'event_mT;Transverse mass of lepton+MET (GeV);Events/5 GeV',
        'bins':range(0,501,5),
    },
    'event_mlj':{
        'title':'event_m3;M3 (GeV);Events/5 GeV',
        'bins':range(0,501,5),
    },
    'event_m3':{
        'title':'event_m3;M3 (GeV);Events/5 GeV',
        'bins':range(0,501,5),
    },
}
cats = {
    'data':{
        'titles':["SingleElectron", "SingleMuon"],
        'onlyfor':["el", "mu"],
    },
    'ttbar':{
        'titles':['t#bar{t}+Jets#rightarrow l^{#pm}'],
    },
    'WJet':{
        'titles':['W+Jets'],
    },
    'bkg':{
        'titles':["t#bar{t}+Jets:Others", "SingleTop",
                  "Dibosons", "Tribosons", "Z/#gamma#rightarrow ll",],
    },
}
rescale = {
    "t#bar{t}+Jets#rightarrow l^{#pm}":0.3/0.428, ## gen filter to select mu/el channels
    "t#bar{t}+Jets:Others":(1-0.3/0.428), ## gen filter to remove mu/el channels
}
dataset = json.loads(open("pass2/dataset.json").read())
for item in cats.itervalues():
    item['files'] = [TFile(dataset[title]['hist']) for title in item['titles']]

## Start to build histograms
for name, item in cats.iteritems():
    fout = TFile("%s/%s.root" % (outdir, name), "recreate")
    item['file'] = fout
    item['hist'] = {}

    for mode in modes:
        dout = fout.mkdir(mode)
        for hName, hItem in hists.iteritems():
            bins = hItem['bins']

            sumw = [0.]*(len(bins)+1)
            err2 = sumw[:]

            hout = TH1F(hName, hItem['title'], len(bins)-1, array('d', bins))
            item['hist'][hName] = hout
            hout.SetDirectory(dout)

            for iin, fin in enumerate(item['files']):
                if 'onlyfor' in item and item['onlyfor'][iin] != mode: continue

                scale = 1
                title = item['titles'][iin]
                if title in rescale: scale *= rescale[title]
                hPath = "%s/%s/%s/%s" % (module, mode, step, hName)
                hin = fin.Get(hPath)
                for ibin in range(hin.GetNbinsX()+2):
                    jbin = min(len(sumw)-1, hout.FindBin(hin.GetXaxis().GetBinCenter(ibin)))
                    sumw[jbin] += scale*hin.GetBinContent(ibin)
                    err2[jbin] += scale*(hin.GetBinError(ibin)**2)

            for ibin, (y, err2) in enumerate(zip(sumw, err2)):
                hout.SetBinContent(ibin, y)
                hout.SetBinError(ibin, sqrt(err2))
        dout.Write()

    fout.Write()

