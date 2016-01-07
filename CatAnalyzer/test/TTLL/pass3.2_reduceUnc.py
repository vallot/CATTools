#!/usr/bin/env python

from ROOT import *
import json
import sys, os
from math import sqrt
from multiprocessing import Pool

#gROOT.ProcessLine(".L submacros/combine.C+")
try:
    res = gSystem.CompileMacro("submacros/combine.C", "k")
    if res != 1:
        raise "Cannot compile"

    def runCombine(fNameCen, fNameUp, fNameDn, fNames, plotNames, combineBy):
        combine(fNameCen, fNameUp, fNameDn,
                vstring(fNames), vstring(plotNames), combineBy)

    def vstring(l):
        out = std.vector("string")()
        for x in l: out.push_back(std.string(x))
        return out

except:
    def runCombine(fNameCen, fNameUp, fNameDn, fNames, plotNames, combineBy):
        fcen = TFile(fNameCen)
        if fcen == None: return
        fcen.SetBit(TFile.kDevNull)

        print "@@ Processing", fNameCen

        fup = TFile(fNameUp, "recreate")
        fdn = TFile(fNameDn, "recreate")

        files = []
        for fName in fNames:
            files.append(TFile(fName))
            files[-1].SetBit(TFile.kDevNull)

        for pName in plotNames:
            hcen = fcen.Get(pName)
            if hcen == None: continue
            nbins = hcen.GetNbinsX()

            diffs = [[]]*(nbins+2)
            for f in files:
                h = f.Get(pName)
                for b in range(nbins+2):
                    diffs[b].append(h.GetBinContent(b)-hcen.GetBinContent(b))

        rdir = os.path.dirame(pName)
        dup = fup.GetDirectory(rdir)
        if dup == None:
            fup.mkdir(rdir)
            dup = fup.GetDirectory(rdir)
        dup.cd()
        hup = hcen.Clone()

        ddn = fdn.GetDirectory(rdir)
        if ddn == None:
            fdn.mkdir(rdir)
            ddn = fdn.GetDirectory(rdir)
        ddn.cd()
        hdn = hcen.Clone()

        if combineBy == "hessian":
            for b in range(nbins+2):
                dqsqr = sum([dyi[b]**2 for dyi in diffs])
                hup.AddBinContent(b, sqrt(dysqr))
                hdn.AddBinContent(b, -sqrt(dysqr))
        elif combineBy == "envelope":
            for b in range(nbins+2):
                dymax = max([dyi[b] for dyi in diffs])
                dymin = min([dyi[b] for dyi in diffs])
                hup.AddBinContent(b, dymax)
                hdn.AddBinContent(b, dymin)
        elif combineBy == "gauss":
            for b in range(nbins+2):
                n = len(diffs)
                dqsqr = sum([dyi[b]**2 for dyi in diffs])
                hup.AddBinContent(b, sqrt(dysqr)/n)
                hdn.AddBinContent(b, -sqrt(dysqr)/n)

        dup.cd()
        hup.Write("", TObject.kOverwrite)

        ddn.cd()
        hdn.Write("", TObject.kOverwrite)

        print "@@ Clean up opened files"
        for f in files:
            gROOT.GetListOfFiles().Remove(f)
        fcen.Close()
        fup.Close()
        fdn.CLose()

        print "@@ Finished", fNameCen

if __name__ == '__main__':
    ## Build sample list to combine uncertainty
    ## Sort by file category/fileName and results will be reduce to up and dn
    plotsJS = json.loads(open("pass2/plots.json").read())
    plotNames = [p["name"] for p in plotsJS["plots"]]
    uncToReduce = {}
    samplesJS = json.loads(open("pass2/samples.json").read())
    for fPath in samplesJS:
        if samplesJS[fPath]["type"] == "data": continue
        if 'central' in fPath: continue

        cat, id, fName = fPath.split('/')[1:]
        if id in ('up', 'dn'): continue

        #fName = fName.replace("pass2", "pass3/hists")
        #uncToReduce[cat]["files"].append(fName)
        key = (cat, fName)
        if key not in uncToReduce: uncToReduce[key] = []
        uncToReduce[key].append(fPath.replace("pass2", "pass3/hists"))

    ## Start to loop over all of them and do the reduction.
    pool = Pool(20)
    for cat, fName in uncToReduce:
        ## Prepare output
        if not os.path.exists("pass3/hists/%s/up" % cat): os.makedirs("pass3/hists/%s/up" % cat)
        if not os.path.exists("pass3/hists/%s/dn" % cat): os.makedirs("pass3/hists/%s/dn" % cat)

        fNameCen = "pass3/hists/central/"+fName
        fNameUp = "pass3/hists/%s/up/%s" % (cat, fName)
        fNameDn = "pass3/hists/%s/dn/%s" % (cat, fName)

        if   cat == 'gen_PDF'  : combineBy = "hessian"
        elif cat == 'gen_scale': combineBy = "envelope"
        else:
            print "!!!! Combine method was not defined for this category,", cat
            print "!!!! Skip this one..."
            continue

        pool.apply_async(runCombine,
                         [fNameCen, fNameUp, fNameDn,
                          uncToReduce[(cat, fName)], plotNames, combineBy])
    pool.close()
    pool.join()

