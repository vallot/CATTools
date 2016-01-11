#!/usr/bin/env python

import json
import os
from multiprocessing import Pool, cpu_count

outbase = "pass3/hists"
try:
    from ROOT import *
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
    import imp
    runCombine = imp.load_source("combine", "submacros/combine.py").combine

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

        key = (cat, fName)
        if key not in uncToReduce: uncToReduce[key] = []
        uncToReduce[key].append(fPath.replace("pass2", outbase))

    ## Start to loop over all of them and do the reduction.
    pool = Pool(cpu_count())
    for cat, fName in uncToReduce:
        ## Prepare output
        if not os.path.exists("%s/%s/up" % (outbase, cat)): os.makedirs("%s/%s/up" % (outbase, cat))
        if not os.path.exists("%s/%s/dn" % (outbase, cat)): os.makedirs("%s/%s/dn" % (outbase, cat))

        fNameCen = "%s/central/%s" % (outbase, fName)
        fNameUp = "%s/%s/up/%s" % (outbase, cat, fName)
        fNameDn = "%s/%s/dn/%s" % (outbase, cat, fName)

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

