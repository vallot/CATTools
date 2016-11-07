#!/usr/bin/env python

from math import sqrt

## Start to print cut flow
def printCutflow(cutflow):
    step = cutflow["step"]
    nstep = cutflow["nstep"]
    count = cutflow["count"]
    error = cutflow["error"]
    modes = count.keys()

    ## Set the field width
    print count[modes[0]].keys()
    fw0 = max(4, max([len(x) for x in count[modes[0]].keys()]))
    fws = [len(x) for x in step]
    fwtot = sum(fws)+len(fws)+12+fw0

    for mode in count.keys():
        print "="*((fwtot-12)/2), "Cutflow for", mode, "="*((fwtot-12)/2)
        print " "*fws[0], "|",
        print " | ".join(step)
        tfmt = "%"+str(fw0)+"s |"
        count_bkg = [0.]*nstep
        errsq_bkg = [0.]*nstep
        for x in count[mode]:
            if x == "Data" or ('t_bar_t' in x and 'Others' not in x): continue
            c, e = count[mode][x], error[mode][x]
            print tfmt % x,
            print " | ".join([("%"+str(fws[i]*2/3)+".1f+-%"+str(fws[i]/3)+".1f") % (c[i], e[i]) for i in range(nstep)])
            for i in range(nstep):
                count_bkg[i] += c[i]
                errsq_bkg[i] += e[i]**2
        print "-"*fwtot
        count_sig = [0.]*nstep
        errsq_sig = [0.]*nstep
        for x in count[mode]:
            if 't_bar_t' not in x or 'Others' in x: continue
            c, e = count[mode][x], error[mode][x]
            print tfmt % x,
            print " | ".join([("%"+str(fws[i]*2/3)+".1f+-%"+str(fws[i]/3)+".1f") % (c[i], e[i]) for i in range(nstep)])
            for i in range(nstep):
                count_sig[i] += c[i]
                errsq_sig[i] += e[i]**2
        print "-"*fwtot
        print tfmt % "All Signal",
        print " | ".join([("%"+str(fws[i]*2/3)+".1f+-%"+str(fws[i]/3)+".1f") % (count_sig[i], sqrt(errsq_sig[i])) for i in range(nstep)])
        print tfmt % "All Bkg",
        print " | ".join([("%"+str(fws[i]*2/3)+".1f+-%"+str(fws[i]/3)+".1f") % (count_bkg[i], sqrt(errsq_bkg[i])) for i in range(nstep)])
        print tfmt % "All MC",
        print " | ".join([("%"+str(fws[i]*2/3)+".1f+-%"+str(fws[i]/3)+".1f") % (count_sig[i]+count_bkg[i], sqrt(errsq_sig[i]+errsq_bkg[i])) for i in range(nstep)])
        print "-"*fwtot
        print tfmt % "Data",
        print " | ".join([("%"+str(fws[i]-3)+"d   ") % count[mode]["Data"][i] for i in range(nstep)])
        print "="*fwtot
        print

if __name__ == '__main__':
    import json
    cutflow = json.loads(open("pass2/cutflow.json").read())
    printCutflow(cutflow)
