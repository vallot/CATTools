#!/usr/bin/env python

## Start to print cut flow
def printCutflow(cutflow):
    step = cutflow["step"]
    nstep = cutflow["nstep"]
    count = cutflow["count"]
    modes = count.keys()

    ## Set the field width
    fw0 = max(4, max([len(x) for x in cutflow["count"][modes[0]]]))
    fws = [len(x) for x in step]
    fwtot = sum(fws)+len(fws)+12+fw0

    for mode in count.keys():
        print "="*((fwtot-12)/2), "Cutflow for", mode, "="*((fwtot-12)/2)
        print " "*fws[0], "|",
        print " | ".join(step)
        tfmt = "%"+str(fw0)+"s |"
        count_bkg = [0.]*nstep
        for x in count[mode]:
            if x == "Data" or 't_bar_t' in x: continue
            print tfmt % x,
            print " | ".join([("%"+str(fws[i])+".2f") % count[mode][x][i] for i in range(nstep)])
            for i in range(nstep): count_bkg[i] += count[mode][x][i]
        print "-"*fwtot
        count_sig = [0.]*nstep
        for x in count[mode]:
            if 't_bar_t' not in x: continue
            print tfmt % x,
            print " | ".join([("%"+str(fws[i])+".2f") % count[mode][x][i] for i in range(nstep)])
            for i in range(nstep): count_sig[i] += count[mode][x][i]
        print "-"*fwtot
        print tfmt % "All Signal",
        print " | ".join([("%"+str(fws[i])+".2f") % count_sig[i] for i in range(nstep)])
        print tfmt % "All Bkg",
        print " | ".join([("%"+str(fws[i])+".2f") % count_bkg[i] for i in range(nstep)])
        print tfmt % "All MC",
        print " | ".join([("%"+str(fws[i])+".2f") % (count_sig[i] + count_bkg[i]) for i in range(nstep)])
        print "-"*fwtot
        print tfmt % "Data",
        print " | ".join([("%"+str(fws[i]-3)+"d   ") % count[mode]["Data"][i] for i in range(nstep)])
        print "="*fwtot
        print

if __name__ == '__main__':
    import json
    cutflow = json.loads(open("pass2/cutflow.json").read())
    printCutflow(cutflow)
