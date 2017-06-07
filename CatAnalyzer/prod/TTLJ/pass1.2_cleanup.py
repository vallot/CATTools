#!/usr/bin/env python

import sys, os
from os.path import isdir as isdir
from os.path import join as pathjoin
from multiprocessing import Pool, cpu_count
import imp

from CATTools.CommonTools.condorTools import jobStatus
from CATTools.CommonTools.haddsplit import haddsplit

def cleanup(srcDir, outDir):
    d = srcDir.rstrip('/')
    if not isdir(d): return 

    nqueue = 0
    for l in open(pathjoin(d, 'submit.jds')).readlines():
        if not l.startswith('queue'): continue
        if len(l.split()) == 1: nqueue = 1
        else: nqueue = int(l.split()[-1])

    rootFiles = []
    for x in os.listdir(d):
        if not x.endswith('.root'): continue
        x = pathjoin(d, x)
        if os.stat(x).st_size <= 200: continue
        rootFiles.append(x)

    if len(rootFiles) == nqueue:
        print
        print "*"*40
        print "*", d, "is complete. Merge and clean up"
        print "*"*40
        print
        isOK = haddsplit(srcDir, rootFiles, outDir)
        if isOK:
            os.system("rm -f %s/job.tar.gz" % d)
            os.system("tar czf %s/log.tgz %s/*.log %s/*.err %s/*.txt" % (d, d, d, d))
            os.system("rm -f %s/*.log %s/*.err %s/*.txt" % (d, d, d))
    else:
        print
        print "+"*40
        print "+", d, "is incomplete. wait for the job completion"
        print "+"*40
        print

if __name__ == '__main__':
    srcBase = "pass1"
    outBase = "pass2"

    nSub = 0
    jobsFinished, jobsFailed = [], []
    for sample, dirName, files in os.walk(srcBase):
        if '.create-batch' not in files or 'job.tar.gz' not in files: continue
        if 'condor.log' not in files:
            print sample, " incomplete submission"
            continue

        stat = jobStatus(sample)
        if len(stat[0])+len(stat[1])+len(stat[2]) == 0:
            print sample, " without jobs"
            continue

        nSub += len(stat[0])
        if len(stat[0])+len(stat[1])+len(stat[2]) == 0: continue
        elif len(stat[0])+len(stat[2]) == 0: jobsFinished.append(sample)
        elif len(stat[0]) == 0 and len(stat[2]) > 0: jobsFailed.append((sample, stat[2][:]))

    pool = Pool(cpu_count())
    for sample in jobsFinished:
        outDir = outBase+'/'+sample.split('/',1)[-1]
        pool.apply_async(cleanup, [sample, outDir])
    pool.close()
    pool.join()

    if len(jobsFailed) > 0:
        print "@@ There are some failed jobs"
        print "@@ Please resubmit following jobs"
        if os.path.exists("%s/resubmit" % srcBase): os.system("rm -rf %s/resubmit" % srcBase)
        os.makedirs("%s/resubmit" % srcBase)
        fsubmit = open("%s/resubmit/submit.sh" % srcBase, "w")
        for sample, jobs in jobsFailed:
            print sample, jobs
            files = []
            for job in jobs:
                os.system("rm -f %s/*_%03d.root" % (sample, job))
                process = imp.load_source("process", "%s/job_%03d_cfg.py" % (sample, job)).process
                for file in process.source.fileNames:
                    files.append(file)
            prefix = sample.replace("%s/" % srcBase, "").replace("/", "_")
            frerun = open("%s/resubmit/%s.txt" % (srcBase, prefix), "w")
            for f in files: print>>frerun, f
            frerun.close()
            frerun = open("%s/resubmit/%s_cfg.py" % (srcBase, prefix), "w")
            frerun.write(process.dumpPython())
            frerun.close()

            print>>fsubmit, ("create-batch --jobName {0} --fileList {0}.txt --cfg {0}_cfg.py --maxFiles 1".format(prefix))

        fsubmit.close()
        os.system("chmod +x %s/resubmit/submit.sh" % srcBase)

    if nSub > 0:
        print "@@ Done, wait for the %d jobs to finish" % nSub
    if nSub == 0 and len(jobsFailed) == 0:
        print "@@ Finished"
