#!/usr/bin/env python

import sys, os
from os.path import isdir as isdir
from os.path import join as pathjoin
from multiprocessing import Pool, cpu_count
import imp

from CATTools.CommonTools.condorTools import jobStatus

prefix = "hist"
outDir = "pass1"

def hadd2(destDir, inputFiles):
    ## First split files by prefix
    inputFiles.sort(key = lambda x : int(x.split('_')[-1][:-5]))

    ## estimate total size and merge not to exceed 1Gbytes
    totalSize = sum([os.stat(f).st_size for f in inputFiles])
    fSize = 1024.*1024*1024 # 1GB
    nOutput = int(totalSize/fSize)+1

    if nOutput == 1:
        cmd = 'hadd -f %s.root ' % destDir
        cmd += ' '.join(inputFiles)
        os.system(cmd)
    else:
        nTotalFiles = len(inputFiles)
        nFiles = nTotalFiles/nOutput

        for i in range(nOutput):
            print i
            cmd = 'hadd -f %s_%d.root ' % (destDir, i)
            if i == nOutput-1: cmd += ' '.join(inputFiles[i*nFiles:])
            else: cmd += ' '.join(inputFiles[i*nFiles:(i+1)*nFiles])
            os.system(cmd)
    for x in inputFiles: os.system("rm -f %s" % (x))

def cleanup(d):
    d = d.rstrip('/')
    if not isdir(d): return 

    nqueue = 0
    for l in open(pathjoin(d, 'submit.jds')).readlines():
        if not l.startswith('queue'): continue
        if len(l.split()) == 1: nqueue = 1
        else: nqueue = int(l.split()[-1])

    rootFiles = []
    for x in os.listdir(d):
        if not x.endswith('.root') or not x.startswith(prefix): continue
        x = pathjoin(d, x)
        if os.stat(x).st_size <= 200: continue
        rootFiles.append(x)

    if len(rootFiles) == nqueue:
        print
        print "*"*40
        print "*", d, "is complete. Merge and clean up"
        print "*"*40
        print
        hadd2(d, rootFiles)
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
    pool = Pool(cpu_count())

    nSub = 0
    jobsFinished, jobsFailed = [], []
    for sample in os.listdir(outDir):
        sample = pathjoin(outDir, sample)
        if not isdir(sample): continue

        stat = jobStatus(sample)
        if len(stat[0])+len(stat[1])+len(stat[2]) == 0: continue

        nSub += len(stat[0])
        if len(stat[0])+len(stat[1])+len(stat[2]) == 0: continue
        elif len(stat[0])+len(stat[2]) == 0: jobsFinished.append(sample)
        elif len(stat[0]) == 0 and len(stat[2]) > 0: jobsFailed.append((sample, stat[2][:]))

    for sample in [pathjoin(outDir, x) for x in os.listdir(outDir) if isdir(pathjoin(outDir, x))]:
        for subsample in [pathjoin(sample, x) for x in os.listdir(sample) if isdir(pathjoin(sample, x))]:
            for syst in [pathjoin(subsample, x) for x in os.listdir(subsample) if isdir(pathjoin(subsample, x))]:
                stat = jobStatus(syst)
                if len(stat[0])+len(stat[1])+len(stat[2]) == 0: continue

                nSub += len(stat[0])
                if len(stat[0])+len(stat[1])+len(stat[2]) == 0: continue
                elif len(stat[0])+len(stat[2]) == 0: jobsFinished.append(syst)
                elif len(stat[0]) == 0 and len(stat[2]) > 0: jobsFailed.append((syst, stat[2][:]))

    outFiles = []
    for sample in jobsFinished:
        pool.apply_async(cleanup, [sample])
        outFiles.append(sample+'.root')

    pool.close()
    pool.join()

    if len(jobsFailed) > 0:
        print "@@ There are some failed jobs"
        print "@@ Please resubmit following jobs"
        if os.path.exists("%s/resubmit" % outDir): os.system("rm -rf %s/resubmit" % outDir)
        os.makedirs("%s/resubmit" % outDir)
        fsubmit = open("%s/resubmit/submit.sh" % outDir, "w")
        for sample, jobs in jobsFailed:
            print sample, jobs
            files = []
            for job in jobs:
                os.system("rm -f %s/%s_%03d.root" % (sample, prefix, job))
                process = imp.load_source("process", "%s/job_%03d_cfg.py" % (sample, job)).process
                for file in process.source.fileNames:
                    files.append(file)
            prefix = sample.replace("%s/" % outDir, "").replace("/", "_")
            frerun = open("%s/resubmit/%s.txt" % (outDir, prefix), "w")
            for f in files: print>>frerun, f
            frerun.close()
            frerun = open("%s/resubmit/%s_cfg.py" % (outDir, prefix), "w")
            frerun.write(process.dumpPython())
            frerun.close()

            print>>fsubmit, ("create-batch --jobName {0} --fileList {0}.txt --cfg {0}_cfg.py --maxFiles 1".format(prefix))

        fsubmit.close()
        os.system("chmod +x %s/resubmit/submit.sh" % outDir)

    if nSub > 0:
        print "@@ Done, wait for the %d jobs to finish" % nSub
    if nSub == 0 and len(jobsFailed) == 0:
        print "@@ Finished"
