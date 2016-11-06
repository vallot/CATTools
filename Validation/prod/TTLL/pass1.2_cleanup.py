#!/usr/bin/env python

import sys, os
from os.path import isdir as isdir
from os.path import join as pathjoin
from multiprocessing import Pool, cpu_count
import imp

def hadd(d):
    d = d.rstrip('/')
    if not isdir(d): return 

    nqueue = 0
    for l in open(pathjoin(d, 'submit.jds')).readlines():
        if not l.startswith('queue'): continue
        if len(l.split()) == 1: nqueue = 1
        else: nqueue = int(l.split()[-1])

    rootFiles = []
    for x in os.listdir(d):
        if not x.endswith('.root') or not x.startswith("hist_"): continue
        x = pathjoin(d, x)
        if os.stat(x).st_size <= 200: continue
        rootFiles.append(x)

    if len(rootFiles) == nqueue:

        cmd = 'hadd -f %s.root ' % d
        cmd += ' '.join(rootFiles)

        print
        print "*"*40
        print "*", d, "is complete. Merge and clean up"
        print "*"*40
        print
        os.system(cmd)
        os.system("rm -f %s/job.tar.gz" % d)
        os.system("tar czf %s/log.tgz %s/*.log %s/*.err %s/*.txt" % (d, d, d, d))
        os.system("rm -f %s/*.log %s/*.err %s/*.txt" % (d, d, d))
        for x in rootFiles: os.system("rm -f %s" % (x))
    else:
        print
        print "+"*40
        print "+", d, "is incomplete. wait for the job completion"
        print "+"*40
        print

def jobStatus(d):
    if not os.path.exists("%s/condor.log" % d): return [], [], []
    jobs = {}
    lines = [x.strip() for x in open("%s/condor.log" % d).readlines()]
    for i, l in enumerate(lines):
        if 'Job submitted' in l:
            section = int(l.split()[1][1:-1].split('.')[1])
            jobs[section] = 'SUBMIT'
        elif 'Job terminated' in l:
            section = int(l.split()[1][1:-1].split('.')[1])
            retVal = int(lines[i+1].split()[-1][:-1])
            if retVal == 0: jobs[section] = 'DONE'
            else: jobs[section] = 'ERROR %d' % retVal

    jobsSub, jobsDone, jobsErr = [], [], []
    for section in jobs:
        if jobs[section] == 'SUBMIT': jobsSub.append(section)
        elif jobs[section] == 'DONE': jobsDone.append(section)
        else: jobsErr.append(section)
    return jobsSub, jobsDone, jobsErr

if __name__ == '__main__':
    pool = Pool(cpu_count())
    pass1Dir = "pass1"

    nSub = 0
    jobsFinished, jobsFailed = [], []
    for sample in os.listdir(pass1Dir):
        sample = pathjoin(pass1Dir, sample, 'nominal')
        if not isdir(sample): continue

        stat = jobStatus(sample)
        if len(stat[0])+len(stat[1])+len(stat[2]) == 0: continue

        nSub += len(stat[0])
        if len(stat[0])+len(stat[1])+len(stat[2]) == 0: continue
        elif len(stat[0])+len(stat[2]) == 0: jobsFinished.append(sample)
        elif len(stat[0]) == 0 and len(stat[2]) > 0: jobsFailed.append((sample, stat[2][:]))

    for sample in [pathjoin(pass1Dir, x) for x in os.listdir(pass1Dir) if isdir(pathjoin(pass1Dir, x))]:
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
        pool.apply_async(hadd, [sample])
        outFiles.append(sample+'.root')

    pool.close()
    pool.join()

    if len(jobsFailed) > 0:
        print "@@ There are some failed jobs"
        print "@@ Please resubmit following jobs"
        #if os.path.exists("pass1/resubmit"): os.system("rm -rf pass1/resubmit")
        if not os.path.exists("pass1/resubmit"): os.makedirs("pass1/resubmit")
        fsubmit = open("pass1/resubmit/submit.sh", "w")
        for sample, jobs in jobsFailed:
            print sample, jobs
            files = []
            for job in jobs:
                os.system("rm -f %s/hist_%03d.root" % (sample, job))
                process = imp.load_source("process", "%s/job_%03d_cfg.py" % (sample, job)).process
                for file in process.source.fileNames:
                    files.append(file)
            prefix = sample.replace("pass1/", "").replace("/", "_")
            frerun = open("pass1/resubmit/%s.txt" % prefix, "w")
            for f in files: print>>frerun, f
            frerun.close()
            frerun = open("pass1/resubmit/%s_cfg.py" % prefix, "w")
            frerun.write(process.dumpPython())
            frerun.close()

            print>>fsubmit, ("create-batch --jobName {0} --fileList {0}.txt --cfg {0}_cfg.py --maxFiles 1".format(prefix))

        fsubmit.close()

    if nSub > 0:
        print "@@ Done, wait for the %d jobs to finish" % nSub
    if nSub == 0 and len(jobsFailed) == 0:
        print "@@ Finished"
