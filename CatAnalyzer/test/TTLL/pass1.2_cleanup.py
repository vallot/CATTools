#!/usr/bin/env python

import sys, os
from os.path import isdir as isdir
from os.path import join as pathjoin
from multiprocessing import Pool, cpu_count

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

if __name__ == '__main__':
    pool = Pool(cpu_count())
    pass1Dir = "pass1"
    outFiles = []
    for sample in os.listdir(pass1Dir):
        sample = pathjoin(pass1Dir, sample)
        if not isdir(sample): continue

        pool.apply_async(hadd, [pathjoin(sample,'central')])
        outFiles.append(pathjoin(sample,'central.root'))

        for category in os.listdir(sample):
            if category == 'central': continue
            category = pathjoin(sample, category)
            if not isdir(category): continue

            for direction in os.listdir(category):
                direction = pathjoin(category, direction)
                if not isdir(direction): continue
                pool.apply_async(hadd, [direction])
                outFiles.append(direction+'.root')

    pool.close()
    pool.join()

    print "@@ Done"
