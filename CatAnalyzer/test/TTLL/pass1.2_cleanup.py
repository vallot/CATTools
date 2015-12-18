#!/usr/bin/env python

import sys, os
from os.path import isdir as isdir
from os.path import join as pathjoin

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
        if not x.endswith('.root'): continue
        x = pathjoin(d, x)
        if os.stat(x).st_size <= 200: continue
        rootFiles.append(x)

    if len(rootFiles) == nqueue:

        cmd = 'hadd -f %s.root ' % d
        cmd += ' '.join(rootFiles)

        print
        print "*"*40
        print "*", d, "is complete. Merging and clean up"
        print "*"*40
        print
        os.system(cmd)
        os.system("rm -rf %s" % d)
    else:
        print
        print "+"*40
        print "+", d, "is incomplete. wait for the job completion"
        print "+"*40
        print

pass1Dir = "pass1"
outFiles = []
for sample in os.listdir(pass1Dir):
    sample = pathjoin(pass1Dir, sample)
    if not isdir(sample): continue

    hadd(pathjoin(sample,'central'))
    outFiles.append(pathjoin(sample,'central.root'))

    for category in os.listdir(sample):
        if category == 'central': continue
        category = pathjoin(sample, category)
        if not isdir(category): continue

        for direction in os.listdir(category):
            direction = pathjoin(category, direction)
            if not isdir(direction): continue
            hadd(direction)
            outFiles.append(direction+'.root')

for l in outFiles:
    print l
