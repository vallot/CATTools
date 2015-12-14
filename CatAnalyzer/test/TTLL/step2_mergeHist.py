#!/usr/bin/env python

import sys, os
from os.path import isdir as isdir
from os.path import join as pathjoin

def hadd(d):
    d = d.rstrip('/')
    if not isdir(d): return 
    rootFiles = [pathjoin(d, x) for x in os.listdir(d) if x.endswith('.root')]

    cmd = 'hadd -f %s.root ' % d
    cmd += ' '.join(rootFiles)

    os.system(cmd)

step1Dir = "step1_makeHist"
outFiles = []
for sample in os.listdir(step1Dir):
    sample = pathjoin(step1Dir, sample)
    if not isdir(sample): continue

    hadd(pathjoin(sample,'central'))
    outFiles.append(pathjoin(sample,'central.root'))

    for category in os.listdir(sample):
        if category == 'central': continue
        category = pathjoin(sample, category)
        hadd(category)
        outFiles.append(category+'.root')

for l in outFiles:
    print l
