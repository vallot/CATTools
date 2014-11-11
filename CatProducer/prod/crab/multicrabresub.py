#!/usr/bin/python

import os,time

def runCommand(commandstring, loopn):
    print "loop", loopn, ",", commandstring
    os.system(commandstring)
    
multicrabdir = []
for fn in os.listdir('.'):
    if os.path.isdir(fn):
        if fn.startswith('multicrab_'):
            multicrabdir.append(fn)
            print "found multicrab job dir:", fn

for n in range(50):
    for d in multicrabdir:
        runCommand("multicrab -submit 500 -c "+d, n)
        runCommand("multicrab -status -c "+d, n)
        runCommand("multicrab -get -c "+d, n)
        runCommand("multicrab -resubmit bad -c "+d, n)
