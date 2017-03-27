#!/usr/bin/env python

import os

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

