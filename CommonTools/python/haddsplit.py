#!/usr/bin/env python

import os

def haddsplit(srcDir, inputFiles, destDir):
    if os.path.exists(destDir) and len(os.listdir(destDir)) > 1:
        print "Directory", destDir, "is not empty. skip."
        return False
    if not os.path.exists(destDir): os.makedirs(destDir)

    ## estimate total size and merge not to exceed 1Gbytes
    totalSize = sum([os.stat(f).st_size for f in inputFiles])
    fSize = 4*1024.*1024*1024 # 4GB
    nOutput = int(totalSize/fSize)+1

    nTotalFiles = len(inputFiles)
    nFiles = nTotalFiles/nOutput

    prefix = '.'.join(os.path.basename(inputFiles[0]).split('.')[:-1])
    if '_' in prefix and prefix.split('_')[-1].isdigit():
        prefix = '_'.join(prefix.split('_')[:-1])

    for i in range(nOutput):
        cmd = 'hadd -f %s/%s_%d.root ' % (destDir, prefix, i)
        if i == nOutput-1: cmd += ' '.join(inputFiles[i*nFiles:])
        else: cmd += ' '.join(inputFiles[i*nFiles:(i+1)*nFiles])
        os.system(cmd)
    for x in inputFiles: os.system("rm -f %s" % (x))

    return True

