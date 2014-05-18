#!/bin/bash

echo "SETUP ENVIRONMENT "
export DCACHE_RA_BUFFER=50000000
export ROOTSYS='/swmgrs/cmss/slc4_ia32_gcc345/lcg/root/5.18.00a-cms17'
export DYLD_LIBRARY_PATH=$ROOTSYS/lib
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$ROOTSYS/lib:`pwd`:/usr/local/lib
echo ROOTSYS = $ROOTSYS
echo "DONE"
chmod +x EXECUTABLE
echo "EXECUTE CODE"
./EXECUTABLE
echo "DONE"
