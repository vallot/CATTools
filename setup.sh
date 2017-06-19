#!/bin/bash
yum -y install glibc-devel
source /cvmfs/cms.cern.ch/cmsset_default.sh
scramv1 project CMSSW CMSSW_8_0_26_patch1
cd CMSSW_8_0_26_patch1/src
eval `scramv1 runtime -sh`
gcc -v
ln -s /CATTools .
cd CATTools

if [ -n $PR_SLUG -a -n $PR_BRANCH ]; then
    GITHUBNAME=`echo $PR_SLUG | cut -f1 -d/`
    #GITHUBREPO=`echo $PULL_REQUEST_SLUG | cut -f2- -d/`

    git remote add $GITHUBNAME https://github.com/$PR_SLUG
    git fetch $GITHUBNAME
    git checkout $PR_BRANCH
fi

git submodule init
git submodule update
cd CommonTools/scripts
git checkout master
git pull
cd ../../..

scram b -j2
