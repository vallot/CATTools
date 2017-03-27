#!/bin/bash

if [ _$1 == _ ]; then
    echo Usage: $0 CATToolsVersion
    echo example: $0 v7-6-6
    exit 1
fi

VERSION=$1

#dumpRoot pass2/preview.root
tar czf pass2.tgz pass2
cp pass2/plots.json preview/
cp pass2/cutflow.json preview/
ssh gate.sscc.uos.ac.kr "mkdir /var/www/html/CATTools/preview/$VERSION"
scp -r preview/* gate.sscc.uos.ac.kr:/var/www/html/CATTools/preview/$VERSION
scp pass2.tgz gate.sscc.uos.ac.kr:/var/www/html/CATTools/preview/$VERSION/
