#!/bin/bash

if [ _$1 == _ ]; then
    echo Usage: $0 CATToolsVersion
    echo example: $0 v7-6-6
    exit 1
fi

VERSION=$1

#dumpRoot pass2/preview.root
cp pass2/plots.json preview/
cp pass2/cutflow.json preview/
ssh -p8022 higgs.hanyang.ac.kr "mkdir /var/www/html/CATTools/preview/$VERSION"
scp -rP8022 preview/* higgs.hanyang.ac.kr:/var/www/html/CATTools/preview/$VERSION
