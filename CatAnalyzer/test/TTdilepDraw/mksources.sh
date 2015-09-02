#!/bin/bash

# ///////////////////////////////////////////////////////////////////////////
# ///    Useage: 1. make a txt file of list of all files                  ///
# ///               in certain dir with full dir at ./texts/              ///
# ///            2. print create-batch running command for each,          ///
# ///               so that make it possible to copy-paste directly.      ///
# ///                                                                     ///
# ///    * ONLY FOR CASE YOU HAVE ONE FOLDER UNDER INPUT DIRECTORY.       ///
# ///////////////////////////////////////////////////////////////////////////

wherefilesin="/cms/scratch/CAT/"
batchdest="/cms/scratch/tt8888tt/cat74/src/CATTools/CatAnalyzer/test/results/"

if [ ! -f dirlist.txt ]; then
	lst=`find $wherefilesin -name catTuple_1.root`
	for l in $lst
	do
		echo ${l%catTuple_1.root} >> dirlist.txt
	done
fi

if [ ! -d texts ]; then
	mkdir texts
fi

while read line
do
	samplename=`echo $line | cut -d'/' -f5`
	`/bin/ls -d $line* > texts/$samplename.txt`
	echo "./create-batch --jobName "$samplename "--fileList texts/"$samplename".txt --maxFiles 10 --cfg run_TtbarDiLeptonAnalyzer.py --transferDest "$batchdest$samplename
done < dirlist.txt

