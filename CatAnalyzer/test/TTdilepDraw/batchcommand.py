import json
import os

jsonfile = os.environ["CMSSW_BASE"]+"/src/CATTools/CatAnalyzer/data/dataset.json"
with open(jsonfile) as data_file:    
    data = json.load(data_file)

mcfilelist = ['TT_powheg',
              'WJets',
              'SingleTbar_tW',
              'SingleTop_tW',
              'ZZ',
              'WW',
              'WZ',
              'DYJets']
rdfilelist = ['MuonEG']

for samplename in mcfilelist:
	fileList = os.environ["CMSSW_BASE"]+"/src/CATTools/CatAnalyzer/data/dataset_"+samplename+".txt"
	print "create-batch --jobName "+samplename+" --fileList "+fileList+" --maxFiles 5 --cfg run_TtbarDiLeptonAnalyzer.py" 

print
for samplename in rdfilelist:
	fileList = os.environ["CMSSW_BASE"]+"/src/CATTools/CatAnalyzer/data/dataset_"+samplename+".txt"
	print "create-batch --jobName "+samplename+" --fileList "+fileList+" --maxFiles 5 --cfg run_TtbarDiLeptonAnalyzer.py" 

print "condor_q"
