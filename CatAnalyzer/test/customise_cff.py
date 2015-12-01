import FWCore.ParameterSet.Config as cms
import os

def customise(process):

    lumiFile = 'Cert_246908-260627_13TeV_PromptReco_Collisions15_25ns_JSON.txt'

    runOnMC = True
    for i in process.source.fileNames:
        if 'Run2015' in i:
            runOnMC=False
   
    isTTbar=False
    for i in process.source.fileNames:
        if '/TT' in i or '/tt' in i:
            isTTbar=True
   
    if not runOnMC:
        from FWCore.PythonUtilities.LumiList import LumiList
        lumiList = LumiList(os.environ["CMSSW_BASE"]+'/src/CATTools/CatProducer/prod/LumiMask/'+lumiFile)    
        #lumiList = LumiList(os.environ["CMSSW_BASE"]+'/src/CATTools/CommonTools/test/ttbb/'+lumiFile)    
        process.source.lumisToProcess = lumiList.getVLuminosityBlockRange()  
   
   
#########################################
        
#customise(process)
