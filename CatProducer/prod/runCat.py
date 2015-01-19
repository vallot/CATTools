from CATTools.CatProducer.catTemplate_cfg import *
## some options
doSecVertex=True # for jpsi candidates
    
## setting up arguements
from FWCore.ParameterSet.VarParsing import VarParsing
options = VarParsing ('python')
options.register('runOnMC', True, VarParsing.multiplicity.singleton, VarParsing.varType.bool, "runOnMC: 1  default")
options.register('useMiniAOD', True, VarParsing.multiplicity.singleton, VarParsing.varType.bool, "useMiniAOD: 1  default")
options.register('globalTag', '', VarParsing.multiplicity.singleton, VarParsing.varType.string, "globalTag: 1  default")

options.parseArguments()
runOnMC = options.runOnMC
useMiniAOD = options.useMiniAOD
globalTag = options.globalTag

print "runOnMC =",runOnMC,"and useMiniAOD =",useMiniAOD

from Configuration.AlCa.GlobalTag import GlobalTag
if useMiniAOD:
    if runOnMC:
        process.GlobalTag = GlobalTag(process.GlobalTag, 'PLS170_V7AN2::All', '')
    else:
        process.GlobalTag = GlobalTag(process.GlobalTag, 'GR_70_V2_AN1::All', '')

####################################################################################################
## from miniAOD/patTuple_mini.py to run miniAOD maker when starting from AOD
if not useMiniAOD:
    if not globalTag:
        print "ERROR!!!! Need correct globalTag to run on AOD"
    process.GlobalTag = GlobalTag(process.GlobalTag, globalTag, '')
    process.load('Configuration.StandardSequences.PAT_cff')
    process.load("RecoVertex.AdaptiveVertexFinder.inclusiveVertexing_cff")
    from PhysicsTools.PatAlgos.slimming.miniAOD_tools import *
    if runOnMC:
        miniAOD_customizeAllMC(process)
    else :
        miniAOD_customizeAllData(process)
####################################################################################################
## setting up catTools
print "process.GlobalTag.globaltag =",process.GlobalTag.globaltag
process.load('PhysicsTools.PatAlgos.slimming.unpackedTracksAndVertices_cfi')

from CATTools.CatProducer.catSetup_cff import *
catSetup(process, runOnMC, doSecVertex)

process.out.outputCommands += [
        'keep *_slimmedPhotons*_*_*',
        'keep *_slimmedElectrons_*_*',
        'keep *_slimmedMuons*_*_*',
        'keep *_slimmedTaus*_*_*',
        'keep *_slimmedJets*_*_*',
        'keep *_slimmedMETs*_*_*',
        'keep *_slimmedSecondaryVertices*_*_*',
        ## add extra METs

        'keep recoPhotonCores_reducedEgamma_*_*',
        'keep recoGsfElectronCores_reducedEgamma_*_*',
        'keep recoConversions_reducedEgamma_*_*',
        'keep recoSuperClusters_reducedEgamma_*_*',
        'keep recoCaloClusters_reducedEgamma_*_*',
        'keep EcalRecHitsSorted_reducedEgamma_*_*',
        

        'drop *_*_caloTowers_*',
        'drop *_*_pfCandidates_*',
        'drop *_*_genJets_*',

        'keep *_offlineBeamSpot_*_*',
        'keep *_offlineSlimmedPrimaryVertices_*_*',
        'keep patPackedCandidates_packedPFCandidates_*_*',

        'keep double_fixedGridRho*__*', 

        'keep *_selectedPatTrigger_*_*',
        'keep patPackedTriggerPrescales_patTrigger__*',
        'keep *_l1extraParticles_*_*',
        'keep L1GlobalTriggerReadoutRecord_gtDigis_*_*',
        'keep *_TriggerResults_*_HLT',
        'keep *_TriggerResults_*_PAT', # for MET filters
	'keep patPackedCandidates_lostTracks_*_*',
        'keep HcalNoiseSummary_hcalnoise__*',
        'keep *_caTopTagInfos_*_*',
        'keep *_slimmedGenJets_*_*',
        'keep patPackedGenParticles_packedGenParticles_*_*',
        'keep recoGenParticles_prunedGenParticles_*_*',
        'keep LHEEventProduct_*_*_*',
        'keep PileupSummaryInfos_*_*_*',
        'keep GenFilterInfo_*_*_*',
        'keep GenEventInfoProduct_generator_*_*',
        # RUN
        'keep LHERunInfoProduct_*_*_*',
        'keep GenRunInfoProduct_*_*_*',
        'keep L1GtTriggerMenuLite_l1GtTriggerMenuLite__*'
]

process.maxEvents.input = options.maxEvents

process.source.fileNames = options.inputFiles

## to suppress the long output at the end of the job
process.MessageLogger.cerr.threshold = ''
if options.maxEvents < 0:
    process.MessageLogger.cerr.FwkReport.reportEvery = 1000
process.options.wantSummary = False
