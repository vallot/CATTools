import FWCore.ParameterSet.Config as cms
process = cms.Process("CATeX")

process.load("FWCore.MessageService.MessageLogger_cfi")
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )
process.options = cms.untracked.PSet( wantSummary = cms.untracked.bool(False) )

process.options.allowUnscheduled = cms.untracked.bool(True)
process.MessageLogger.cerr.FwkReport.reportEvery = 50000

process.source = cms.Source("PoolSource", fileNames = cms.untracked.vstring())
process.source.fileNames = [
#    'file:/state/partition1/store/user/jhgoh/FCNC/Synchronization/201611/v802/TT_TuneCUETP8M1_13TeV-powheg-pythia8.root',
#    'file:catTuple_seed81.root',
    'file:catTuple.root',
]

process.load("CATTools.CatAnalyzer.filters_cff")
#process.load("CATTools.Validation.ttllEventSelector_cfi")
#process.load("CATTools.Validation.validation_cff")
process.load("CATTools.Validation.eventsTopFCNC_cff")
process.filterTrigMU.triggersToMatch = ['HLT_IsoMu24_v', 'HLT_IsoTkMu24_v',]
process.filterTrigEL.triggersToMatch = ['HLT_Ele32_eta2p1_WPTight_Gsf_v']
process.eventsTopFCNC.electron.idName = "cutBasedElectronID-Spring15-25ns-V1-standalone-medium"
process.eventsTopFCNC.electron.vetoIdName = "cutBasedElectronID-Spring15-25ns-V1-standalone-loose"
#process.eventsTopFCNC.vertex.src = "offlineSlimmedPrimaryVertices"
process.el = process.eventsTopFCNC.clone(
    channel = cms.string("electron"),
    eventFile = cms.untracked.string("eventList_electron.txt"))
process.mu = process.eventsTopFCNC.clone(
    channel = cms.string("muon"),
    eventFile = cms.untracked.string("eventList_muon.txt"))
delattr(process, 'eventsTopFCNC')

process.TFileService = cms.Service("TFileService",
    fileName = cms.string("hist.root"),
)

process.p_el = cms.Path(process.el)
process.p_mu = cms.Path(process.mu)

