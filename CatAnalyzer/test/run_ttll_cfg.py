import FWCore.ParameterSet.Config as cms
import os

process = cms.Process("Ana")

process.load("FWCore.MessageLogger.MessageLogger_cfi")
#process.load("Configuration.StandardSequences.Services_cff")
#process.load("SimGeneral.HepPDTESSource.pythiapdt_cfi")

process.options = cms.untracked.PSet(
    wantSummary = cms.untracked.bool(True),
    allowUnscheduled = cms.untracked.bool(True),
)
#process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(100) )
process.MessageLogger.cerr.FwkReport.reportEvery = 10000

process.source = cms.Source("PoolSource", fileNames = cms.untracked.vstring())
process.source.fileNames = [
'root://cms-xrdr.sdfarm.kr:1094//xrd/store/group/CAT/TT_TuneCUETP8M1_13TeV-powheg-pythia8/v7-3-4_RunIISpring15DR74-Asympt50ns_MCRUN2_74_V9A-v4/150810_215031/0000/catTuple_1.root',
]

process.out = cms.OutputModule("PoolOutputModule",
    fileName = cms.untracked.string("out.root"),
    outputCommands = cms.untracked.vstring(
        "drop *",
        "keep *_*_*_Ana",
    )
)

#process.outPath = cms.EndPath(process.out)

process.ttbar = cms.EDProducer("TTbarDileptonProducer",
#    solver = cms.string("Default"),
    solver = cms.string("CMSKIN"),
#    solver = cms.string("NUWGT"),
#    solver = cms.string("MT2"),
    muons = cms.InputTag("catMuons"),
    electrons = cms.InputTag("catElectrons"),
    jets = cms.InputTag("catJets"),
    mets = cms.InputTag("catMETs"),
)

process.filterRECO = cms.EDProducer("CATTriggerPacker",
    srcs = cms.VInputTag(
        cms.InputTag("catTrigger", "CSCTightHaloFilter"),
        cms.InputTag("catTrigger", "EcalDeadCellTriggerPrimitiveFilter"),
        cms.InputTag("catTrigger", "HBHENoiseFilter"),
        cms.InputTag("catTrigger", "eeBadScFilter"),
        cms.InputTag("catTrigger", "goodVertices"),
    ),
)

process.HLTMu = cms.EDProducer("CATTriggerPacker",
    srcs = cms.VInputTag(
        cms.InputTag("catTrigger", ""),
    ),
)

process.HLTEl = cms.EDProducer("CATTriggerPacker",
    srcs = cms.VInputTag(
        cms.InputTag("catTrigger", "HLTEle12CaloIdLTrackIdLIsoVL"),
        cms.InputTag("catTrigger", "HLTEle17CaloIdLTrackIdLIsoVL"),
        cms.InputTag("catTrigger", "HLTEle27eta2p1WPLooseGsfTriCentralPFJet30"),
    ),
)

process.HLTMuMu = cms.EDProducer("CATTriggerPacker",
    srcs = cms.VInputTag(
        cms.InputTag("catTrigger", "HLTMu17Mu8DZ"),
        cms.InputTag("catTrigger", "HLTMu17TkMu8DZ"),
        cms.InputTag("catTrigger", "HLTMu17TrkIsoVVLMu8TrkIsoVVL"),
        cms.InputTag("catTrigger", "HLTMu17TrkIsoVVLMu8TrkIsoVVLDZ"),
        cms.InputTag("catTrigger", "HLTMu17TrkIsoVVLTkMu8TrkIsoVVL"),
    ),
)

process.HLTElEl = cms.EDProducer("CATTriggerPacker",
    srcs = cms.VInputTag(
        cms.InputTag("catTrigger", "HLTDoubleEle33CaloIdLGsfTrkIdVL"),
        cms.InputTag("catTrigger", "HLTEle17Ele12CaloIdLTrackIdLIsoVLDZ"),
        cms.InputTag("catTrigger", "HLTEle23Ele12CaloIdLTrackIdLIsoVL"),
        cms.InputTag("catTrigger", "HLTEle23Ele12CaloIdLTrackIdLIsoVLDZ"),
    ),
)

process.HLTMuEl = cms.EDProducer("CATTriggerPacker",
    srcs = cms.VInputTag(
        cms.InputTag("catTrigger", "HLTMu17TrkIsoVVLEle12CaloIdLTrackIdLIsoVL"),
        cms.InputTag("catTrigger", "HLTMu8TrkIsoVVLEle17CaloIdLTrackIdLIsoVL"),
    ),
)

process.ntuple = cms.EDAnalyzer("GenericNtupleMaker",
    failureMode = cms.untracked.string("keep"), # choose one among keep/skip/error
    #failureMode = cms.untracked.string("error"), # choose one among keep/skip/error
    eventCounters = cms.vstring("nEventsTotal"), #"nEventsTotal", "nEventsClean", "nEventsPAT"),
    int = cms.PSet(
        nVertex   = cms.PSet(src = cms.InputTag("catVertex:nGoodPV")),
        filterRECO = cms.PSet(src = cms.InputTag("filterRECO:and")),
        HLTMuMu = cms.PSet(src = cms.InputTag("HLTMuMu:or")),
        HLTElEl = cms.PSet(src = cms.InputTag("HLTElEl:or")),
        HLTMuEl = cms.PSet(src = cms.InputTag("HLTMuEl:or")),
        #HLTMu = cms.PSet(src = cms.InputTag("recoEventInfo","HLTSingleMu")),
        HLTEl = cms.PSet(src = cms.InputTag("HLTEl:or")),
    ),
    double = cms.PSet(
        puWeight   = cms.PSet(src = cms.InputTag("pileupWeight")),
        puWeightUp = cms.PSet(src = cms.InputTag("pileupWeight", "up")),
        puWeightDn = cms.PSet(src = cms.InputTag("pileupWeight", "dn")),
    ),
    doubles = cms.PSet(
        pdfWeight = cms.PSet(src = cms.InputTag("pdfWeight")),
    ),
    cands = cms.PSet(
        pseudoTop = cms.PSet(
            src = cms.InputTag("pseudoTop"),
            exprs = cms.untracked.PSet(
                pt  = cms.string("pt"),
                eta = cms.string("eta"),
                phi = cms.string("phi"),
                m   = cms.string("mass"),
                pdgId = cms.string("pdgId"),
                q = cms.string("charge"),
                #status = cms.string("status"),
            ),
            selections = cms.untracked.PSet(),
        ),
        partonTop = cms.PSet(
            src = cms.InputTag("partonTop"),
            #index = cms.untracked.int32(0),
            exprs = cms.untracked.PSet(
                pt  = cms.string("pt"),
                eta = cms.string("eta"),
                phi = cms.string("phi"),
                m   = cms.string("mass"),
                pdgId = cms.string("pdgId"),
                q = cms.string("charge"),
                #status = cms.string("status"),
            ),
            selections = cms.untracked.PSet(),
        ),
        ttbarLep = cms.PSet(
            src = cms.InputTag("ttbar"),
            exprs = cms.untracked.PSet(
                pt  = cms.string("pt"),
                eta = cms.string("eta"),
                phi = cms.string("phi"),
                m   = cms.string("mass"),
                pdgId = cms.string("pdgId"),
                q = cms.string("charge"),
            ),
            selections = cms.untracked.PSet(),
        ),
    ),
)

process.load("CATTools.CatProducer.pseudoTop_cff")
delattr(process, 'pseudoTop')
process.p = cms.Path(
    process.filterRECO*
    process.ntuple
)

process.TFileService = cms.Service("TFileService",
    fileName = cms.string("ntuple.root"),
)

