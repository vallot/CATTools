import FWCore.ParameterSet.Config as cms
import os

process = cms.Process("Ana")

process.load("FWCore.MessageLogger.MessageLogger_cfi")
#process.load("Configuration.StandardSequences.Services_cff")
#process.load("SimGeneral.HepPDTESSource.pythiapdt_cfi")

process.options = cms.untracked.PSet(
#    wantSummary = cms.untracked.bool(True),
    allowUnscheduled = cms.untracked.bool(True),
)
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )
process.MessageLogger.cerr.FwkReport.reportEvery = 10000

process.source = cms.Source("PoolSource", fileNames = cms.untracked.vstring())
process.source.fileNames = [
'root://cms-xrdr.sdfarm.kr:1094//xrd/store/group/CAT/TT_TuneCUETP8M1_13TeV-powheg-pythia8/v7-4-0_RunIISpring15DR74-Asympt25ns_MCRUN2_74_V9-v2/150909_163325/0000/catTuple_2.root'
#'file:///store1/jhgoh/CAT/catTuple__TT_TuneCUETP8M1_13TeV-powheg-pythia8__V7-3-6.root',
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
#    solver = cms.string("MAOS"),
#    solver = cms.string("DESYSmeared"),
#    solver = cms.string("DESYMassLoop"),
    muons = cms.InputTag("catMuons"),
    electrons = cms.InputTag("catElectrons"),
    jets = cms.InputTag("catJets"),
    mets = cms.InputTag("catMETs"),
)

process.filterRECO = cms.EDProducer("CATTriggerPacker",
    src = cms.InputTag("catTrigger"),
    triggersToMatch = cms.vstring(
        "CSCTightHaloFilter",
        "EcalDeadCellTriggerPrimitiveFilter",
        "HBHENoiseFilter",
        "eeBadScFilter",
        "goodVertices",
    ),
)

process.HLTMu = cms.EDProducer("CATTriggerPacker",
    src = cms.InputTag("catTrigger"),
    triggersToMatch = cms.vstring(
    ),
)

process.HLTEl = cms.EDProducer("CATTriggerPacker",
    src = cms.InputTag("catTrigger"),
    triggersToMatch = cms.vstring(
        "HLT_Ele12_CaloIdL_TrackIdL_IsoVL",
        "HLT_Ele17_CaloIdL_TrackIdL_IsoVL",
        "HLT_Ele27eta2p1_WPLooseGsf_TriCentralPFJet30",
    ),
)

process.HLTMuMu = cms.EDProducer("CATTriggerPacker",
    src = cms.InputTag("catTrigger"),
    triggersToMatch = cms.vstring(
        "HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ",
        "HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ",
        "HLT_Mu17_Mu8_DZ",
        "HLT_Mu17_TkMu8_DZ",
    ),
)

process.HLTElEl = cms.EDProducer("CATTriggerPacker",
    src = cms.InputTag("catTrigger"),
    triggersToMatch = cms.vstring(
        "HLT_DoubleEle33_CaloIdL_GsfTrkIdVL",
        "HLT_Ele17_Ele12_CaloIdL_TrackIdL_IsoVL_DZ",
        "HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL",
        "HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ",
    ),
)


process.HLTMuEl = cms.EDProducer("CATTriggerPacker",
    src = cms.InputTag("catTrigger"),
    triggersToMatch = cms.vstring(
        "HLT_Mu17_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL",
        "HLT_Mu8_TrkIsoVVL_Ele17_CaloIdL_TrackIdL_IsoVL",
    ),
)

process.ntuple = cms.EDAnalyzer("GenericNtupleMaker",
    failureMode = cms.untracked.string("keep"), # choose one among keep/skip/error
    #failureMode = cms.untracked.string("error"), # choose one among keep/skip/error
    eventCounters = cms.vstring("nEventsTotal"), #"nEventsTotal", "nEventsClean", "nEventsPAT"),
    int = cms.PSet(
        nVertex   = cms.InputTag("catVertex:nGoodPV"),
        filterRECO = cms.InputTag("filterRECO:and"),
        HLTMuMu = cms.InputTag("HLTMuMu:or"),
        HLTElEl = cms.InputTag("HLTElEl:or"),
        HLTMuEl = cms.InputTag("HLTMuEl:or"),
        #HLTMu = cms.InputTag("recoEventInfo","HLTSingleMu"),
        HLTEl = cms.InputTag("HLTEl:or"),
    ),
    double = cms.PSet(
        puWeight   = cms.InputTag("pileupWeight"),
        puWeightUp = cms.InputTag("pileupWeight", "up"),
        puWeightDn = cms.InputTag("pileupWeight", "dn"),
    ),
    doubles = cms.PSet(
        pdfWeight = cms.InputTag("pdfWeight"),
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

