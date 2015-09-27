import FWCore.ParameterSet.Config as cms
import os

process = cms.Process("Ana")

process.load("FWCore.MessageLogger.MessageLogger_cfi")
process.load("Configuration.StandardSequences.Services_cff")
#process.load("SimGeneral.HepPDTESSource.pythiapdt_cfi")

process.options = cms.untracked.PSet(
#    wantSummary = cms.untracked.bool(True),
    allowUnscheduled = cms.untracked.bool(True),
)
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )
process.MessageLogger.cerr.FwkReport.reportEvery = 10000

process.source = cms.Source("PoolSource", fileNames = cms.untracked.vstring())
process.source.fileNames = [
#'root://cms-xrdr.sdfarm.kr:1094//xrd/store/group/CAT/TT_TuneCUETP8M1_13TeV-powheg-pythia8/v7-4-0_RunIISpring15DR74-Asympt25ns_MCRUN2_74_V9-v2/150909_163325/0000/catTuple_2.root'
#'file:///store1/jhgoh/CAT/catTuple__TT_TuneCUETP8M1_13TeV-powheg-pythia8__V7-3-6.root',
'file:../../CatProducer/prod/catTuple.root',
]

process.out = cms.OutputModule("PoolOutputModule",
    fileName = cms.untracked.string("out.root"),
    outputCommands = cms.untracked.vstring(
        "drop *",
#        "keep *_*_*_Ana",
        "keep *_ttbarEvent_*_*",
    )
)

#process.outPath = cms.EndPath(process.out)

process.ttbarEvent = cms.EDProducer("TTbarDileptonEventSelector",
    muons = cms.InputTag("catMuons"),
    electrons = cms.InputTag("catElectrons"),
    jets = cms.InputTag("catJets"),
    mets = cms.InputTag("catMETs"),
    keepVetoLeptons = cms.bool(False),
    checkOverlapFromVetoLepton = cms.bool(False),
    sortByBtag = cms.bool(False),
    eleIdName = cms.string("cutBasedElectronID-Spring15-25ns-V1-standalone-medium"),
    eleVetoIdName = cms.string("cutBasedElectronID-Spring15-25ns-V1-standalone-veto"),
    bTagName = cms.string("combinedInclusiveSecondaryVertexV2BJetTags"),
)

process.load("CATTools.CatAnalyzer.ttbarDileptonKinSolutionProducer_cfi")
process.ttbar = process.ttbarDileptonKin.clone(
    leptons = cms.InputTag("ttbarEvent:leptons"),
    jets = cms.InputTag("ttbarEvent:jets"),
    met = cms.InputTag("ttbarEvent:met"),
    metphi = cms.InputTag("ttbarEvent:metphi"),
)

process.filterRECO = cms.EDProducer("CATTriggerPacker",
    triggerResults = cms.InputTag("TriggerResults::CAT"),
    triggerPrescales = cms.InputTag("patTrigger"),
    combineByOr = cms.bool(False),
    triggersToMatch = cms.vstring(
        "CSCTightHaloFilter",
        "EcalDeadCellTriggerPrimitiveFilter",
        "HBHENoiseFilter",
        "eeBadScFilter",
        "goodVertices",
    ),
)

process.HLTMu = cms.EDProducer("CATTriggerPacker",
    triggerResults = cms.InputTag("TriggerResults::HLT"),
    triggerPrescales = cms.InputTag("patTrigger"),
    combineByOr = cms.bool(True),
    triggersToMatch = cms.vstring(
    ),
)

process.HLTEl = cms.EDProducer("CATTriggerPacker",
    triggerResults = cms.InputTag("TriggerResults::HLT"),
    triggerPrescales = cms.InputTag("patTrigger"),
    combineByOr = cms.bool(True),
    triggersToMatch = cms.vstring(
        "HLT_Ele12_CaloIdL_TrackIdL_IsoVL",
        "HLT_Ele17_CaloIdL_TrackIdL_IsoVL",
        "HLT_Ele27eta2p1_WPLooseGsf_TriCentralPFJet30",
    ),
)

process.HLTMuMu = cms.EDProducer("CATTriggerPacker",
    triggerResults = cms.InputTag("TriggerResults::HLT"),
    triggerPrescales = cms.InputTag("patTrigger"),
    combineByOr = cms.bool(True),
    triggersToMatch = cms.vstring(
        "HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ",
        "HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ",
        "HLT_Mu17_Mu8_DZ",
        "HLT_Mu17_TkMu8_DZ",
    ),
)

process.HLTElEl = cms.EDProducer("CATTriggerPacker",
    triggerResults = cms.InputTag("TriggerResults::HLT"),
    triggerPrescales = cms.InputTag("patTrigger"),
    combineByOr = cms.bool(True),
    triggersToMatch = cms.vstring(
        "HLT_DoubleEle33_CaloIdL_GsfTrkIdVL",
        "HLT_Ele17_Ele12_CaloIdL_TrackIdL_IsoVL_DZ",
        "HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL",
        "HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ",
    ),
)


process.HLTMuEl = cms.EDProducer("CATTriggerPacker",
    triggerResults = cms.InputTag("TriggerResults::HLT"),
    triggerPrescales = cms.InputTag("patTrigger"),
    combineByOr = cms.bool(True),
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
        filterRECO = cms.InputTag("filterRECO"),
        HLTMuMu = cms.InputTag("HLTMuMu"),
        HLTElEl = cms.InputTag("HLTElEl"),
        HLTMuEl = cms.InputTag("HLTMuEl"),
        #HLTMu = cms.InputTag("recoEventInfo","HLTSingleMu"),
        HLTEl = cms.InputTag("HLTEl"),
        nBjets = cms.InputTag("ttbarEvent:nBjets"),
    ),
    float = cms.PSet(
        puWeight   = cms.InputTag("pileupWeight"),
        #puWeightUp = cms.InputTag("pileupWeight", "up"),
        #puWeightDn = cms.InputTag("pileupWeight", "dn"),
        met = cms.InputTag("ttbarEvent:met"),
        metphi = cms.InputTag("ttbarEvent:metphi"),
    ),
    floats = cms.PSet(
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
        ),
    ),
)

process.load("CATTools.CatAnalyzer.ttbarDileptonKinSolutionAlgos_cff")
for algo in ["CMSKin", "MT2", "DESYMassLoop", "DESYSmeared"]:
    setattr(process, 'ttbar'+algo, process.ttbar.clone(solver = getattr(process, 'ttbarDileptonKinAlgoPSet'+algo)))
    setattr(process.ntuple.floats, 'ttbar%s_mLL' % algo, cms.InputTag("ttbar%s:mLL" % algo))
    setattr(process.ntuple.floats, 'ttbar%s_mLB' % algo, cms.InputTag("ttbar%s:mLB" % algo))
    setattr(process.ntuple.floats, 'ttbar%s_dphi' % algo, cms.InputTag("ttbar%s:dphi" % algo))
    setattr(process.ntuple.cands, 'ttbar'+algo, cms.PSet(
        src = cms.InputTag('ttbar'+algo),
        exprs = cms.untracked.PSet(
            pt  = cms.string("pt"),
            eta = cms.string("eta"),
            phi = cms.string("phi"),
            m   = cms.string("mass"),
            pdgId = cms.string("pdgId"),
            q = cms.string("charge"),
        )
    ))
    if 'DESYSmeared' in algo:
        process.RandomNumberGeneratorService.ttbarDESYSmeared = cms.PSet(
            initialSeed = cms.untracked.uint32(123456),
            engineName = cms.untracked.string('TRandom3')
        )

process.load("CATTools.CatProducer.mcTruthTop.mcTruthTop_cff")
delattr(process, 'pseudoTop')
process.p = cms.Path(
    process.ntuple
)

process.TFileService = cms.Service("TFileService",
    fileName = cms.string("ntuple.root"),
)

