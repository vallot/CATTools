import FWCore.ParameterSet.Config as cms
from CATTools.CatAnalyzer.leptonSF_cff import *
from CATTools.CatAnalyzer.flatGenWeights_cfi import *

eventsTTLL = cms.EDFilter("TTLLEventSelector",
    isMC = cms.bool(True),
    ## alwaysAcceptAfter : Accept event even though selection may fail _AFTER_ this step
    ## Use case: store ntuple only for events that passes step4
    applyFilterAt = cms.int32(4),

    # Physics objects
    muon = cms.PSet(
        src = cms.InputTag("catMuons"),
        scaleDirection = cms.int32(0),
        #scaleDirection = cms.int32(-1),
        #scaleDirection = cms.int32(+1),
        efficiencySF = muonSFTight,
        efficiencySFDirection = cms.int32(0),
    ),

    electron = cms.PSet(
        src = cms.InputTag("catElectrons"),
        idName = cms.string("cutBasedElectronID-Spring15-25ns-V1-standalone-medium"),
        #idName = cms.string("mvaEleID-Spring15-25ns-Trig-V1-wp90"),
        scaleDirection = cms.int32(0),
        #scaleDirection = cms.int32(-1),
        #scaleDirection = cms.int32(+1),
        efficiencySF = electronSFCutBasedIDMediumWP,
        efficiencySFDirection = cms.int32(0),
    ),

    jet = cms.PSet(
        src = cms.InputTag("catJets"),
        bTagName = cms.string("pfCombinedInclusiveSecondaryVertexV2BJetTags"),
        bTagWP = cms.string("CSVL"),
        scaleDirection = cms.int32(0),
        #scaleDirection = cms.int32(-1),
        #scaleDirection = cms.int32(+1),
        resolDirection = cms.int32(0),
        #resolDirection = cms.int32(-1),
        #resolDirection = cms.int32(+1),
        skipJER = cms.bool(False), # Needed for synchronization
    ),

    met = cms.PSet(
        src = cms.InputTag("catMETs"),
    ),

    vertex = cms.PSet(
        nVertex = cms.InputTag("catVertex", "nGoodPV"),
        #src = cms.InputTag("catVertex"),
        pileupWeight = cms.InputTag("pileupWeight"),
    ),

    # Filters
    filters = cms.PSet(
        filterRECO = cms.InputTag("filterRECO"),
        trigMUEL = cms.InputTag("filterTrigMUEL"),
        trigMUMU = cms.InputTag("filterTrigMUMU"),
        trigELEL = cms.InputTag("filterTrigELEL"),
        ignoreTrig = cms.bool(False), # Accept event even if it does not pass HLT. Needed for synchronization
    ),

    # Event weights
    genWeight = cms.PSet(
        index = cms.int32(-1),
        src = cms.InputTag("flatGenWeights"),
    ),
    extWeights = cms.VInputTag(),
)

