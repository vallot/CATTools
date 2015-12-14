import FWCore.ParameterSet.Config as cms

ttll = cms.EDFilter("TTLLEventSelector",
    isMC = cms.bool(True),
    filterCutStepBefore = cms.int32(4),

    # Physics objects
    muon = cms.PSet(
        src = cms.InputTag("catMuons"),
        scale = cms.string("0"),
        #scale = cms.string("+"),
        #scale = cms.string("-"),
    ),

    electron = cms.PSet(
        src = cms.InputTag("catElectrons"),
        idName = cms.string("cutBasedElectronID-Spring15-25ns-V1-standalone-medium"),
        scale = cms.string("0"),
        #scale = cms.string("+"),
        #scale = cms.string("-"),
    ),

    jet = cms.PSet(
        src = cms.InputTag("catJets"),
        bTagName = cms.string("pfCombinedInclusiveSecondaryVertexV2BJetTags"),
        bTagMin = cms.double(0.605),
        scale = cms.string("0"),
        #scale = cms.string("+"),
        #scale = cms.string("-"),
    ),

    met = cms.PSet(
        src = cms.InputTag("catMETs"),
        scale = cms.string("0"),
        #scale = cms.string("+"),
        #scale = cms.string("-"),
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
    ),

    # Event weights
    genWeight = cms.PSet(
        index = cms.uint32(0),
        src = cms.InputTag("genWeight", "genWeight"),
    ),
)

