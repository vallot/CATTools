import FWCore.ParameterSet.Config as cms

ntuple = cms.EDAnalyzer("GenericNtupleMaker",
    failureMode = cms.untracked.string("error"),
    eventCounters = cms.vstring(),
    int = cms.PSet(
        nPV = cms.InputTag("catVertex", "nGoodPV"),
    ),
    ints = cms.PSet(
        modes = cms.InputTag("partonTop", "modes"),
    ),
    float = cms.PSet(
        weight = cms.InputTag("eventsTTLL", "weight"),
        met_pt = cms.InputTag("eventsTTLL", "met"),
        met_phi = cms.InputTag("eventsTTLL", "metphi"),
    ),
    floats = cms.PSet(
    ),
    cands = cms.PSet(
        leptons = cms.PSet(
            src = cms.InputTag("eventsTTLL", "leptons"),
            exprs = cms.untracked.PSet(
                pt  = cms.string("pt"),
                eta = cms.string("eta"),
                phi = cms.string("phi"),
                m   = cms.string("mass"),
                q = cms.string("charge"),
            ),
        ),
        jets = cms.PSet(
            src = cms.InputTag("eventsTTLL", "jets"),
            exprs = cms.untracked.PSet(
                pt  = cms.string("pt"),
                eta = cms.string("eta"),
                phi = cms.string("phi"),
                m   = cms.string("mass"),
                bTag = cms.string("bDiscriminator('pfCombinedInclusiveSecondaryVertexV2BJetTags')"),
                q = cms.string("charge"),
            ),
        ),
    ),
)

