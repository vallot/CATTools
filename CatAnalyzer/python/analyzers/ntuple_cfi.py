import FWCore.ParameterSet.Config as cms

ntuple = cms.EDAnalyzer("GenericNtupleMaker",
    failureMode = cms.untracked.string("error"),
    #failureMode = cms.untracked.string("skip"),
    eventCounters = cms.vstring("nEventsTotal"),
    int = cms.PSet(
        vertex_n = cms.InputTag("catVertex", "nGoodPV"),
        cutstep = cms.InputTag("eventSel", "cutstep"),
    ),
    float = cms.PSet(
        met_pt = cms.InputTag("eventSel", "met"),
        met_phi = cms.InputTag("eventSel", "metphi"),
    ),
    ints = cms.PSet(),
    floats = cms.PSet(),
    cands = cms.PSet(
        electrons = cms.PSet(
            src = cms.InputTag("eventSel", "electrons"),
            exprsF = cms.untracked.PSet(
                pt  = cms.string("pt"),
                eta = cms.string("eta"),
                phi = cms.string("phi"),
                m   = cms.string("mass"),
                scEta = cms.string("scEta"),
                relIso = cms.string("relIso(0.3)"),
            ),
            exprsI = cms.untracked.PSet(
                q = cms.string("charge"),
                pdgId = cms.string("pdgId"),
            ),
        ),
        muons = cms.PSet(
            src = cms.InputTag("eventSel", "muons"),
            exprsF = cms.untracked.PSet(
                pt  = cms.string("pt"),
                eta = cms.string("eta"),
                phi = cms.string("phi"),
                m   = cms.string("mass"),
                relIso = cms.string("relIso(0.4)"),
            ),
            exprsI = cms.untracked.PSet(
                q = cms.string("charge"),
                pdgId = cms.string("pdgId"),
            ),
        ),
        jets = cms.PSet(
            src = cms.InputTag("eventSel", "jets"),
            exprsF = cms.untracked.PSet(
                pt  = cms.string("pt"),
                eta = cms.string("eta"),
                phi = cms.string("phi"),
                m   = cms.string("mass"),
                CSV = cms.string("bDiscriminator('pfCombinedInclusiveSecondaryVertexV2BJetTags')"),
                CvsL = cms.string("bDiscriminator('pfCombinedCvsLJetTags')"),
                CvsB = cms.string("bDiscriminator('pfCombinedCvsBJetTags')"),
            ),
            exprsI = cms.untracked.PSet(
                #q = cms.string("charge"),
            ),
        ),
    ),
)

