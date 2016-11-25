import FWCore.ParameterSet.Config as cms

def enableAdditionalMETFilters(process, runOnMC=True):
    #######################################################################
    ## Event filters from MET https://twiki.cern.ch/twiki/bin/viewauth/CMS/MissingETOptionalFiltersRun2
    ## New muon filters to be run on the fly
    process.load('RecoMET.METFilters.BadPFMuonFilter_cfi')
    process.BadPFMuonFilter.muons = cms.InputTag("slimmedMuons")
    process.BadPFMuonFilter.PFCandidates = cms.InputTag("packedPFCandidates")

    process.load('RecoMET.METFilters.BadChargedCandidateFilter_cfi')
    process.BadChargedCandidateFilter.muons = cms.InputTag("slimmedMuons")
    process.BadChargedCandidateFilter.PFCandidates = cms.InputTag("packedPFCandidates")
 
    if not hasattr(process, 'nEventsFiltered'):
        process.nEventsFiltered = cms.EDProducer("EventCountProducer")
        process.p += process.nEventsFiltered
    process.p += (process.BadPFMuonFilter*process.BadChargedCandidateFilter*process.nEventsFiltered)

    return process
