import FWCore.ParameterSet.Config as cms

#Scraping event filter
from CATTools.CatProducer.eventCleaning.scrapingFilter_cfi import *

#good primary vertex filter
from CATTools.CatProducer.eventCleaning.primaryVertexFilter_cfi import *
 
# CSC Beam Halo Filter
# Is this really required?
# from RecoMET.METAnalyzers.CSCHaloFilter_cfi import *

# HB + HE noise filter
from CommonTools.RecoAlgos.HBHENoiseFilter_cfi import *

# HCAL laser events filter
from RecoMET.METFilters.hcalLaserEventFilter_cfi import *
hcalLaserEventFilter.vetoByRunEventNumber=cms.untracked.bool(False)
hcalLaserEventFilter.vetoByHBHEOccupancy=cms.untracked.bool(True)

# ECAL dead cell filter
from RecoMET.METFilters.EcalDeadCellTriggerPrimitiveFilter_cfi import *
EcalDeadCellTriggerPrimitiveFilter.tpDigiCollection = cms.InputTag("ecalTPSkimNA")

# Tracking failure filter
from RecoMET.METFilters.trackingFailureFilter_cfi import *
trackingFailureFilter.VertexSource = cms.InputTag("goodOfflinePrimaryVertices", "", "")

# The EE bad SuperCrystal filter
from RecoMET.METFilters.eeBadScFilter_cfi import *

eventCleaning = cms.Sequence(
#		CSCTightHaloFilter*
                HBHENoiseFilter*
                scrapingFilter*
                hcalLaserEventFilter*
                EcalDeadCellTriggerPrimitiveFilter*
#                primaryVertexFilter * # this will be applied before this sequence. do we need to apply this only for data?
                trackingFailureFilter*
                eeBadScFilter
)
