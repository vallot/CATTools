from PhysicsTools.PatAlgos.patTemplate_cfg import *

catJetsSource = "selectedPatJetsPFlow"
catGenJetsSource = "ak5GenJets"
catMuonsSource = "selectedPatMuonsPFlow"
catElectronsSource = "selectedPatElectronsPFlow"
catPhotonsSource = "selectedPatPhotons"
catMETsSource = "patMETsPFlow"
catVertexSource = "offlinePrimaryVertices"
catMCsource = "genParticles"
catBeamSpot = "offlineBeamSpot"

def catTool(process):
    process.load("CATTools.CatProducer.eventCleaning.eventCleaning_cff")
    process.load("CATTools.CatProducer.catCandidates_cff")

    print "process.catJets.src",process.catJets.src
    process.p += process.eventCleaning+process.makeCatCandidates

    ## electron ID tool
    process.load('EgammaAnalysis.ElectronTools.electronIdMVAProducer_cfi')
    process.eidMVASequence = cms.Sequence(  process.mvaTrigV0 + process.mvaNonTrigV0 )
    process.patElectronsPFlow.electronIDSources.mvaTrigV0    = cms.InputTag("mvaTrigV0")
    process.patElectronsPFlow.electronIDSources.mvaNonTrigV0 = cms.InputTag("mvaNonTrigV0")
    process.patDefaultSequencePFlow.replace( process.patElectronsPFlow, process.eidMVASequence * process.patElectronsPFlow )

    ## cuts on selected Pat objects
    getattr(process,catJetsSource).cut = cms.string("pt > 20")
    getattr(process,catMuonsSource).cut = cms.string("pt > 5 || isPFMuon || (pt > 3 && (isGlobalMuon || isStandAloneMuon || numberOfMatches > 0 || muonID('RPCMuLoose')))") 
    getattr(process,catElectronsSource).cut = cms.string("pt > 5") 
    getattr(process,catPhotonsSource).cut = cms.string("pt > 5")

    # no taus for now...
    #process.selectedPatTausPFlow.cut = cms.string("pt > 18. && tauID('decayModeFinding')> 0.5")

    process.patJetPartonMatchPFlow.mcStatus = [ 3, 23 ]
    process.patPFParticlesPFlow.embedGenMatch = cms.bool(True)

    ## adding in user variables ## temp to add into objects
    process.patMuonsPFlow.isolationValues.user = cms.VInputTag("muPFIsoValueCharged03PFlow","muPFIsoValueNeutral03PFlow","muPFIsoValueGamma03PFlow","muPFIsoValuePU03PFlow","muPFIsoValueChargedAll03PFlow")

    process.patElectronsPFlow.isolationValues.user = cms.VInputTag("elPFIsoValueCharged03PFIdPFlow","elPFIsoValueNeutral03PFIdPFlow","elPFIsoValueGamma03PFIdPFlow","elPFIsoValuePU03PFIdPFlow","elPFIsoValueChargedAll03PFIdPFlow")

    process.patJetsPFlow.addTagInfos = cms.bool(True)
    process.patJetsPFlow.userData.userFunctions = cms.vstring( "? hasTagInfo('secondaryVertex') && tagInfoSecondaryVertex('secondaryVertex').nVertices() > 0 ? "
"tagInfoSecondaryVertex('secondaryVertex').secondaryVertex(0).p4().mass() : 0",
"? hasTagInfo('secondaryVertex') && tagInfoSecondaryVertex('secondaryVertex').nVertices() > 0 ? "
"tagInfoSecondaryVertex('secondaryVertex').flightDistance(0).value() : 0",
"? hasTagInfo('secondaryVertex') && tagInfoSecondaryVertex('secondaryVertex').nVertices() > 0 ? "
"tagInfoSecondaryVertex('secondaryVertex').flightDistance(0).error() : 0",
)
    process.patJetsPFlow.userData.userFunctionLabels = cms.vstring('secvtxMass','Lxy','LxyErr')

