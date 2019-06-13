import FWCore.ParameterSet.Config as cms

## for egamma pid https://twiki.cern.ch/twiki/bin/viewauth/CMS/CutBasedElectronIdentificationRun2#Recipe_for_regular_users_for_8_0
def enableElectronVID(process):
    from PhysicsTools.SelectorUtils.tools.vid_id_tools import switchOnVIDElectronIdProducer, DataFormat, setupAllVIDIdsInModule, setupVIDElectronSelection
    switchOnVIDElectronIdProducer(process, DataFormat.MiniAOD)

    electron_ids = [
        #'mvaElectronID_Fall17_noIso_V2_cff',
        #'mvaElectronID_Fall17_iso_V2_cff',
        'cutBasedElectronID_Fall17_94X_V2_cff',
    ]
    for idmod in electron_ids:
        idmod = "RecoEgamma.ElectronIdentification.Identification."+idmod
        setupAllVIDIdsInModule(process,idmod,setupVIDElectronSelection)

    electron_idNames = [
        #"egmGsfElectronIDs:mvaEleID-Fall17-noIso-V2-wp90",
        #"egmGsfElectronIDs:mvaEleID-Fall17-noIso-V2-wp80",
        #"egmGsfElectronIDs:mvaEleID-Fall17-noIso-V2-wpLoose",
        #"egmGsfElectronIDs:mvaEleID-Fall17-iso-V2-wp90",
        #"egmGsfElectronIDs:mvaEleID-Fall17-iso-V2-wp80",
        #"egmGsfElectronIDs:mvaEleID-Fall17-iso-V2-wpLoose",
        #"egmGsfElectronIDs:mvaEleID-Fall17-iso-V2-wpHZZ",
        "egmGsfElectronIDs:cutBasedElectronID-Fall17-94X-V2-veto",
        "egmGsfElectronIDs:cutBasedElectronID-Fall17-94X-V2-loose",
        "egmGsfElectronIDs:cutBasedElectronID-Fall17-94X-V2-medium",
        "egmGsfElectronIDs:cutBasedElectronID-Fall17-94X-V2-tight",
    ]
    process.catElectrons.electronIDSources = cms.PSet()
    for idName in electron_idNames:
        setattr(process.catElectrons.electronIDSources,
                idName.split(':',1)[1].replace('-','_'),
                cms.InputTag(idName))
        process.catElectrons.electronIDs.append(idName.split(':',1)[1])

    process.egmGsfElectronIDs.physicsObjectSrc = process.catElectrons.unsmaredElectrons
    #process.electronMVAValueMapProducer.srcMiniAOD = process.catElectrons.unsmaredElectrons

    return process

def enablePhotonVID(process):
    from PhysicsTools.SelectorUtils.tools.vid_id_tools import switchOnVIDPhotonIdProducer, DataFormat, setupAllVIDIdsInModule, setupVIDPhotonSelection
    switchOnVIDPhotonIdProducer(process, DataFormat.MiniAOD)

    #Use IDs from https://github.com/cms-sw/cmssw/blob/CMSSW_7_4_X/RecoEgamma/PhotonIdentification/python/Identification/ 
    photon_ids = [
        'mvaPhotonID_Fall17_94X_V2_cff',
    ]
    for idmod in photon_ids:
        idmod = 'RecoEgamma.PhotonIdentification.Identification.'+idmod
        setupAllVIDIdsInModule(process,idmod,setupVIDPhotonSelection)

    photon_idNames = [
        "egmPhotonIDs:mvaPhoID-RunIIFall17-v2-wp80",
        "egmPhotonIDs:mvaPhoID-RunIIFall17-v2-wp90",
    ]
    process.catPhotons.photonIDSources = cms.PSet()
    for idName in photon_idNames:
        setattr(process.catPhotons.photonIDSources,
                idName.split(':',1)[1].replace('-','_'),
                cms.InputTag(idName))
        process.catPhotons.photonIDs.append(idName.split(':',1)[1])

    process.egmPhotonIDs.physicsObjectSrc = process.catPhotons.unsmearedPhotons
    process.photonMVAValueMapProducer.srcMiniAOD = process.catPhotons.unsmearedPhotons
    process.photonIDValueMapProducer.srcMiniAOD = process.catPhotons.unsmearedPhotons

    return process

