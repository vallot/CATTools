import FWCore.ParameterSet.Config as cms

## for egamma pid https://twiki.cern.ch/twiki/bin/viewauth/CMS/CutBasedElectronIdentificationRun2#Recipe_for_regular_users_for_8_0
def enableElectronNoIsoID(process):
    cutBasedIds = []
    for idSet in process.egmGsfElectronIDs.physicsObjectIDs:
        idName = idSet.idDefinition.idName.value()
        if not idName.startswith('cutBasedElectronID-'): continue
        cutBasedIds.append(idSet.clone())

    for idSet in cutBasedIds:
        newIdName = idSet.idDefinition.idName.value()+'-noiso'
        newCutFlow = cms.VPSet()

        for cut in idSet.idDefinition.cutFlow:
            newCut = cut.clone()
            if newCut.cutName.value().endswith('IsoCut'): newCut.isIgnored = True
            newCutFlow.append(newCut)

        newIdSet = cms.PSet(
            idDefinition = cms.PSet(
                cutFlow = newCutFlow,
                idName = cms.string(newIdName),
            ),
            idMD5 = cms.string(''),
            isPOGApproved = cms.untracked.bool(False),
        )

        process.egmGsfElectronIDs.physicsObjectIDs.append(newIdSet)

        setattr(process.catElectrons.electronIDSources,
                newIdName.replace('-','_'),
                cms.InputTag('egmGsfElectronIDs', newIdName))
        process.catElectrons.electronIDs.append(newIdName)

    return process

