import ROOT
ROOT.gSystem.Load("libFWCoreFWLite.so");
ROOT.gSystem.Load("libDataFormatsFWLite.so");
ROOT.AutoLibraryLoader.enable()
from DataFormats.FWLite import Events, Handle
import os,sys
    
events = Events("$CMSSW_BASE/src/CATTools/CatProducer/prod/catTuple.root")

cut_lep_pt = 20.
cut_lep_eta = 2.1
cut_lep_iso = 10.12
cut_jet_pt = 30.
cut_jet_eta = 4.7

metsLabel, mets = "catMETs", Handle("std::vector<cat::MET>")
jetsLabel, jets = "catJets", Handle("std::vector<cat::Jet>")
muonsLabel, muons = "catMuons", Handle("std::vector<cat::Muon>")
electronsLabel, electrons = "catElectrons", Handle("std::vector<cat::Electron>")
vertexsLabel, vertexs = "offlineSlimmedPrimaryVertices", Handle("std::vector<reco::Vertex>")
nmuons = 0
njets = 0
npass = 0
chems=[0, 0, 0, 0, 0]
chmms=[0, 0, 0, 0, 0]
chees=[0, 0, 0, 0, 0]
totaliev = 0
for iev,event in enumerate(events):
    if iev%10000 == 0:
        print round(iev/50000.*100), "%"
    selectedelecs=[]
    selectedmuons=[]
    selectedjets=[]

    event.getByLabel(vertexsLabel, vertexs)
    vert = vertexs.product()[0]
    
    event.getByLabel(electronsLabel,electrons)
    for g,m in enumerate(electrons.product()):
        if not m.electronID("cutBasedElectronID-PHYS14-PU20bx25-V2-standalone-medium"):
            continue
        if not m.passConversionVeto():
            continue
        if not m.isPF():
            continue
#        if m.gsfTrack.hitPattern.numberOfLostHits('MISSING_INNER_HITS') > 0:
#            continue
        if m.pt() <= cut_lep_pt:
            continue
        if abs(m.scEta()) <= 1.4442:
            if m.relIso(0.3) >= 0.1649:
                continue
        if abs(m.scEta()) >= 1.566:
            if m.relIso(0.3) >= 0.2075:
                continue
        if abs(m.scEta()) > 1.4442 and abs(m.scEta()) < 1.566:
            continue
        if abs(m.eta()) >= 2.5:
            continue
        selectedelecs.append(m)
       
    event.getByLabel(muonsLabel,muons)
    for g,m in enumerate(muons.product()):
        if not m.isTightMuon():
            continue
        if m.pt() <= cut_lep_pt:
            continue
        if abs(m.eta()) >= 2.4:
            continue
        if m.relIso(0.4) >= 0.12:
            continue
        selectedmuons.append(m)
        nmuons +=1


#step1
    lep1 = ROOT.TLorentzVector()
    lep2 = ROOT.TLorentzVector()
    ch=[]
    charge = 1
    if len(selectedmuons) + len(selectedelecs) == 2:
        if len(selectedmuons) == 1:
            lep1.SetPtEtaPhiM(selectedmuons[0].pt(), selectedmuons[0].eta(), selectedmuons[0].phi(), 0.1056)
            lep2.SetPtEtaPhiM(selectedelecs[0].pt(), selectedelecs[0].eta(), selectedelecs[0].phi(), 0.00051)
            charge = selectedmuons[0].charge()*selectedelecs[0].charge()
            ch = "chem"
        elif len(selectedmuons) == 2:
            lep1.SetPtEtaPhiM(selectedmuons[0].pt(), selectedmuons[0].eta(), selectedmuons[0].phi(), 0.1056)
            lep2.SetPtEtaPhiM(selectedmuons[1].pt(), selectedmuons[1].eta(), selectedmuons[1].phi(), 0.1056)
            charge = selectedmuons[0].charge()*selectedmuons[1].charge()
            ch = "chmm"
        else:
            lep1.SetPtEtaPhiM(selectedelecs[0].pt(), selectedelecs[0].eta(), selectedelecs[0].phi(), 0.00051)
            lep2.SetPtEtaPhiM(selectedelecs[1].pt(), selectedelecs[1].eta(), selectedelecs[1].phi(), 0.00051)
            charge = selectedelecs[0].charge()*selectedelecs[1].charge()
            ch = "chee"
        
    if (lep1+lep2).M() <= 20:
        continue
    if charge > 0:
        continue

    if ch == "chem": chems[0] +=1
    if ch == "chmm": chmms[0] +=1
    if ch == "chee": chees[0] +=1

#step2
    if ch == []:
        continue
    if ch != "chem":
        if (lep1+lep2).M() > 76 and (lep1+lep2).M() < 106:
            continue

    if ch == "chem": chems[1] +=1
    if ch == "chmm": chmms[1] +=1
    if ch == "chee": chees[1] +=1

#step3    
    event.getByLabel(jetsLabel,jets)
    for g,j in enumerate(jets.product()):
        if not j.LooseId():
            continue
        #if j.pileupJetId() < 0.9:
        #    continue
        if j.pt() <= cut_jet_pt:
            continue
        if abs(j.eta()) >= 2.4:
            continue 
        jet = ROOT.TLorentzVector(j.px(), j.py(), j.pz(), j.energy())
        dr1 = jet.DeltaR(lep1)
        dr2 = jet.DeltaR(lep2)
        if dr1 <= 0.4:
            continue
        if dr2 <= 0.4:
            continue
        selectedjets.append(j)
        njets += 1

    if len(selectedjets) < 2:
        continue

    if ch == "chem": chems[2] +=1
    if ch == "chmm": chmms[2] +=1
    if ch == "chee": chees[2] +=1

#step4
    event.getByLabel(metsLabel,mets)
    met = mets.product()[0]

    if ch == "chem": chems[3] +=1

    if ch == "chmm" or ch == "chee":
        if met.pt() <= 40.:
            continue
    if ch == "chmm": chmms[3] +=1
    if ch == "chee": chees[3] +=1

#step5
    btag = 0
    for g,j in enumerate(selectedjets):
        jets_CSVInclV2 = j.bDiscriminator("pfCombinedInclusiveSecondaryVertexV2BJetTags")
        if jets_CSVInclV2 <= 0.814:
            continue
        btag += 1

    if btag == 0:
        continue
 
    if ch == "chem": chems[4] +=1
    if ch == "chmm": chmms[4] +=1
    if ch == "chee": chees[4] +=1

    """
    print "met pt", met.pt(), "met eta", met.eta(), "met phi", met.phi()

    if len(selectedmuons) < 2:
        continue

    for m in selectedmuons:
        print "muon pt", m.pt(), "muon eta", m.eta(), "muon charge", m.charge()
            
    if selectedmuons[0].charge()*selectedmuons[1].charge() == 1:
        continue
    mu1 = selectedmuons[0]
    mu2 = selectedmuons[1]
    
    higgs = ROOT.TLorentzVector(mu1.px(), mu1.py(), mu1.pz(), mu1.energy()) + ROOT.TLorentzVector(mu2.px(), mu2.py(), mu2.pz(), mu2.energy())

    print higgs.M()
    
    npass +=1
    """

print
print "<num of passed>"
print "    ", "s1 ", "s2 ", "s3 ", "s4 ", "s5 "
print "chee" , 
for i in range(0,5):
    print chees[i],
print
print "chem" , 
for i in range(0,5):
    print chems[i],
print
print "chmm" , 
for i in range(0,5):
    print chmms[i],
print
print

