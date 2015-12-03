import ROOT
def matchFill(matched, idx, genJpsi, catJpsi) :
  gen = ROOT.TLorentzVector()
  old = ROOT.TLorentzVector()
  new = ROOT.TLorentzVector()
  gen.SetPtEtaPhiM( genJpsi.pt(), genJpsi.eta(), genJpsi.phi(), genJpsi.mass())
  if ( idx in matched ) :
    print "already filled"
    old_reco = matched[idx]
    new_reco = catJpsi
    old.SetPtEtaPhiM( old_reco.pt(), old_reco.eta(), old_reco.phi(), old_reco.mass())
    new.SetPtEtaPhiM( new_reco.pt(), new_reco.eta(), new_reco.phi(), new_reco.mass())
    if ( gen.DeltaR( old) > gen.DeltaR( new) ) : matched[idx] = catJpsi
  else : matched[idx] = catJpsi
  return matched

def Matching( genJpsis, catJpsis ) :
  matched = {}
  for gen_idx, genJpsi in enumerate(genJpsis) :
    gen = ROOT.TLorentzVector()
    gen.SetPtEtaPhiM( genJpsi.pt(), genJpsi.eta(), genJpsi.phi(), genJpsi.mass())
    for catJpsi in catJpsis :
      cat = ROOT.TLorentzVector()
      cat.SetPtEtaPhiM( catJpsi.pt(), catJpsi.eta(), catJpsi.phi(), catJpsi.mass())
      dR = gen.DeltaR( cat)
      dPt = abs(cat.Pt() - gen.Pt())/gen.Pt()
      if ( dR<0.15 and dPt<0.05) :
        matched = matchFill( matched, gen_idx, genJpsi, catJpsi)
  return matched
