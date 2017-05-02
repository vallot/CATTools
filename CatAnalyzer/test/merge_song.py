import sys,os

#mclist = ['TT_powheg', 'WJets', 'SingleTbar_tW', 'SingleTop_tW', 'ZZ', 'WW', 'WZ', 'DYJets', 'DYJets_10to50']
mclist = ['TT_powheg', 'WJets', 'ZZ', 'WW', 'WZ', 'DYJets', 'DYJets_10to50']
datalist = ['DoubleEG','DoubleMuon','MuonEG']
period = ['B','C','D','E','F','G','H_v2','H_v3']


for d in mclist:
    haddMC="hadd /xrootd/store/user/king11kr/ntuples_TtbarDstar_v806/CMSSW_8_0_26_patch1/%s.root /xrootd/store/user/king11kr/ntuples_TtbarDstar_v806/CMSSW_8_0_26_patch1/%s/cattree*.root"%(d,d)
    #print haddMC
    os.system(haddMC)

for d in datalist:
    for p in period:
        haddRD="hadd /xrootd/store/user/king11kr/ntuples_TtbarDstar_v806/CMSSW_8_0_26_patch1/%s_Run2016%s.root /xrootd/store/user/king11kr/ntuples_TtbarDstar_v806/CMSSW_8_0_26_patch1/%s_Run2016%s/cattree*.root"%(d,p,d,p)
        #print haddRD
        os.system(haddRD)
