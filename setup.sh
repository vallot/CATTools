scram project -n cattools CMSSW_9_4_9
cd cattools/src
cmsenv
git-cms-init -q
git checkout -b cat90x

git cms-addpkg PhysicsTools/PatAlgos
git cms-merge-topic cms-met:METFixEE2017_949_v2
#git-cms-addpkg RecoEgamma/ElectronIdentification
#git-cms-addpkg EgammaAnalysis/ElectronTools
#git-cms-addpkg RecoMET/METFilters
#git-cms-addpkg RecoBTag/DeepFlavour

git clone https://github.com/vallot/CATTools -b cat90x --recurse-submodules -j8
cd CATTools/CommonTools/scripts
git checkout master
cd ../../..

## Production only - remove large files because of limitation in crab job file size
rm -f CATTools/CatAnalyzer/data/KinReco_input.root
rm -f CATTools/CatAnalyzer/data/KoreaDesyKinRecoInput.root
rm -f CATTools/CatAnalyzer/data/KoreaKinRecoInput_pseudo.root
rm -f CATTools/CatAnalyzer/data/desyKinRecoInput.root
#rm -rf RecoEgamma/*Identification/data/Spring15
#rm -rf RecoEgamma/*Identification/data/PHYS14

scram b -j30

catGetDatasetInfo v9-0-0

## Production only - do the unit test
#cd CATTools/CatProducer
#voms-proxy-init -voms cms
#scram b runtests
