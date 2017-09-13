scram project -n cattools CMSSW_8_0_29
cd cattools/src
cmsenv
git-cms-init -q
git checkout -b cat80x

git-cms-addpkg RecoEgamma/ElectronIdentification
git-cms-addpkg EgammaAnalysis/ElectronTools
git-cms-addpkg RecoMET/METFilters
git-cms-addpkg RecoBTag/DeepFlavour

sed -i 's/badGlobalMuonTagger.clone/badGlobalMuonTaggerMAOD.clone/g' RecoMET/METFilters/python/badGlobalMuonTaggersMiniAOD_cff.py

#brings in HEEP V70 into VID
git cms-merge-topic -u Sam-Harper:HEEPV70VID_8010_ReducedCheckout

#for other E/gamma IDs in VID if you wish to have them
git cms-merge-topic -u ikrav:egm_id_80X_v3

#only necessary to run HEEP V70 on AOD
git cms-merge-topic -u Sam-Harper:PackedCandNoPuppi 
git-cms-merge-topic -u cms-egamma:EGM_gain_v1
git-cms-merge-topic -u ikrav:egm_id_80X_v3_photons

cd EgammaAnalysis/ElectronTools/data
wget https://github.com/cms-egamma/RegressionDatabase/blob/master/SQLiteFiles/GED_80X_Winter2016/ged_regression_20170114.db?raw=true -O ged_regression_20170114.db
git clone -b Moriond17_gainSwitch_unc https://github.com/ECALELFS/ScalesSmearings.git
cd ../../..

mkdir -p RecoBTag/DeepFlavour/data
wget http://home.fnal.gov/~verzetti//DeepFlavour/training/DeepFlavourNoSL.json -O  RecoBTag/DeepFlavour/data/DeepFlavourNoSL.json
wget http://mon.iihe.ac.be/~smoortga/DeepFlavour/CMSSW_implementation_DeepCMVA/Model_DeepCMVA.json -O RecoBTag/DeepFlavour/data/Model_DeepCMVA.json

git clone -b egm_id_80X_v1 https://github.com/ikrav/RecoEgamma-ElectronIdentification.git data1 
git clone -b egm_id_80X_v1 https://github.com/ikrav/RecoEgamma-PhotonIdentification.git data2
git clone https://github.com/cms-data/RecoEgamma-ElectronIdentification data3

rsync -avz data1/* RecoEgamma/ElectronIdentification/data/
rsync -avz data2/* RecoEgamma/ElectronIdentification/data/
rsync -avz data3/* RecoEgamma/PhotonIdentification/data/
rm -rf data1 data2 data3

git clone https://github.com/vallot/CATTools
cd CATTools
git checkout -b v8-0-7 v8-0-7
git submodule init
git submodule update
cd ..

## Production only - remove large files because of limitation in crab job file size
rm -f CATTools/CatAnalyzer/data/KinReco_input.root
rm -f CATTools/CatAnalyzer/data/KoreaDesyKinRecoInput.root
rm -f CATTools/CatAnalyzer/data/KoreaKinRecoInput_pseudo.root
rm -f CATTools/CatAnalyzer/data/desyKinRecoInput.root
rm -rf RecoEgamma/*Identification/data/Spring15
rm -rf RecoEgamma/*Identification/data/PHYS14

scram b -j30

catGetDatasetInfo

## Production only - do the unit test
#cd CATTools/CatProducer
#voms-proxy-init -voms cms
#scram b runtests
