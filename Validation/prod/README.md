## Histogram based ttbar dilepton analysis chain
The collection in this directory is aiming to enable most of analysis chains
based on static fine-binned histograms.

### Histogram definitions and selection ###
cuts are almost hard-coded and controled from a single analyzer module.
The corresponding run configuration file is also kept simple as possible.
But for the future extensions, the analyzer is actually inherited from the
EDFilter and the cfg is configured to store events which passes certain
cut flow step. All histograms are very fine binned, 1GeV for mass, momentum
distributions and ~0.1 or 0.01 for eta, phi like plots. Users can rebin
the histograms after own binning study.

For a details of the analyzer module see the source codes:
  * Analysis module : CATTools/CatAnalyzer/plugins/TTLLEventSelector.cc
  * Initial configuration : CATTools/CatAnalyzer/python/ttll/ttllEventSelector\_cfi.py
  * Working run configuration : "CATTools/CatAnalyzer/test/TTLL/analyze\_\*\_cfg.py"

### Starting the workflow ###
You can find set of scripts to run step-by-step.

Start from the pass1.1\_eventSel.py to prepare job submission. The
script will indicate what to do in the next step - submit jobs to KISTI
cluster using the auto-generated submittor shell script. Follow pass1.2,
pass2.1, etc after checking each steps are finished without serious problems.
Users are allows to modify run scripts on demand. At the final step of
the workflow using the script, user will have necessary control plots
including the systematic error calculation.

#### The Histogram (and selected event data) production ####
Start the histogram production and control plots.
Go to the working directory, CATTools/CatAnalyzer/test/TTLL,
and run pass1.1\_eventSel.py script to prepare job submission.
  * CATTools/CatAnalyzer/test/TTLL/pass1.1\_eventSel.py

Running this script will create a directory "pass1" which contains simple bash script
to call create-batch tool for each samples times all available systematic uncertainty variations.
```bash
cd pass1
cat submit*.sh | sed -e 's;create-batch;-n ;g' | xargs -L1 -P20 create-batch
for i in */central; do pushd $i; ./submit.sh; popd; done
for i in */*/up; do pushd $i; ./submit.sh; popd; done
for i in */*/dn; do pushd $i; ./submit.sh; popd; done
for i in */gen_scale/*; do pushd $i; ./submit.sh; popd; done
for i in */gen_PDF/*;do pushd $i; ./submit.sh; popd; done
#./submit_sig_central.sh &
#./submit_bkg_central.sh &
#./submit_data_central.sh &

#./submit_sig_unc.sh &
#./submit_bkg_unc.sh &
#./submit_data_central.sh &
cd ..
```

The number of jobs leaches ~5000, so you need to be patient.
Please also note that you need large disk spaces to submit everything in one shot,
more than 100GBytes to store all output and temporary files.

You can clean up temporary files if a job cluster is finished by runing the second script,
```bash
./pass1.2_cleanup.py
```
Histograms are merged into one and temporary files for job submission are cleaned up.

For analysis dependent parts, event data with CAT data formats are stored for the events
which passes a cut step (Step 4 nJet in the current default), for the "central" samples.
You can write your own analyzer and configuration files to run on these "reduced" CATTuples
which can be browsed under the "pass1/\*/central" directory with the file name "out\_\*.root".
```bash
find pass1 -name 'central/out\*.root'
```

#### Normalization and combination of physics processes ####
The histogram output from the first pass are normalized to the cross section and event statistics
at the second steps. In addition, some samples that are splited during the production are merged
to one - DYJets10to50 and ZJets samples are merged to one file for Drell-Yan process.
You can skip some samples if you want to ignore by setting "Blacklist" in the script.
```bash
./pass2.1_selectSamples.py
## Edit pass2/samples.json manually if needed.
./pass2.2_scaleMerge.py
```

You can find scaled histograms under the pass2 directory, separated by their variations.

Some data-driven background estimations can be done after the pass2.2 step.

The Drell-Yan scale factor can be extracted using the pass2.3 script,
```bash
./pass2.3_DYEstimation.py
```
will print out the DY scale factors and store the results under the pass2/scaler\_DY.json.

Run the preview step to debug the plots and check the very first cut flow table.
```bash
./pass2.4_preview.py
```

This script produces all Data-MC comparison plots into the pass2/preview.root file
in TCanvas formats. You can dump the plots to image files using the dumpRoot script
which is shipped in the hep-tools package.

#### Common control plot production and basic measurements ####
The production-level plots are produced within the pass3. In this step, historam ranges
and binnings are adjusted to fit to user analyses. All systematic uncertainty variations
are also configured to be ready for the measurement.

Start from the rebining. Edit the pass3.1_rebin.py to be consistent with your analysis.
Run the script when ready.
```bash
./pass3.1_rebin.py
```

There are some systematic uncertainty variations which are left to be combined after the
binning.
```bash
./pass3.2_reduceUnc.py
```

Currently, gen_PDF and gen_scale categories are known to the author and they
are combined by "hessian" and "envelope" method. The results are stored under "up" and
"dn" directory with the same file structure. From this step, all samples for the study
are placed under "central", "up" or "dn".

Now everything is ready. We can proceed to run the RooFit/RooStats to extract
measurements with systematic uncertainty calculation - the output root files
must have HistFactory familiar format.

#### Extended study ####
The initial cfg file is configured to store events in the EDM format.
By default, CAT-objects which are directly used
in the reference event selection are kept in the event content.
Event files are not stored for the systematic uncertainty variation
samples. Users can do the remaining analysis based on these objects
such as kinematic reconstruction, complicated rescailing, etc.

For the kinematic reconstruction, see the producer module
  * Producer module : CATTools/CatAnalyzer/plugins/TTLLKinSolutionProducer.cc
  * corresponding initial config : CATTools/CatAnalyzer/python/ttll/ttbarDileptonKinSolutionProducer\_cfi.py

