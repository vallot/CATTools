## Histogram based ttbar dilepton analysis chain
The collection in this directory is aiming to enable most of analysis chains
based on static fine-binned histograms.

**Analyzer module** Histogram definitions and selection
cuts are almost hard-coded and controled from a single analyzer module.
The corresponding run configuration file is also kept simple as possible.
But for the future extensions, the analyzer is actually inherited from the
EDFilter and the cfg is configured to store events which passes certain
cut flow step. All histograms are very fine binned, 1GeV for mass, momentum
distributions and ~0.1 or 0.01 for eta, phi like plots. Users can rebin
the histograms after own binning study.

**Workflow** The helper scripts are named to stand for the steps.
Start from the pass1.1_eventSel.py to prepare job submission. The
script will indicate what to do in the next step - submit jobs to KISTI
cluster using the auto-generated submittor shell script. Follow pass1.2,
pass2.1, etc after checking each steps are finished without serious problems.
Users are allows to modify run scripts on demand. At the final step of
the workflow using the script, user will have necessary control plots
including the systematic error calculation.

**Extended study** The initial cfg file is configured to store events
in the EDM format. By default, CAT-objects which are directly used
in the reference event selection are kept in the event content.
Event files are not stored for the systematic uncertainty variation
samples. Users can do the remaining analysis based on these objects
such as kinematic reconstruction, complicated rescailing, etc.

