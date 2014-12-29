CATTools/TopAnalysis
========

plugins, scripts and macro files for standard top quark analysis with CATTools

## Instruction to analyse ttbar-dilepton channel
(Hopefully) Top quark differential cross section measurement and related topics can be done.
An example cfg file can be found in the test directory.

```bash
cd test
cmsRun run_ttll_cfg.py
```


The "TTbarDileptonProducer" module allows to do standard object combinations with some
kinematic reconstructions. The output of this module is the best solution of the reconstruction
algorithm. Leptons and jets used in the solution is linked as daughter particles (access via edm::Ptr).

Kinematic reconstruction algorithm have to be added by user.
As a result of kinematic reconstructions, each algorithm should store its quality variable.
Note that positive value is treated to be good one, so Chi-square should be stored with additional
negative sign, for example.

By default, dummy algorithm is implemented but other widely used algorithms will be arrived.
- Simple four vector sum
- Neutrino weighting (not added)
- MT2 solution (to be implemented)
- MAOS solution (to be implemented)
