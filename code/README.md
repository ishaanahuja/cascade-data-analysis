# PhD thesis: Multi-Strange particle production in p-Pb collisions at √s = 8.16 TeV
- Author: Ishaan Ahuja
- Supervisor: Marek Bombara

## Running Locally

Clone project and navigate to this folder.

1. Start `cvmfs`
```shell
systemctl start cvmfs
```

2. Enter CentOS 7 toolbox
```shell
toolbox enter c7
```

3. Login to AliROOT and setup ALICE environment
```shell
source /cvmfs/alice.cern.ch/etc/login.sh && alienv enter VO_ALICE@AliPhysics::vAN-20210223_ROOT6-1
```
4.  Run analysis locally
```shell
root runAnalysis_ROOT6.C
```
The results are stored in `AnalysisResults.root`


### Note
1. For more info on running analysis tasks, please see [here](https://gitlab.science.upjs.sk/jf/support/-/blob/master/doc/ALICE/README.md).
