# QuickStart

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
aliroot runAnalysis_ROOT6.C
```
The results are stored in `CorrelationTask.root`

5. To generate sibling/mixing ratio and projection histograms, run:
```shell
root CorrelationProjections.C
``` 
The output is stored in `CorrelationProjections.root`

### Note
1. By default, caching is enabled. This will download and store ~8GB PbPbLHC11h data to your local storage for re-running the macro. If you wish to disable caching, comment out line #10 in `runAnalysis_ROOT6.C`:
```shell
    TFile::SetCacheFileDir(gSystem->HomeDirectory(), 1, 1);
```
2. For more info on running analysis tasks, please see [here](https://gitlab.science.upjs.sk/jf/support/-/blob/master/doc/ALICE/README.md).