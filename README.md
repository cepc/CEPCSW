# CEPCSW

CEPC offline software prototype based on Gaudi.

The core part of this project (CEPCSW/FWCore) is taken from FCCSW/FWCore.
Please refer to https://github.com/HEP-FCC/FCCSW

## Quick start

```
$ source /cvmfs/cepcsw.ihep.ac.cn/prototype/setup.sh
$ git clone git@cepcgit.ihep.ac.cn:cepc-prototype/CEPCSW.git
$ cd CEPCSW
$ mkdir build && cd build
$ cmake ..
$ make
$ ./run gaudirun.py '$EXAMPLESROOT/options/helloalg.py'
```

## Packages

* Core: Core Software Based on Gaudi

* Generator: Physics Generator

* Simulation: Detector Simulation

* Reconstruction: Reconstruction

* Examples: For new comers and users

* cepcenv: tools to help developers build packages

