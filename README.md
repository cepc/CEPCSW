# CEPCSW

CEPC offline software prototype based on [Key4hep](https://github.com/key4hep).

## Quick start

Before run following commands, please make sure you setup the CVMFS:

```
$ source /cvmfs/cepcsw.ihep.ac.cn/prototype/releases/externals/97.0.2/setup.sh
$ git clone git@github.com:cepc/CEPCSW.git
$ cd CEPCSW
$ git checkout master # branch name
$ mkdir build && cd build
$ cmake .. -DHOST_BINARY_TAG=${BINARY_TAG}
$ make
$ ./run gaudirun.py '$EXAMPLESROOT/options/helloalg.py'
```

## Packages

* Examples: For new comers and users

* Detector: Geometry

* Generator: Physics Generator

* Simulation: Detector Simulation

* Reconstruction: Reconstruction

* Service: Common Service

## Full Chain

Detector simulation: 
```
$ ./run gaudirun.py '$EXAMPLESROOT/options/tut_detsim.py'
```

## Conventions for collections
Keep the collection names compatible between the prototype and the existing CEPC software.

* MCParticle
* VXDCollection
* SITCollection
* TPCCollection
* SETCollection
