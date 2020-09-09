# [CEPCSW](https://cepc.github.io/CEPCSW/)

CEPC offline software prototype based on [Key4hep](https://github.com/key4hep).

## Quick start

Start the container in lxslc7 (OS: CentOS7):
```
$ /cvmfs/container.ihep.ac.cn/bin/hep_container shell SL6
```

Before run following commands, please make sure you setup the CVMFS:

```
$ git clone git@github.com:cepc/CEPCSW.git
$ cd CEPCSW
$ git checkout master # branch name
$ source setup.sh
$ ./build.sh
$ ./run.sh Examples/options/helloalg.py
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
