# [CEPCSW](https://cepc.github.io/CEPCSW/)

[![Build Status](https://www.travis-ci.com/cepc/CEPCSW.svg?branch=master)](https://www.travis-ci.com/cepc/CEPCSW)

CEPC offline software prototype based on [Key4hep](https://github.com/key4hep).

## Quick start

Start an SL6 container in lxslc7 (OS: CentOS7):
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

* Digitization: Digitization

* Reconstruction: Reconstruction


## Conventions for collections
Keep the collection names compatible between the prototype and the existing CEPC software.

* MCParticle
* VXDCollection
* SITCollection
* TPCCollection
* SETCollection
