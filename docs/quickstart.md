# Quick start

## Start development environment in a Singularity Container

Start the container in lxslc7 (OS: CentOS7):
```
$ /cvmfs/container.ihep.ac.cn/bin/hep_container shell SL6
```

## Get code using git

```
$ git clone git@github.com:cepc/CEPCSW.git
$ cd CEPCSW
```
## Setup and build 

```
$ source setup.sh
$ ./build.sh
```

## Run a simple test

```
$ ./run.sh Examples/options/tut_detsim_SDT.py
```

