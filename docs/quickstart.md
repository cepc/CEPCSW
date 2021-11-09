# Quick start

## Start development environment in lxslc7 (CentOS 7)

```
$ ssh -Y USERNAME@lxslc7.ihep.ac.cn
```

## Manage code using git

Fork the CEPCSW into your own repo. For an example:
* https://github.com/*USERNAME*/CEPCSW

Get the source code from your own repo:
```
$ git clone git@github.com:USERNAME/CEPCSW.git
$ cd CEPCSW
```

Add the upstream repo:
```
$ git remote add cepc https://github.com/cepc/CEPCSW.git
```

Sync and merge the upstream repo:
```
$ git remote update
$ git merge cepc/master
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

