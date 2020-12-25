#!/bin/bash
##############################################################################
# Run script for CEPCSW:
# - run a simple job
#
# Usage:
# $ ./run.sh Examples/options/helloalg.py
# or:
# $ 
#
# Author: Tao Lin <lintao@ihep.ac.cn>
##############################################################################

function info:() {
    echo "INFO: $*" 1>&2
}

function error:() {
    echo "ERROR: $*" 1>&2
}

function check-cepcsw-envvar() {
    # CEPCSWEXTERNAL is defined in /cvmfs/cepcsw.ihep.ac.cn/prototype/releases/externals/
    if [ -z "${CEPCSWEXTERNAL}" ]; then
        error: "The CEPCSW is not setup. Please source setup.sh."
        return 1
    fi
}

function build-dir() {
    local blddir=build

    # If detect the extra env var, append it to the build dir
    if [ -n "${lcg_version}" ]; then
        blddir=${blddir}.${lcg_version}
    fi
    if [ -n "${lcg_platform}" ]; then
        blddir=${blddir}.${lcg_platform}
    fi

    echo $blddir
}

function check-working-builddir() {
    local blddir=$(build-dir)
    if [ ! -d "$blddir" ]; then
        mkdir $blddir || {
            error: "Failed to create $blddir"
            return 1
        }
    fi
}

function run-job() {
    local blddir=$(build-dir)

    $blddir/run gaudirun.py $*
}

##############################################################################
# Parse the command line options
##############################################################################

# The current default platform
lcg_platform=x86_64-centos7-gcc8-opt
lcg_version=98.0.0

check-cepcsw-envvar || exit -1

check-working-builddir || exit -1

run-job $* || exit -1
