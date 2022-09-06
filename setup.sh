#!/bin/bash
##############################################################################
# Setup script for CEPCSW:
# - setup the external libraries
#
# Usage:
# $ source setup.sh
#
# Author: Tao Lin <lintao@ihep.ac.cn>
##############################################################################

THISSCRITDIR=$(dirname $(readlink -e "${BASH_SOURCE[0]}" 2>/dev/null) 2>/dev/null) # Darwin readlink doesnt accept -e

function info:() {
    echo "INFO: $*" 1>&2
}

function error:() {
    echo "ERROR: $*" 1>&2
}

function lcg-version() {
    echo ${CEPCSW_LCG_VERSION}
}

function cepcsw-base() {
    echo /cvmfs/cepcsw.ihep.ac.cn/prototype
}

function cepcsw-external() {
    echo $(cepcsw-base)/releases/externals/$(lcg-version)
}

function setup-external() {
    local f=$(cepcsw-external)/setup-${CEPCSW_LCG_VERSION}.sh
    if [ ! -f $f ]; then
        error: "Failed to find setup script $f"
        return
    fi

    info: "Setup CEPCSW externals: $f"
    source $f

}

function setup-install-area() {
    local installarea=$THISSCRITDIR/InstallArea
    if [ ! -d "$installarea" ]; then
        return
    fi

    export PATH=$installarea/bin:$PATH
    export LD_LIBRARY_PATH=$installarea/lib:$LD_LIBRARY_PATH
    export PYTHONPATH=$installarea/lib:$PYTHONPATH
    export PYTHONPATH=$installarea/python:$PYTHONPATH
    export ROOT_INCLUDE_PATH=$installarea/include:$ROOT_INCLUDE_PATH
    info: "Setup CEPCSW: $installarea"
}

##############################################################################
# Parse the command line options
##############################################################################

# CEPCSW_LCG_VERSION=${1}; shift

if [ -z "$CEPCSW_LCG_VERSION" ]; then
    CEPCSW_LCG_VERSION=101.0.1
fi
export CEPCSW_LCG_VERSION

setup-external
setup-install-area
