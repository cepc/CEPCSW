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


##############################################################################
# Parse the command line options
##############################################################################

# CEPCSW_LCG_VERSION=${1}; shift

if [ -z "$CEPCSW_LCG_VERSION" ]; then
    CEPCSW_LCG_VERSION=98.0.0
fi
export CEPCSW_LCG_VERSION

setup-external
