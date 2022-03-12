#!/bin/bash
##############################################################################
# Setup script for CEPCSW:
# - build the cepcsw
#
# Usage:
# $ bash build.sh
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

    if [ -n "${bldtool}" ]; then
        blddir=${blddir}.${bldtool}
    fi

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

function run-cmake() {
    local blddir=$(build-dir)

    if [ ! -d "$blddir" ]; then
        error: "Failed to create $blddir"
        return 1
    fi
    local bldgen
    if [ -n "$bldtool" ] ; then
        case $bldtool in
            ninja)
                bldgen="-G Ninja"
                ;;
            *)
                ;;
        esac
    fi

    cd $blddir

    cmake .. -DHOST_BINARY_TAG=${lcg_platform} ${bldgen} || {
        error: "Failed to cmake"
        return 1
    }

}

function run-make() {
    cmake --build .
}

##############################################################################
# Parse the command line options
##############################################################################

# The current default platform
lcg_platform=x86_64-centos7-gcc8-opt
lcg_version=101.0.0

bldtool=${CEPCSW_BLDTOOL} # make, ninja # set in env var

check-cepcsw-envvar || exit -1

check-working-builddir || exit -1

run-cmake || exit -1

run-make || exit -1
