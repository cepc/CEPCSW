#!/bin/bash
############################################
# Description:
#   Manage the github runners in singularity
# Usage:
#   $ ./setup-github-runner new <TOKEN>
#   $ ./setup-github-runner start
# Author: Tao Lin <lintao AT ihep.ac.cn>
############################################


#############################################
# Configuration
#############################################
export RUNNER_TOP_DIR=/tmp/$USER/github-runner
export SINGULARITY_BINDPATH=/cvmfs
export RUNNER_REPO=https://github.com/cepc/CEPCSW

[ -d "$RUNNER_TOP_DIR" ] || mkdir $RUNNER_TOP_DIR

#############################################
# Create a new github action runner (gar)
#############################################

function gar-new-id() {
    local currentid="$(find $RUNNER_TOP_DIR -maxdepth 1 -name github-runner-\* -type d | rev | cut -d- -f 1 | rev | sort -n | tail -n1)"
    if [ -z "$currentid" ]; then
        echo 1
    else
        echo $((currentid+1))
    fi
}

function gar-new-name() {
    echo github-runner-$(gar-new-id)
}

function gar-download-url() {
    echo https://github.com/actions/runner/releases/download/v2.274.2/actions-runner-linux-x64-2.274.2.tar.gz
}
function gar-download-filename() {
    echo actions-runner-linux-x64-2.274.2.tar.gz
}

function gar-new() {
    local dn=$(gar-new-name)
    local fdn=$RUNNER_TOP_DIR/$dn
    if [ -d "$fdn" ]; then
        echo "ERROR: $dn already exists" 1>&2
        exit -1
    fi

    mkdir $fdn || {
        echo "ERROR: Failed to create $fdn" 1>&2
        exit -1
    }


    pushd $RUNNER_TOP_DIR
    if [ ! -f "$(gar-download-filename)" ]; then
        curl -O -L $(gar-download-url) || exit -1
    fi
    popd

    pushd $fdn

    tar xzf $RUNNER_TOP_DIR/$(gar-download-filename) || exit -1

    # start singularity instance
    singularity instance start ~/github-runner.sif ${dn}
    
    singularity run instance://${dn} ./config.sh --url ${RUNNER_REPO} --token ${token} || exit -1
    singularity run instance://${dn} bash -c "./run.sh &"

    popd
    
}

function new() {

    token=$1; shift
    if [ -z "$token" ]; then
        echo "Please pass the token to this script" 1>&2
        exit -1
    fi
    gar-new
}

#############################################
# Start github action runners (gar)
#############################################

function gar-lists() {
    find $RUNNER_TOP_DIR -maxdepth 1 -name github-runner-\* -type d -exec basename {} \;
}

function gar-check() {
    local gar=$1;
    local result=$(singularity instance list $gar | grep $gar)
    if [ -n "$result" ]; then
        echo Y
    else
        echo N
    fi
}

function gar-start() {
    local gar=$1;

    local isrunning=$(gar-check $gar)
    if [ "$isrunning" = "Y" ]; then
        echo "WARNING: $gar is already running. skip it."
        return
    fi

    pushd $RUNNER_TOP_DIR/$gar
    singularity instance start ~/github-runner.sif ${gar}
    singularity run instance://${gar} bash -c "./run.sh &"
    popd
}

function start() {
    local gars="$*"
    if [ -z "$gars" ]; then
        echo "All the github action runners will be started"
        gars="$(gar-lists)"
    fi
    local gar
    for gar in $gars; do
        gar-start $gar
    done
}

#############################################
# Command line options
#############################################

cmd=$1; shift
if [ -z "$cmd" ]; then
    echo "Please specify the command to be invoked" 1>&2
    exit -1
fi

case $cmd in
    new)
        new $*
        ;;
    start)
        start $*
        ;;
    *)
        echo "Unknown command '$cmd'" 1>&2
        ;;
esac
