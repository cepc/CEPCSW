#!/bin/bash
############################################
# Description:
#   Manage the gitlab runners
# Usage:
#   $ ./setup-gitlab-runner new <TOKEN>
#   $ ./setup-gitlab-runner start
#
# Register in crontab:
#
#   */10 * * * * $HOME/setup-github-runner.sh start >& /tmp/$USER/github-runner/start.log
#
# Author: Tao Lin <lintao AT ihep.ac.cn>
############################################

#############################################
# Configuration
#############################################
export RUNNER_TOP_DIR=/tmp/$USER/gitlab-runner
export SINGULARITY_BINDPATH=/cvmfs
export RUNNER_URL=https://code.ihep.ac.cn

[ -d "$RUNNER_TOP_DIR" ] || mkdir $RUNNER_TOP_DIR

#############################################
# Create a new gitlab runner (glr)
#############################################

# ./gitlab-runner register  --url https://code.ihep.ac.cn  --token XXXXXX

function glr-preq() {
    # if $HOME/gitlab-runner exists
    if [ -f "$HOME/gitlab-runner" ]; then
        cp $HOME/gitlab-runner .
    else
        curl -L --output gitlab-runner https://gitlab-runner-downloads.s3.amazonaws.com/latest/binaries/gitlab-runner-linux-amd64
    fi

    chmod +x gitlab-runner

}

function glr-new() {
    local runner_url=$1; shift
    local token=$1; shift
    local executor=${1:-shell}; shift
    local shell=${1:-bash}; shift

    pushd $RUNNER_TOP_DIR

    # check if gitlab-runner exists
    if [ ! -f gitlab-runner ]; then
        glr-preq
    fi

    ./gitlab-runner register --url $runner_url --token $token --executor $executor --shell $shell

    popd
}

function new() {
    local token=$1; shift
    if [ -z "$token" ]; then
        echo "Please pass the token to this script" 1>&2
        exit -1
    fi
    glr-new $RUNNER_URL $token

}

#############################################
# Create a new gitlab runner (glr)
#############################################

function glr-start() {
    local glr=gitlab-runner

    pushd $RUNNER_TOP_DIR

    apptainer instance start ~/github-runner.sif ${glr}
    apptainer run instance://${glr} bash -c "./gitlab-runner run -c ./config.toml &"

    popd
}

function start() {
    glr-start
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
