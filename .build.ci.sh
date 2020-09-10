#!/bin/bash
# This is wrapper to run the build.sh on CI

source setup.sh
./build.sh >& mylog.sh &
buildpid=$!

while ps -p $buildpid 2>/dev/null ; do
    sleep 60
done

tail -n100 mylog.sh

