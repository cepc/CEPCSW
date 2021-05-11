#!/bin/bash
# This is wrapper to run the build.sh on CI

echo "LCG_RELEASE: ${LCG_RELEASE}"
echo "CEPCSW_BLDTOOL: ${CEPCSW_BLDTOOL}"
buildpid=
logfile=mylog.txt

if [ "$LCG_RELEASE" = "KEY4HEP_STACK" ]; then
    logfile=mylog-k4.sh
    source /cvmfs/sw.hsf.org/key4hep/setup.sh
    ./build-k4.sh >& ${logfile} &
    buildpid=$!
else
    source setup.sh
    ./build.sh >& ${logfile} &
    buildpid=$!
fi

while ps -p $buildpid 2>/dev/null ; do
    sleep 60
done &
echoer=$!

trap 'kill $echoer' 0

wait $buildpid
statuspid=$?

tail -n100 ${logfile}

exit $statuspid
