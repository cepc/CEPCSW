#!/bin/bash

#source /besfs/groups/higgs/Software/v01-17-05/root/5.34.07/bin/thisroot.sh
#
#export INSDIR="/afs/ihep.ac.cn/soft/common/gcc/whizard227/install"
#export GCCPATH=${INSDIR}
#export PATH_TO_GFORTRAN=$GCCPATH/bin
#alias gcc=$GCCPATH/bin/gcc
#alias g++=$GCCPATH/bin/g++
#export LD_GFORTRAN=${GCCDIR}/lib
#
#export LCIO=/afs/ihep.ac.cn/soft/common/gcc/v01-17-05/lcio/v02-04-03/
#export DRUIDDIR=/afs/ihep.ac.cn/users/m/manqi/Software/Druid_SL6
#export PATH=$DRUIDDIR/bin:/afs/ihep.ac.cn/users/m/manqi/ArborV3Deve/SL6/Marlin/v01-05/bin:/afs/ihep.ac.cn/users/m/manqi/Software/ilcsoft/v01-16/CMake/2.8.5/bin:/usr/local/bin:/usr/bin:/bin:/usr/local/sbin:/usr/sbin:/sbin:$ROOTSYS/bin:/usr/X11R6/bin:$LCIO/bin:/opt/glite/bin/:/opt/lcg/bin/:/afs/ihep.ac.cn/soft/common/gcc/v01-17-05/CMake/2.8.5/bin
#export LD_LIBRARY_PATH=$ROOTSYS/lib:/home/manqi/Softwares/ilcsoft/v01-08-01/CLHEP/2.0.4.2/lib:$LCIO/lib:/usr/lib:/lib:${GCCPATH}/lib64:${GCCPATH}/lib:

unset MARLIN_DLL

#export MARLIN_DLL=/afs/ihep.ac.cn/soft/common/gcc/v01-17-05/PandoraAnalysis/v00-05/lib/libPandoraAnalysis.so:/afs/ihep.ac.cn/soft/common/gcc/v01-17-05/MarlinPandora/v00-11/lib/libMarlinPandora.so:/afs/ihep.ac.cn/users/y/yudan/yudan/Arbor_v3.3/ArborMag/lib/libRanger.so:/afs/ihep.ac.cn/soft/common/gcc/v01-17-05/MarlinTPC/v00-16/lib/libMarlinTPC.so:/afs/ihep.ac.cn/soft/common/gcc/v01-17-05/MarlinTrkProcessors/v01-10/lib/libMarlinTrkProcessors.so:/afs/ihep.ac.cn/soft/common/gcc/v01-17-05/Clupatra/v00-10/lib/libClupatra.so:/afs/ihep.ac.cn/soft/common/gcc/v01-17-05/ForwardTracking/v01-07/lib/libForwardTracking.so:/afs/ihep.ac.cn/soft/common/gcc/v01-17-05/MarlinReco/v01-09/lib/libMarlinReco.so:/afs/ihep.ac.cn/soft/common/gcc/v01-17-05/LCFIPlus/v00-05-02/lib/libLCFIPlus.so:/afs/ihep.ac.cn/soft/common/gcc/v01-17-05/Overlay/v00-13/lib/libOverlay.so
export MARLIN_DLL=/cvmfs/cepc.ihep.ac.cn/software/cepcsoft/x86_64-sl6-gcc49/cepcsoft/0.1.0-rc9/Reconstruction/Tracking/MarlinTrkProcessors/01-10/lib/libMarlinTrkProcessors.so:/cvmfs/cepc.ihep.ac.cn/software/cepcsoft/x86_64-sl6-gcc49/cepcsoft/0.1.0-rc9/Reconstruction/Tracking/ForwardTracking/01-07/lib/libForwardTracking.so:/cvmfs/cepc.ihep.ac.cn/software/cepcsoft/x86_64-sl6-gcc49/cepcsoft/0.1.0-rc9/Reconstruction/Tracking/Clupatra/00-10/lib/libClupatra.so:/cvmfs/cepc.ihep.ac.cn/software/cepcsoft/x86_64-sl6-gcc49/cepcsoft/0.1.0-rc9/Analysis/SubDetectorStudy/TPC/MarlinTPC/00-16/lib/libMarlinTPC.so:/cvmfs/cepc.ihep.ac.cn/software/cepcsoft/x86_64-sl6-gcc49/cepcsoft/0.1.0-rc9/Reconstruction/HighLevelObjectFinding/Jets/LCFIPlus/00-05-02/lib/libLCFIPlus.so:/cvmfs/cepc.ihep.ac.cn/software/cepcsoft/x86_64-sl6-gcc49/cepcsoft/0.1.0-rc9/Reconstruction/PFA/Pandora/MarlinPandora/00-11/lib/libMarlinPandora.so:/cvmfs/cepc.ihep.ac.cn/software/cepcsoft/x86_64-sl6-gcc49/cepcsoft/0.1.0-rc9/Reconstruction/HighLevelObjectFinding/Jets/FastJetClustering/00-02/lib/libFastJetClustering.so:/cvmfs/cepc.ihep.ac.cn/software/cepcsoft/x86_64-sl6-gcc49/cepcsoft/0.1.0-rc9/Reconstruction/PFA/Garlic/2.10.1/lib/libGarlic.so:/cvmfs/cepc.ihep.ac.cn/software/cepcsoft/x86_64-sl6-gcc49/cepcsoft/0.1.0-rc9/Reconstruction/HighLevelObjectFinding/Jets/LCFIVertex/00-06-01/lib/libLCFIVertex.so:/cvmfs/cepc.ihep.ac.cn/software/cepcsoft/x86_64-sl6-gcc49/cepcsoft/0.1.0-rc9/Reconstruction/Digitization/MarlinReco/01-09/lib/libMarlinReco.so:/cvmfs/cepc.ihep.ac.cn/software/cepcsoft/x86_64-sl6-gcc49/cepcsoft/0.1.0-rc9/EventDisplay/CEDViewer/01-07-02/lib/libCEDViewer.so:/cvmfs/cepc.ihep.ac.cn/software/cepcsoft/x86_64-sl6-gcc49/cepcsoft/0.1.0-rc9/Reconstruction/PFA/Pandora/PandoraAnalysis/00-05/lib/libPandoraAnalysis.so:/cvmfs/cepc.ihep.ac.cn/software/cepcsoft/x86_64-sl6-gcc49/cepcsoft/0.1.0-rc9/Analysis/Overlay/00-13/lib/libOverlay.so:/afs/ihep.ac.cn/users/y/yudan/yudan/Arbor_v3.3/ArborMag/lib/libRanger.so:/cvmfs/cepc.ihep.ac.cn/software/cepcsoft/x86_64-sl6-gcc49/cepcsoft/0.1.0-rc9/Reconstruction/HighLevelObjectFinding/Jets/MarlinFastJet/00-02/lib/libMarlinFastJet.so:/cvmfs/cepc.ihep.ac.cn/software/cepcsoft/x86_64-sl6-gcc49/cepcsoft/0.1.0-rc9/Framework/LCTuple/01-03/lib/libLCTuple.so


