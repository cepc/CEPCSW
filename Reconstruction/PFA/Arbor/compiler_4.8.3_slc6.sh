

export INSDIR="/afs/ihep.ac.cn/soft/common/gcc/whizard227/install"
export PATH=/afs/ihep.ac.cn/soft/common/gcc/v01-17-05/CMake/2.8.5/bin:$PATH

export GCCPATH=${INSDIR}
export PATH_TO_GFORTRAN=$GCCPATH/bin
alias gcc=$GCCPATH/bin/gcc
alias g++=$GCCPATH/bin/g++

export LD_GFORTRAN=${GCCDIR}/lib
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:${GCCPATH}/lib64:${GCCPATH}/lib:

