gaudi_subdir(DataHelper v0r0)

find_package(EDM4HEP REQUIRED)
set (gsl_include "/cvmfs/cepc.ihep.ac.cn/software/cepcsoft/x86_64-sl6-gcc49/external/GSL/1.14/install/include")
set (gsl_lib1 "/cvmfs/cepc.ihep.ac.cn/software/cepcsoft/x86_64-sl6-gcc49/external/GSL/1.14/install/lib/libgsl.so")
set (gsl_lib2 "/cvmfs/cepc.ihep.ac.cn/software/cepcsoft/x86_64-sl6-gcc49/external/GSL/1.14/install/lib/libgslcblas.so")
INCLUDE_DIRECTORIES(${gsl_include})

gaudi_depends_on_subdirs()

set(DataHelperLib_srcs src/*.cc src/*.cpp)

#gaudi_install_headers(DataHelper)

gaudi_add_library(DataHelperLib ${DataHelperLib_srcs}
    PUBLIC_HEADERS DataHelper
    LINK_LIBRARIES EDM4HEP::edm4hep EDM4HEP::edm4hepDict ${gsl_lib1} ${gsl_lib2}
)
