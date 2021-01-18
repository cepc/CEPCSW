# This variable will be used by GaudiToolbox.cmake to generate the cepcswenv.sh

set(RUN_SCRIPT_EXTRA_COMMANDS "${RUN_SCRIPT_EXTRA_COMMANDS}
export CEPCSW_ROOT=${CMAKE_SOURCE_DIR}

export DETCEPCV4ROOT=${CMAKE_SOURCE_DIR}/Detector/DetCEPCv4
export DETCRDROOT=${CMAKE_SOURCE_DIR}/Detector/DetCRD
export DETDRIFTCHAMBERROOT=${CMAKE_SOURCE_DIR}/Detector/DetDriftChamber
")
