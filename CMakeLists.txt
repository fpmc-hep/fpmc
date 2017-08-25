cmake_minimum_required(VERSION 2.6 FATAL_ERROR)
project(FPMC)
set(PROJECT_VERSION 1)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wno-long-long -Wno-deprecated-declarations -pedantic -std=c++11 -g")
set(CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}/cmake")

set(MODULES Fpmc Herwig External)

#----- enable Fortran

enable_language(Fortran)

#----- external dependencies

include(UseEnvironment)
include_directories(${CLHEP_INCLUDE})

#----- build all the intermediate objects

foreach(_mod ${MODULES})
  include_directories(${PROJECT_SOURCE_DIR}/${_mod})
  add_subdirectory(${_mod})
endforeach()

#----- link everything in a single shared FPMC library (FPMC and HERWIG altogether)

add_library(FPMCBare SHARED EXCLUDE_FROM_ALL $<TARGET_OBJECTS:FPMCCore> $<TARGET_OBJECTS:HerwigCore>)
set_target_properties(FPMCBare PROPERTIES LINKER_LANGUAGE CXX)

add_library(FPMCDeps SHARED EXCLUDE_FROM_ALL $<TARGET_OBJECTS:CHIDe> $<TARGET_OBJECTS:KMR_lumi> $<TARGET_OBJECTS:KMR_lumi_2> $<TARGET_OBJECTS:excl_aaaa> $<TARGET_OBJECTS:softc> $<TARGET_OBJECTS:miscellaneous>)
target_link_libraries(FPMCDeps comphep)
set_target_properties(FPMCDeps PROPERTIES LINKER_LANGUAGE CXX)

add_library(FPMC SHARED $<TARGET_OBJECTS:FPMCCore> $<TARGET_OBJECTS:HerwigCore> $<TARGET_OBJECTS:CHIDe> $<TARGET_OBJECTS:KMR_lumi> $<TARGET_OBJECTS:KMR_lumi_2> $<TARGET_OBJECTS:excl_aaaa> $<TARGET_OBJECTS:softc> $<TARGET_OBJECTS:miscellaneous>)
target_link_libraries(FPMC comphep)
set_target_properties(FPMC PROPERTIES LINKER_LANGUAGE CXX)

#----- all useful libraries

set(LIBRARIES FPMC ${GSL_LIB} ${GSL_CBLAS_LIB} ${LHAPDF_LIB} ${CLHEP_LIB})

#----- test executable

if(CERNLIB_LIB)
  include_directories(Examples)
  add_executable(fpmc-paw Examples/fpmc_main.f Examples/fpmc_lhe.f External/ntuple.f External/reco_ok/reco_ok.f External/reco_ok/pxcone_mod_new.f)
  target_link_libraries(fpmc-paw ${LIBRARIES} ${CERNLIB_LIB} ${CERNLIB_MATHLIB_LIB} ${CERNLIB_PAWLIB_LIB} ${CERNLIB_PACKLIB_LIB})
endif()

if(HEPMC_LIB)
  include_directories(${HEPMC_INCLUDE})
  message(STATUS "HepMC found at " ${HEPMC_LIB})
  include_directories(Examples)
  add_executable(fpmc-hepmc Examples/fpmc_hepmc.cc Examples/dummy_hwaend.f)
  target_link_libraries(fpmc-hepmc ${LIBRARIES} ${HEPMC_LIB} ${HEPMC_FIO_LIB})
endif()

#----- copies

file(COPY ${PROJECT_SOURCE_DIR}/External/CHIDe/Data/ DESTINATION External/CHIDe/Data/)
