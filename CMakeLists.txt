cmake_minimum_required(VERSION 2.6 FATAL_ERROR)
project(FPMC)
set(PROJECT_VERSION 1)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wno-long-long -Wno-deprecated-declarations -pedantic -std=c++11 -g")

set(MODULES Fpmc Herwig External HepMCWrapper Examples Tests)

#----- enable Fortran

enable_language(Fortran)

#----- external dependencies

if($ENV{HOSTNAME} MATCHES "^lxplus[0-9]+.cern.ch")
  set(BASE_DIR "/cvmfs/sft.cern.ch/lcg/external")
  set(CERNLIB_DIR "${BASE_DIR}/cernlib/2006a/x86_64-slc6-gcc47-opt")
  set(GSL_DIR "${BASE_DIR}/GSL/1.14/x86_64-slc5-gcc44-opt")
  set(CLHEP_DIR "${BASE_DIR}/clhep/2.2.0.4/x86_64-slc6-gcc48-opt")
  set(HEPMC_DIR "${BASE_DIR}/HepMC/2.06.08/x86_64-slc6-gcc48-opt")
  set(LHAPDF_DIR "${BASE_DIR}/MCGenerators/lhapdf/5.8.9/x86_64-slc6-gcc46-opt")

  message(STATUS "Compiling on LXPLUS. Do not forget to source the environment variables!")

  find_library(CERNLIB_LIB NAMES cernlib HINTS "${CERNLIB_DIR}/lib")
  find_library(CERNLIB_MATHLIB_LIB NAMES mathlib HINTS "${CERNLIB_DIR}/lib")
  find_library(CERNLIB_PACKLIB_LIB NAMES packlib HINTS "${CERNLIB_DIR}/lib")
  find_library(CERNLIB_PAWLIB_LIB NAMES pawlib HINTS "${CERNLIB_DIR}/lib")
  message(STATUS "PAW found in ${CERNLIB_PAWLIB_LIB}")

  find_library(GSL_LIB gsl HINTS "${GSL_DIR}/lib")
  find_library(GSL_CBLAS_LIB gslcblas HINTS "${GSL_DIR}/lib")
  #--- searching for CLHEP
  find_library(CLHEP_LIB CLHEP HINTS "${CLHEP_DIR}/lib")
  find_path(CLHEP_INCLUDE CLHEP HINTS "${CLHEP_DIR}/include")
  #--- searching for LHAPDF
  find_library(LHAPDF_LIB LHAPDF HINTS "${LHAPDF_DIR}/lib")
  #--- searching for HepMC
  find_library(HEPMC_LIB HepMC HINTS "${HEPMC_DIR}/lib")
  find_library(HEPMC_FIO_LIB HepMCfio HINTS "${HEPMC_DIR}/lib")
  find_path(HEPMC_INCLUDE HepMC HINTS "${HEPMC_DIR}/include")
else()
  find_library(GSL_LIB gsl)
  find_library(GSL_CBLAS_LIB gslcblas)
  find_library(CLHEP_LIB CLHEP)
  find_path(CLHEP_INCLUDE CLHEP)
  find_library(LHAPDF_LIB LHAPDF)
  find_library(HEPMC_LIB HepMC)
  find_path(HEPMC_INCLUDE HepMC)
endif()

include_directories(${CLHEP_INCLUDE})

#----- all useful libraries

set(LIBRARIES FPMC ${GSL_LIB} ${CLHEP_LIB} ${GSL_CBLAS_LIB} ${LHAPDF_LIB})

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

#----- link the HepMC wrapper into a dedicated library

if(HEPMC_LIB)
  add_library(HepMCWrapper SHARED $<TARGET_OBJECTS:Wrapper>)
  target_link_libraries(HepMCWrapper ${CLHEP_LIB})
endif()

#----- test executable

if(CERNLIB_LIB)
  include_directories(Examples)
  add_executable(fpmc-paw Examples/fpmc_main.f Examples/fpmc_lhe.f External/ntuple.f External/reco_ok/reco_ok.f External/reco_ok/pxcone_mod_new.f Examples/ffcard.f Examples/fpmc_welcome.f)
  target_link_libraries(fpmc-paw ${LIBRARIES} ${CERNLIB_LIB} ${CERNLIB_MATHLIB_LIB} ${CERNLIB_PAWLIB_LIB} ${CERNLIB_PACKLIB_LIB})
endif()

if(HEPMC_LIB)
  include_directories(${HEPMC_INCLUDE})
  message(STATUS "HepMC found at " ${HEPMC_LIB})
  include_directories(HepMCWrapper Examples)
  add_executable(fpmc-hepmc HepMCWrapper/main.cc Examples/dummy_hwaend.f)
  target_link_libraries(fpmc-hepmc ${LIBRARIES} HepMCWrapper ${HEPMC_LIB} ${HEPMC_FIO_LIB})
endif()

#----- copies

file(COPY ${PROJECT_SOURCE_DIR}/External/CHIDe/Data/ DESTINATION External/CHIDe/Data/)
