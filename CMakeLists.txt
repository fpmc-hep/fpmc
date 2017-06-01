cmake_minimum_required(VERSION 2.6 FATAL_ERROR)
project(FPMC)
set(PROJECT_VERSION 1)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wno-long-long -pedantic -std=c++11 -g")
set(MODULES Fpmc Herwig External HepMCWrapper Examples)

#----- enable Fortran

enable_language(Fortran)

#----- external dependencies

if($ENV{HOSTNAME} MATCHES "^lxplus[0-9]+.cern.ch")
  set(BASE_DIR "/cvmfs/sft.cern.ch/lcg/external")
  set(GSL_DIR "${BASE_DIR}/GSL/1.14/x86_64-slc5-gcc44-opt")
  set(CLHEP_DIR "${BASE_DIR}/clhep/2.2.0.4/x86_64-slc6-gcc48-opt")
  set(LHAPDF_DIR "${BASE_DIR}/MCGenerators/lhapdf/5.8.9/x86_64-slc6-gcc46-opt")

  message(STATUS "Compiling on LXPLUS")

  find_library(GSL_LIB gsl HINTS "${GSL_DIR}/lib")
  find_library(GSL_CBLAS_LIB gslcblas HINTS "${GSL_DIR}/lib")
  #--- searching for CLHEP
  find_library(CLHEP_LIB CLHEP HINTS "${CLHEP_DIR}/lib")
  find_path(CLHEP_INCLUDE CLHEP HINTS "${CLHEP_DIR}/include")
  #--- searching for LHAPDF
  find_library(LHAPDF_LIB LHAPDF HINTS "${LHAPDF_DIR}/lib")
  #--- searching for HepMC
  find_library(HEPMC_LIB HepMC HINTS "${BASE_DIR}/HepMC/2.06.08/x86_64-slc6-gcc48-opt/lib")
  find_library(HEPMC_FIO_LIB HepMCfio HINTS "${BASE_DIR}/HepMC/2.06.08/x86_64-slc6-gcc48-opt/lib")
  find_path(HEPMC_INCLUDE HepMC HINTS "${BASE_DIR}/HepMC/2.06.08/x86_64-slc6-gcc48-opt/include")
  #find_library(HEPMC_LIB HepMC HINTS "${HEPMC_DIR}/lib")
  #find_path(HEPMC_INCLUDE HepMC HINTS "${HEPMC_DIR}/include")
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
find_package(ROOT)
#if(LHAPDF_LIB)
#  message(STATUS "LHAPDF found at " ${LHAPDF_LIB})
#else()
#  message(FATAL_ERROR "Did not find LHAPDF in the common paths!")
#endif()

#----- build all the intermediate objects

foreach(_mod ${MODULES})
  include_directories(${PROJECT_SOURCE_DIR}/${_mod})
  add_subdirectory(${_mod})
endforeach()

#----- link everything in a single shared FPMC library

add_library(FPMC SHARED $<TARGET_OBJECTS:FPMCCore>
                        $<TARGET_OBJECTS:HerwigCore>)
set_target_properties(FPMC PROPERTIES LINKER_LANGUAGE CXX)

#----- link the HepMC wrapper into an export library

if(HEPMC_LIB)
  add_library(FPMCHepMC SHARED $<TARGET_OBJECTS:FPMCCore>
                               $<TARGET_OBJECTS:HerwigCore>
                               #$<TARGET_OBJECTS:ff>
                               $<TARGET_OBJECTS:Export>)
  target_link_libraries(FPMCHepMC ${CLHEP_LIB})
endif()

#----- test executable

set(LIBRARIES FPMC Externals comphep ${GSL_LIB} ${GSL_CBLAS_LIB} ${LHAPDF_LIB})
include_directories(Examples)
add_executable(fpmc Examples/fpmc_main.f)
target_link_libraries(fpmc ${LIBRARIES})

if(HEPMC_LIB)
  include_directories(${HEPMC_INCLUDE})
  message(STATUS "HepMC found at " ${HEPMC_LIB})
  if(ROOT_LIBRARIES)
    message(STATUS "ROOT found in " ${ROOT_INCLUDE_DIRS})
    include_directories(HepMCWrapper Examples)
    add_executable(fpmc-hepmc HepMCWrapper/main.cc)
    #target_link_libraries(fpmc-hepmc ${LIBRARIES} FPMCHepMC ${HEPMC_LIB} ${HEPMC_FIO_LIB} ${ROOT_LIBRARIES})
    target_link_libraries(fpmc-hepmc ${LIBRARIES} FPMCHepMC ${HEPMC_LIB} ${HEPMC_FIO_LIB})
    #set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -DHEPMC_NO_DEPRECATED")
  endif()
endif()
