cmake_minimum_required(VERSION 2.6 FATAL_ERROR)
project(FPMC)
set(PROJECT_VERSION 1)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wno-long-long -Wno-deprecated-declarations -pedantic -std=c++11 -lgfortran -g")
set(CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}/cmake")

set(MODULES Fpmc Herwig)

#----- enable Fortran

enable_language(Fortran)

#----- external dependencies

include(UseEnvironment)
include_directories(${PROJECT_SOURCE_DIR})
include_directories(${CLHEP_INCLUDE})
include_directories(${HEPMC_INCLUDE})
message(STATUS "HepMC found at " ${HEPMC_LIB})

#----- build all the intermediate objects

foreach(_mod ${MODULES})
  include_directories(${PROJECT_SOURCE_DIR}/${_mod})
  add_subdirectory(${_mod})
endforeach()

#----- link everything in a single shared FPMC library (FPMC and HERWIG altogether)

add_library(Fpmc SHARED $<TARGET_OBJECTS:FPMCCore> $<TARGET_OBJECTS:HerwigCore> $<TARGET_OBJECTS:CHIDe> $<TARGET_OBJECTS:KMR_lumi> $<TARGET_OBJECTS:KMR_lumi_2> $<TARGET_OBJECTS:excl_aaaa> $<TARGET_OBJECTS:softc> $<TARGET_OBJECTS:miscellaneous>)
target_link_libraries(Fpmc comphep)
set_target_properties(Fpmc PROPERTIES LINKER_LANGUAGE CXX)

#----- all useful libraries

set(LIBRARIES Fpmc ${GSL_LIB} ${GSL_CBLAS_LIB} ${LHAPDF_LIB} ${CLHEP_LIB} ${HEPMC_LIB} ${HEPMC_FIO_LIB})

#----- test executable

if(CERNLIB_LIB)
  include_directories(Examples)
  add_executable(fpmc-paw Examples/fpmc_main.f Examples/fpmc_lhe.f Examples/ntuple.f Examples/reco_ok/reco_ok.f Examples/reco_ok/pxcone_mod_new.f)
  target_link_libraries(fpmc-paw ${LIBRARIES} ${CERNLIB_LIB} ${CERNLIB_MATHLIB_LIB} ${CERNLIB_PAWLIB_LIB} ${CERNLIB_PACKLIB_LIB})
endif()

include_directories(${HEPMC_INCLUDE})
include_directories(Examples)

add_executable(fpmc Examples/fpmc_bare.cc Examples/dummy_hwaend.f)
target_link_libraries(fpmc ${LIBRARIES})

add_executable(fpmc-hepmc Examples/fpmc_hepmc.cc Examples/dummy_hwaend.f)
target_link_libraries(fpmc-hepmc ${LIBRARIES})

#----- installation procedure

install(TARGETS Fpmc DESTINATION lib)
install(DIRECTORY Fpmc DESTINATION include FILES_MATCHING PATTERN "*.h")
install(DIRECTORY Fpmc DESTINATION include FILES_MATCHING PATTERN "*.inc")

#----- copies

file(COPY ${PROJECT_SOURCE_DIR}/Datacards/ DESTINATION Datacards/)
file(COPY ${PROJECT_SOURCE_DIR}/Fpmc/External/CHIDe/Data/ DESTINATION Fpmc/External/CHIDe/Data/)

