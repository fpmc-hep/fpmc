cmake_minimum_required(VERSION 2.6 FATAL_ERROR)
project(FPMC)
set(PROJECT_VERSION 1)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wno-long-long -Wno-deprecated-declarations -pedantic -lgfortran -std=c++11 -g")
#----- small fix for gcc-gfortran versions >= 8.1
#----- (while CHIDe/KMR codes are not fixed)
if(CMAKE_CXX_COMPILER_ID MATCHES "Clang")
else()
  if(CMAKE_CXX_COMPILER_VERSION VERSION_GREATER 8.1)
    set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -std=legacy")
  endif()
endif()
set(CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}/cmake")

set(MODULES Fpmc Herwig)

#----- enable Fortran

enable_language(Fortran)

#----- external dependencies

include(UseEnvironment)
include_directories(${PROJECT_SOURCE_DIR})
if(CLHEP_LIB)
  message(STATUS "CLHEP found at ${CLHEP_LIB}")
  include_directories(${CLHEP_INCLUDE})
else()
  message(FATAL_ERROR "CLHEP is required to build ${PROJECT_NAME}")
endif()
if(HEPMC_LIB)
  message(STATUS "HepMC found at ${HEPMC_LIB}")
  include_directories(${HEPMC_INCLUDE})
else()
  message(FATAL_ERROR "HepMC is required to build ${PROJECT_NAME}")
endif()

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

include_directories(Fpmc/interface)
include_directories(Fpmc/External)
include_directories(Examples)
include_directories(${HEPMC_INCLUDE})

if(CERNLIB_PAWLIB_LIB)
  add_executable(fpmc-lhe Examples/fpmc_main_lhe.f Examples/fpmc_lhe.f Examples/fpmc_load_params.f Examples/fpmc_var_ini.f)
  target_link_libraries(fpmc-lhe ${LIBRARIES} ${CERNLIB_PACKLIB_LIB})

  add_executable(fpmc-paw Examples/fpmc_main.f Examples/fpmc_load_params.f Examples/fpmc_var_ini.f Examples/fpmc_lhe.f Examples/ntuple.f)
  target_link_libraries(fpmc-paw ${LIBRARIES} ${CERNLIB_MATHLIB_LIB} ${CERNLIB_PAWLIB_LIB} ${CERNLIB_PACKLIB_LIB})
else()
  add_executable(fpmc-lhe Examples/fpmc_main_lhe.f Examples/fpmc_lhe.f Examples/fpmc_load_params.f Examples/fpmc_var_ini.f Examples/fpmc_legacy_ffcard.cc)
  target_link_libraries(fpmc-lhe ${LIBRARIES})
endif()

add_executable(fpmc Examples/fpmc_bare.cc Examples/dummy_hwaend.f)
target_link_libraries(fpmc ${LIBRARIES})

add_executable(fpmc-hepmc Examples/fpmc_hepmc.cc Examples/dummy_hwaend.f)
target_link_libraries(fpmc-hepmc ${LIBRARIES})

#find_package(ROOT QUIET)
#if(ROOT_FOUND)
#  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++14")
#  add_executable(fpmc-root Examples/fpmc_root.cc Examples/dummy_hwaend.f)
#  include_directories(${ROOT_INCLUDE_DIRS})
#  target_link_libraries(fpmc-root ${LIBRARIES} ${ROOT_LIBRARIES})
#endif()

#----- installation procedure

install(TARGETS Fpmc DESTINATION lib)
install(DIRECTORY Fpmc DESTINATION include FILES_MATCHING PATTERN "*.h")
install(DIRECTORY Fpmc DESTINATION include FILES_MATCHING PATTERN "*.inc")

#----- copies

file(COPY ${PROJECT_SOURCE_DIR}/Datacards/ DESTINATION Datacards/)
file(COPY ${PROJECT_SOURCE_DIR}/Fpmc/External/CHIDe/Data/ DESTINATION External/CHIDe/Data/)
file(COPY ${PROJECT_SOURCE_DIR}/Fpmc/External/excl_aaaa/MG_standalone_eft_ttbar/Cards DESTINATION External/excl_aaaa/MG_standalone_eft_ttbar/)
