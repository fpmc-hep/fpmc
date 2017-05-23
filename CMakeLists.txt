cmake_minimum_required(VERSION 2.6 FATAL_ERROR)
project(FPMC)
set(PROJECT_VERSION 1)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wno-long-long -pedantic -g -fPIC")

set(MODULES Fpmc Herwig External)

#----- enable Fortran

enable_language(Fortran)

#----- external dependencies

find_package(HepMC)
find_package(GSL)
#find_package(LHAPDF)
find_library(LHAPDF_LIB LHAPDF)
if(LHAPDF_LIB)
  message(STATUS "LHAPDF found in " ${LHAPDF_LIB})
else()
  message(FATAL_ERROR "Did not find LHAPDF in the common paths!")
endif()

#----- build all the intermediate objects

foreach(_mod ${MODULES})
  include_directories(${PROJECT_SOURCE_DIR}/${_mod})
  add_subdirectory(${_mod})
endforeach()

#----- link everything in a single shared FPMC library

add_library(FPMC SHARED $<TARGET_OBJECTS:FPMCCore>
                        $<TARGET_OBJECTS:HerwigCore>)
set_target_properties(FPMC PROPERTIES LINKER_LANGUAGE CXX)

#----- test executable

add_executable(fpmc fpmc.cpp)
target_link_libraries(fpmc FPMC Externals ${GSL_LIBRARIES} ${LHAPDF_LIB})

if (HEPMC_LIBRARIES)
  message(STATUS "HepMC found in " ${HEPMC_INCLUDE_DIR})
  add_executable(fpmc-hepmc fpmc-hepmc.cpp)
  target_link_libraries(fpmc-hepmc FPMC Externals ${HEPMC_LIBRARIES} ${GSL_LIBRARIES} ${LHAPDF_LIB})
endif()
