cmake_minimum_required(VERSION 2.6 FATAL_ERROR)
project(FPMC)
set(PROJECT_VERSION 1)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wno-long-long -pedantic -g")
set(MODULES Fpmc Herwig External HepMCWrapper Examples)

#----- enable Fortran

enable_language(Fortran)

#----- external dependencies

find_package(ROOT)
find_package(GSL)
find_library(HEPMC_LIB HepMC)
find_library(LHAPDF_LIB LHAPDF)
find_library(CLHEP_LIB CLHEP)
if(LHAPDF_LIB)
  message(STATUS "LHAPDF found at " ${LHAPDF_LIB})
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

#----- link the HepMC wrapper into an export library

add_library(FPMCHepMC SHARED $<TARGET_OBJECTS:FPMCCore>
                             $<TARGET_OBJECTS:HerwigCore>
                             #$<TARGET_OBJECTS:ff>
                             $<TARGET_OBJECTS:Export>)
target_link_libraries(FPMCHepMC ${CLHEP_LIB})

#----- test executable

set(LIBRARIES FPMC Externals comphep ${GSL_LIBRARIES} ${LHAPDF_LIB})
include_directories(Examples)
add_executable(fpmc Examples/fpmc_main.f)
target_link_libraries(fpmc ${LIBRARIES})

if (HEPMC_LIB)
  message(STATUS "HepMC found at " ${HEPMC_LIB})
  if (ROOT_LIBRARIES)
    message(STATUS "ROOT found in " ${ROOT_INCLUDE_DIRS})
    include_directories(HepMCWrapper Examples)
    add_executable(fpmc-hepmc HepMCWrapper/main.cc)
    target_link_libraries(fpmc-hepmc ${LIBRARIES} FPMCHepMC ${HEPMC_LIB} ${ROOT_LIBRARIES})
    set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -DHEPMC_NO_DEPRECATED")
  endif()
endif()
