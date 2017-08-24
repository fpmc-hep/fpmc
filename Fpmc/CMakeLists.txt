add_subdirectory(External)
include_directories(External ${PROJECT_SOURCE_DIR}/Herwig)

add_library(FpmcExternals SHARED $<TARGET_OBJECTS:CHIDe> $<TARGET_OBJECTS:KMR_lumi> $<TARGET_OBJECTS:KMR_lumi_2> $<TARGET_OBJECTS:excl_aaaa> $<TARGET_OBJECTS:softc> $<TARGET_OBJECTS:miscellaneous>)
target_link_libraries(FpmcExternals comphep)

add_library(FpmcCore SHARED src/Fpmc.cc src/fpmc.f src/FpmcParameters.cc ${PROJECT_SOURCE_DIR}/Examples/fpmc_welcome.f)
set_target_properties(FpmcCore PROPERTIES LINKER_LANGUAGE CXX POSITION_INDEPENDENT_CODE ON)

if(HEPMC_LIB)
  add_library(FpmcHepMC SHARED src/HepMCWrapper.cc)
endif()
