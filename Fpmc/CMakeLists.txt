add_subdirectory(External)
include_directories(External ${PROJECT_SOURCE_DIR}/Herwig)

#----- core of the generator

add_library(FpmcCore SHARED src/Fpmc.cc src/fpmc.f src/FpmcParameters.cc ${PROJECT_SOURCE_DIR}/Examples/fpmc_welcome.f)
set_target_properties(FpmcCore PROPERTIES LINKER_LANGUAGE CXX POSITION_INDEPENDENT_CODE ON)

#----- installation procedure

install(TARGETS FpmcCore DESTINATION lib)

if(HEPMC_LIB)
  add_library(FpmcHepMC SHARED src/HepMCWrapper.cc)
  install(TARGETS FpmcHepMC DESTINATION lib)
endif()
