add_subdirectory(External)

include_directories(External)
include_directories(${PROJECT_SOURCE_DIR}/Herwig)
include_directories(interface)

set(fpmc_sources src/fpmc.f src/fpmc_welcome.f src/Fpmc.cc src/FpmcParameters.cc)
if(CERNLIB)
  list(APPEND fpmc_sources src/ffkey.f)
endif()
if(HEPMC_LIB)
  list(APPEND fpmc_sources src/HepMCWrapper.cc)
endif()

#----- core of the generator

add_library(FpmcCore SHARED ${fpmc_sources})
if(HEPMC_LIB)
  target_link_libraries(FpmcCore ${HEPMC_LIB} ${HEPMC_FIO_LIB})
endif()
set_target_properties(FpmcCore PROPERTIES LINKER_LANGUAGE CXX POSITION_INDEPENDENT_CODE ON)

#----- installation procedure

install(TARGETS FpmcCore DESTINATION lib)

