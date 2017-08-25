include_directories(${PROJECT_SOURCE_DIR}/Herwig)
include_directories(${PROJECT_SOURCE_DIR}/External)
include_directories(interface)

set(fpmc_src src/fpmc.f src/fpmc_welcome.f src/Fpmc.cc src/FpmcParameters.cc)
if(CERNLIB)
  list(APPEND fpmc_src src/ffcard.f)
endif()

add_library(FPMCCore OBJECT ${fpmc_src})
set_target_properties(FPMCCore PROPERTIES LINKER_LANGUAGE CXX POSITION_INDEPENDENT_CODE ON)
