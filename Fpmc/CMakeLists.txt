include_directories(${PROJECT_SOURCE_DIR}/Herwig)
include_directories(${PROJECT_SOURCE_DIR}/External)
include_directories(interface)

add_library(FPMCCore OBJECT src/fpmc.f src/fpmc_welcome.f src/Fpmc.cc src/FpmcParameters.cc)
set_target_properties(FPMCCore PROPERTIES LINKER_LANGUAGE CXX POSITION_INDEPENDENT_CODE ON)
