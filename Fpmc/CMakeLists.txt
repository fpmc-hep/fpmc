include_directories(${PROJECT_SOURCE_DIR}/External ${PROJECT_SOURCE_DIR}/Herwig)

add_library(FPMCCore OBJECT fpmc.f)
set_target_properties(FPMCCore PROPERTIES LINKER_LANGUAGE CXX POSITION_INDEPENDENT_CODE ON)
