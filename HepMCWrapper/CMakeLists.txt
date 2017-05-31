include_directories(${PROJECT_SOURCE_DIR}/Examples)
add_library(Export OBJECT Fpmc.cc ${PROJECT_SOURCE_DIR}/Examples/fpmc_welcome.f)
set_target_properties(Export PROPERTIES LINKER_LANGUAGE CXX POSITION_INDEPENDENT_CODE ON)
