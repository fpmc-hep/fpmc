#include_directories(${PROJECT_SOURCE_DIR}/External ${PROJECT_SOURCE_DIR}/Herwig)
add_library(ff OBJECT ffcard.f)
set_target_properties(ff PROPERTIES LINKER_LANGUAGE CXX POSITION_INDEPENDENT_CODE ON)
