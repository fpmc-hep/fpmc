#include_directories(${PROJECT_SOURCE_DIR}/External ${PROJECT_SOURCE_DIR}/Herwig)

add_library(HerwigCore OBJECT herwig6500.f)
set_target_properties(HerwigCore PROPERTIES LINKER_LANGUAGE CXX POSITION_INDEPENDENT_CODE ON)

