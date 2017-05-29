#include_directories(${PROJECT_SOURCE_DIR}/External ${PROJECT_SOURCE_DIR}/Herwig)

add_library(Export OBJECT hwaend_dummy.f)
set_target_properties(Export PROPERTIES LINKER_LANGUAGE CXX POSITION_INDEPENDENT_CODE ON)
