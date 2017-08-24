file(GLOB_RECURSE excl_aaaa_sources *.cpp)

include_directories(commons)

add_library(excl_aaaa OBJECT excl_aaaa_wraper.cpp ${excl_aaaa_sources})
set_target_properties(excl_aaaa PROPERTIES LINKER_LANGUAGE CXX POSITION_INDEPENDENT_CODE ON)
