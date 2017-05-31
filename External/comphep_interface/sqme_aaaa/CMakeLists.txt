include_directories(comphep/include)
file(GLOB sources src/*.cpp)

add_library(sqme_aaaa OBJECT ${sources})
set_target_properties(sqme_aaaa PROPERTIES LINKER_LANGUAGE CXX POSITION_INDEPENDENT_CODE ON)

