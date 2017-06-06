include_directories(comphep/include)
file(GLOB sources src/*.cpp)

add_library(sqme_aaww OBJECT ${sources})
set_target_properties(sqme_aaww PROPERTIES LINKER_LANGUAGE CXX POSITION_INDEPENDENT_CODE ON)

