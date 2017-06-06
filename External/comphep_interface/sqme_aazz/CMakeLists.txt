include_directories(comphep/include)
file(GLOB sources src/*.cpp)

add_library(sqme_aazz OBJECT ${sources})
set_target_properties(sqme_aazz PROPERTIES LINKER_LANGUAGE CXX POSITION_INDEPENDENT_CODE ON)

