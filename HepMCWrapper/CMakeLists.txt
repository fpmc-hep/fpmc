include_directories(${PROJECT_SOURCE_DIR}/Examples)
add_library(Export OBJECT hwaend_dummy.f Fpmc.cc fostream.cc f77out.f dummy.f)
set_target_properties(Export PROPERTIES LINKER_LANGUAGE CXX POSITION_INDEPENDENT_CODE ON)
