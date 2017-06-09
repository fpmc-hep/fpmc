include_directories(${PROJECT_SOURCE_DIR}/Examples ${HEPMC_INCLUDE})
add_library(Wrapper OBJECT Wrapper.cc)
set_target_properties(Wrapper PROPERTIES LINKER_LANGUAGE CXX POSITION_INDEPENDENT_CODE ON)
