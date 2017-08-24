#file(GLOB_RECURSE chide_sources *.f)
file(GLOB chide_sources "Common/*.f" "Processes/*.f" "Diphoton/*.f" "GG/*.f" "Higgs/*.f")

add_library(CHIDe OBJECT ${chide_sources})
set_target_properties(CHIDe PROPERTIES LINKER_LANGUAGE CXX POSITION_INDEPENDENT_CODE ON)
