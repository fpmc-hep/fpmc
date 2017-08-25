file(GLOB kmr_lumi_sources *.f)

add_library(KMR_lumi OBJECT ${kmr_lumi_sources})
set_target_properties(KMR_lumi PROPERTIES LINKER_LANGUAGE CXX POSITION_INDEPENDENT_CODE ON)
