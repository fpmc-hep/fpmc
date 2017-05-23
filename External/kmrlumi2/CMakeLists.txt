file(GLOB kmr_lumi_2_sources *.f)

add_library(KMR_lumi_2 OBJECT ${kmr_lumi_2_sources})
set_target_properties(KMR_lumi_2 PROPERTIES LINKER_LANGUAGE CXX POSITION_INDEPENDENT_CODE ON)
