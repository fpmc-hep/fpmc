file(GLOB softc_sources *.f)

add_library(softc OBJECT ${softc_sources})
set_target_properties(softc PROPERTIES LINKER_LANGUAGE CXX POSITION_INDEPENDENT_CODE ON)
