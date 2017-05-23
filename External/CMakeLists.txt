#----- list of dependencies to build

set(EXTERNAL_MODULES CHIDe
                     kmrlumi kmrlumi2
                     excl_aaaa comphep_interface
		     softc)

#----- build all dependencies

foreach(_mod ${EXTERNAL_MODULES})
  add_subdirectory(${_mod})
endforeach()

#----- link all libraries together

file(GLOB misc_sources sudakov.f tools.f pdf/*.f lininter.f)
add_library(miscellaneous OBJECT ${misc_sources})
set_target_properties(miscellaneous PROPERTIES LINKER_LANGUAGE CXX POSITION_INDEPENDENT_CODE ON)


add_library(Externals SHARED $<TARGET_OBJECTS:CHIDe>
                             $<TARGET_OBJECTS:KMR_lumi>
                             $<TARGET_OBJECTS:KMR_lumi_2>
                             $<TARGET_OBJECTS:excl_aaaa>
                             $<TARGET_OBJECTS:comphep>
                             $<TARGET_OBJECTS:softc>
                             $<TARGET_OBJECTS:miscellaneous>)
#set_target_properties(Externals PROPERTIES LINKER_LANGUAGE CXX POSITION_INDEPENDENT_CODE ON)
