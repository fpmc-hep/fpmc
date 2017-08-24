#----- list of dependencies to build

set(EXTERNAL_MODULES CHIDe kmrlumi kmrlumi2 excl_aaaa softc comphep_interface/sqme_aaaa comphep_interface/sqme_aaww comphep_interface/sqme_aazz)

#----- build all dependencies

foreach(_mod ${EXTERNAL_MODULES})
  add_subdirectory(${_mod})
endforeach()

#----- other dependencies

file(GLOB misc_sources sudakov.f tools.f pdf/*.f lininter.f)
add_library(miscellaneous OBJECT ${misc_sources})
set_target_properties(miscellaneous PROPERTIES LINKER_LANGUAGE CXX POSITION_INDEPENDENT_CODE ON)

#----- wrap all external objects

add_library(FpmcExternals SHARED
  $<TARGET_OBJECTS:CHIDe>
  $<TARGET_OBJECTS:KMR_lumi> $<TARGET_OBJECTS:KMR_lumi_2>
  $<TARGET_OBJECTS:softc>
  $<TARGET_OBJECTS:miscellaneous>
  comphep_interface/comphep_wraper.cpp $<TARGET_OBJECTS:sqme_aaaa> $<TARGET_OBJECTS:sqme_aaww> $<TARGET_OBJECTS:sqme_aazz>
  $<TARGET_OBJECTS:excl_aaaa>
)
set_target_properties(FpmcExternals PROPERTIES LINKER_LANGUAGE CXX)

#----- installation procedure

install(TARGETS FpmcExternals DESTINATION lib)

