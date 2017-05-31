set(processes sqme_aaaa sqme_aazz sqme_aaww)
#set(comphep_sources "")
#foreach(_proc ${processes})
#  file(GLOB sources "${_proc}/src/*.cpp")
#  include_directories(${_proc}/comphep/include)
#  add_library(${_proc} STATIC ${sources})
#  #set_target_properties(${_proc} PROPERTIES LINKER_LANGUAGE CXX POSITION_INDEPENDENT_CODE ON)                                                                                       
#  #message(STATUS ${_proc} " contains " ${sources})
#  #set(comphep_sources ${comphep_sources} ${sources})
#endforeach()
#message(STATUS ${comphep_sources})
#add_library(comphep SHARED comphep_wraper.cpp ${comphep_sources})
#set(GCC_COVERAGE_LINK_FLAGS "-DNUM_DOUBLE -DCOMPHEP")
#set_target_properties(comphep PROPERTIES LINKER_LANGUAGE CXX POSITION_INDEPENDENT_CODE ON)

#----- build all dependencies

foreach(_proc ${processes})
  add_subdirectory(${_proc})
endforeach()

#----- build a big library out of all sub-processes

add_library(comphep SHARED comphep_wraper.cpp
                           $<TARGET_OBJECTS:sqme_aaaa>
                           $<TARGET_OBJECTS:sqme_aaww>
                           $<TARGET_OBJECTS:sqme_aazz>)
