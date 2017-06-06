set(processes sqme_aaaa sqme_aazz sqme_aaww)

#----- build all dependencies

foreach(_proc ${processes})
  add_subdirectory(${_proc})
endforeach()

#----- build a big library out of all sub-processes

add_library(comphep SHARED comphep_wraper.cpp
                           $<TARGET_OBJECTS:sqme_aaaa>
                           $<TARGET_OBJECTS:sqme_aaww>
                           $<TARGET_OBJECTS:sqme_aazz>)
