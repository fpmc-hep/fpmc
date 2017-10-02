#----- list all test executables

file(GLOB executables_noroot RELATIVE ${PROJECT_SOURCE_DIR}/Tests *.cc)
file(GLOB executables_root RELATIVE ${PROJECT_SOURCE_DIR}/Tests *.cxx)

foreach(exec_src ${executables_noroot})
  string(REPLACE ".cc" "" exec_bin ${exec_src})
  add_executable(${exec_bin} ${exec_src})
  set_target_properties(${exec_bin} PROPERTIES EXCLUDE_FROM_ALL true)
  target_link_libraries(${exec_bin} ${LIBRARIES})
  #target_link_libraries(${exec_bin} FPMC)
  #add_test(${exec_bin} ${exec_bin})
  #set_tests_properties(${exec_bin} PROPERTIES PASS_REGULAR_EXPRESSION "Passed")
endforeach()

#----- specify the tests requiring ROOT

find_package(ROOT QUIET)
if(${ROOT_FOUND})
  message(STATUS "ROOT found in ${ROOT_INCLUDE_DIRS}")
  include_directories(${ROOT_INCLUDE_DIRS})
  link_directories(${ROOT_LIBRARY_DIR})

  foreach(exec_src ${executables_root})
    string(REPLACE ".cxx" "" exec_bin ${exec_src})
    add_executable(${exec_bin} ${exec_src})
    set_target_properties(${exec_bin} PROPERTIES EXCLUDE_FROM_ALL true)
    target_link_libraries(${exec_bin} FPMC ${ROOT_LIBRARIES})
  endforeach()
else()
  message(STATUS "ROOT not found! skipping these tests!")
endif()
