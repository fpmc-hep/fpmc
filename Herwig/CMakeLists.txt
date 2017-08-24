add_library(FpmcHerwig SHARED herwig6500.f)
set_target_properties(FpmcHerwig PROPERTIES LINKER_LANGUAGE CXX POSITION_INDEPENDENT_CODE ON)

#----- installation procedure

install(TARGETS FpmcHerwig DESTINATION lib)
