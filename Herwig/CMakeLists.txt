add_library(FpmcHerwig SHARED herwig6500.f)
set_target_properties(FpmcHerwig PROPERTIES LINKER_LANGUAGE CXX POSITION_INDEPENDENT_CODE ON)

add_library(HerwigCore OBJECT herwig6500.f)
set_target_properties(HerwigCore PROPERTIES LINKER_LANGUAGE CXX POSITION_INDEPENDENT_CODE ON)

#----- installation procedure

install(TARGETS HerwigCore DESTINATION lib)
