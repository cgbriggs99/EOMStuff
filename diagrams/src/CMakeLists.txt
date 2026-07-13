target_sources(diagrams PRIVATE diagrams/diagram.cpp diagrams/theory-input.cpp 
	diagrams/utils.cpp diagrams/interaction_generators.cpp diagrams/cc_builder.cpp)

add_executable(diagrams_exe main.cpp)

target_link_libraries(diagrams_exe PUBLIC diagrams)

set_target_properties(diagrams_exe PROPERTIES INSTALL_RPATH_LINK_PATH ON)
set_target_properties(
    diagrams_exe PROPERTIES BUILD_RPATH "${CMAKE_BUILD_RPATH}"
                         INSTALL_RPATH "${CMAKE_INSTALL_RPATH};${CMAKE_INSTALL_LIBDIR}"
  )

install(TARGETS diagrams_exe EXPORT diagrams_exe)

install(EXPORT diagrams_exe NAMESPACE diagrams DESTINATION ${CMAKE_INSTALL_DATAROOTDIR}/diagrams)