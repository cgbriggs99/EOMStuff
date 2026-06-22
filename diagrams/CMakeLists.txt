cmake_minimum_required(VERSION 3.25)

project(diagrams CXX)

set(CMAKE_CXX_STANDARD 20)

include(GNUInstallDirs)

add_library(diagrams SHARED diagram.cpp)

target_sources(diagrams PUBLIC FILE_SET headers TYPE HEADERS FILES diagram.hpp)

install(TARGETS diagrams EXPORT diagrams FILE_SET headers)

export(TARGETS diagrams NAMESPACE diagrams FILE diagrams.cmake)

add_executable(diagrams_exe main.cpp)

target_link_libraries(diagrams_exe PUBLIC diagrams)

set_target_properties(diagrams_exe PROPERTIES INSTALL_RPATH_LINK_PATH ON)
set_target_properties(
    diagrams_exe PROPERTIES BUILD_RPATH "${CMAKE_BUILD_RPATH}"
                         INSTALL_RPATH "${CMAKE_INSTALL_RPATH};${CMAKE_INSTALL_LIBDIR}"
  )

install(TARGETS diagrams_exe EXPORT diagrams_exe)

install(EXPORT diagrams NAMESPACE diagrams DESTINATION share/diagrams)

install(EXPORT diagrams_exe NAMESPACE diagrams DESTINATION share/diagrams)