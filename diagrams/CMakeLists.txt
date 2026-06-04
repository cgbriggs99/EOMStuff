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

install(TARGETS diagrams_exe EXPORT diagrams_exe)

install(EXPORT diagrams NAMESPACE diagrams DESTINATION share/diagrams)

install(EXPORT diagrams_exe NAMESPACE diagrams DESTINATION share/diagrams)