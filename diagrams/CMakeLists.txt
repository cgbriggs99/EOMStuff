cmake_minimum_required(VERSION 3.25)

project(diagrams CXX)

set(CMAKE_CXX_STANDARD 20)

include(GNUInstallDirs)
include(CMakePackageConfigHelpers)

add_library(diagrams SHARED)

add_subdirectory(include)

install(TARGETS diagrams EXPORT diagrams FILE_SET headers)

export(TARGETS diagrams NAMESPACE diagrams:: FILE diagrams.cmake)

install(EXPORT diagrams NAMESPACE diagrams:: DESTINATION ${CMAKE_INSTALL_DATAROOTDIR}/cmake/diagrams)

add_subdirectory(src)

configure_package_config_file(cmake/diagrams-config.cmake.in ${CMAKE_CURRENT_BINARY_DIR}/diagrams-config.cmake INSTALL_DESTINATION ${CMAKE_INSTALL_DATAROOTDIR}/cmake/diagrams)

write_basic_package_version_file(${CMAKE_CURRENT_BINARY_DIR}/diagrams-config-version.cmake VERSION 1.0.0 COMPATIBILITY SameMajorVersion)

install(FILES ${CMAKE_CURRENT_BINARY_DIR}/diagrams-config.cmake ${CMAKE_CURRENT_BINARY_DIR}/diagrams-config-version.cmake DESTINATION ${CMAKE_INSTALL_DATAROOTDIR}/cmake/diagrams)