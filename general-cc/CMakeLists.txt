#
# @BEGIN LICENSE
#
# general_cc by Psi4 Developer, a plugin to:
#
# Psi4: an open-source quantum chemistry software package
#
# Copyright (c) 2007-2025 The Psi4 Developers.
#
# The copyrights for code used from other parties are included in
# the corresponding files.
#
# This file is part of Psi4.
#
# Psi4 is free software; you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as published by
# the Free Software Foundation, version 3.
#
# Psi4 is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public License along
# with Psi4; if not, write to the Free Software Foundation, Inc.,
# 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
#
# @END LICENSE
#

cmake_minimum_required(VERSION 3.19...4.0)
if (POLICY CMP0167)
    cmake_policy(SET CMP0167 OLD)
endif ()

set(CMAKE_CXX_STANDARD 20)

project(general_cc CXX)

if("${CMAKE_BUILD_TYPE}" STREQUAL "Release")
	set(testing_default OFF)
else()
	set(testing_default ON)
endif()

option(GeneralCC_ENABLE_TESTING "Whether to build unit tests." ${testing_default})

if(${GeneralCC_ENABLE_TESTING})
	enable_testing()
endif()

find_package(psi4 1 REQUIRED)

find_package(Einsums REQUIRED)

find_package(diagrams REQUIRED)

include(FetchContent)

fetchcontent_declare(
  Catch2
  URL https://github.com/catchorg/Catch2/archive/v3.15.1.tar.gz
  FIND_PACKAGE_ARGS
  3
)
fetchcontent_makeavailable(Catch2)


include(GNUInstallDirs)
include(CMakePackageConfigHelpers)

add_library(general_cc_lib SHARED)
add_psi4_plugin(general_cc plugin.cc)

target_link_libraries(general_cc_lib PUBLIC diagrams::diagrams Einsums::Einsums)

if(GeneralCC_ENABLE_TESTING)
	target_link_libraries(general_cc_lib PUBLIC Catch2::Catch2)
endif()

target_link_libraries(general_cc PRIVATE general_cc_lib)

add_subdirectory(include)
add_subdirectory(src)

# Install the binary so that the Python can find it.
install(TARGETS general_cc EXPORT general_cc FILE_SET headers LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}/general_cc)

export(TARGETS general_cc NAMESPACE general_cc:: FILE general_cc.cmake)

install(EXPORT general_cc NAMESPACE general_cc:: DESTINATION ${CMAKE_INSTALL_DATAROOTDIR}/cmake/general_cc)

# Install the plugin for Psi4 to use.
install(DIRECTORY general_cc DESTINATION ${CMAKE_INSTALL_LIBDIR})

configure_package_config_file(cmake/general_cc-config.cmake.in ${CMAKE_CURRENT_BINARY_DIR}/general_cc-config.cmake INSTALL_DESTINATION ${CMAKE_INSTALL_DATAROOTDIR}/cmake/general_cc)

write_basic_package_version_file(${CMAKE_CURRENT_BINARY_DIR}/general_cc-config-version.cmake VERSION 1.0.0 COMPATIBILITY SameMajorVersion)

install(FILES ${CMAKE_CURRENT_BINARY_DIR}/general_cc-config.cmake ${CMAKE_CURRENT_BINARY_DIR}/general_cc-config-version.cmake DESTINATION ${CMAKE_INSTALL_DATAROOTDIR}/cmake/general_cc)

