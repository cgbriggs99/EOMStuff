
set(TESTS test_binomial.cpp test_strings.cpp)


foreach(test IN LISTS TESTS)
	cmake_path(GET test STEM name)
	
	add_executable(${name} ${test})
	
	target_link_libraries(${name} PUBLIC psi4::core Einsums::Einsums Catch2::Catch2WithMain general_cc_lib diagrams::diagrams)
	
	add_test(NAME ${name} COMMAND ${name})
endforeach()