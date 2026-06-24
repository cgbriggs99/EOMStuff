
add_subdirectory(general_cc)

if(${GeneralCC_ENABLE_TESTING})
	add_subdirectory(tests)
endif()	