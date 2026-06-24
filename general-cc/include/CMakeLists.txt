
target_sources(general_cc_lib PUBLIC FILE_SET headers TYPE HEADERS FILES general_cc/orbit-string.hpp general_cc/utils.hpp)

if(GeneralCC_ENABLE_TESTING)
	target_sources(general_cc_lib PUBLIC FILE_SET headers TYPE HEADERS FILES testing/strings.hpp)
endif()