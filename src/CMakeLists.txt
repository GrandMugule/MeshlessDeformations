add_library(utils
	STATIC
	adjacency.cpp clustering.cpp integration.cpp shapematching.cpp)

target_link_libraries(utils igl::core)

install(TARGETS utils
	DESTINATION lib)

file(
	GLOB
	headers
	*.h
)

install(FILES ${headers}
	DESTINATION include/${CMAKE_PROJECT_NAME})
