include_directories(${CMAKE_SOURCE_DIR}/src)

file(
	GLOB
	animations_files
	*.cpp
)

foreach(f ${animations_files})
	get_filename_component(animation_name ${f} NAME_WE)
	add_executable(${animation_name} ${f})
	target_link_libraries(${animation_name} utils igl::core igl::opengl_glfw)
	install(PROGRAMS ${CMAKE_CURRENT_BINARY_DIR}/${animation_name}
		DESTINATION bin
		RENAME ${CMAKE_PROJECT_NAME}-${animation_name})
endforeach(f)
