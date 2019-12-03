cmake_minimum_required(VERSION 3.1)

# project configuration
project(MeshlessDeformations)

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${CMAKE_CURRENT_SOURCE_DIR}/cmake)
set(EXECUTABLE_OUTPUT_PATH ${CMAKE_CURRENT_SOURCE_DIR}/bin)
set(CMAKE_FILES_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/build)

# libigl
option(LIBIGL_WITH_OPENGL            "Use OpenGL"         ON)
option(LIBIGL_WITH_OPENGL_GLFW       "Use GLFW"           ON)

find_package(LIBIGL REQUIRED QUIET)

# source files
file(GLOB_RECURSE SRCFILES ${CMAKE_CURRENT_SOURCE_DIR}/src/*)

# executables
add_executable(stretching ${SRCFILES})

target_link_libraries(stretching igl::core igl::opengl_glfw)
