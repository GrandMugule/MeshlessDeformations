cmake_minimum_required(VERSION 3.1)

set(CMAKE_BUILD_TYPE Release)

# project configuration
project(MeshlessDeformations)

set(MSVC_RUNTIME "dynamic")

add_subdirectory(src)
add_subdirectory(animations)

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${CMAKE_CURRENT_SOURCE_DIR}/cmake)
# set(CMAKE_FILES_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/build)

# libigl
option(LIBIGL_WITH_OPENGL            "Use OpenGL"         ON)
option(LIBIGL_WITH_OPENGL_GLFW       "Use GLFW"           ON)

find_package(LIBIGL REQUIRED QUIET)
