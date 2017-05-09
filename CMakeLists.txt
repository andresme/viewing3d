cmake_minimum_required(VERSION 3.6)
project(viewing3d)
find_package(OpenGL REQUIRED)
find_package(GLUT REQUIRED)
include_directories( ${OPENGL_INCLUDE_DIRS}  ${GLUT_INCLUDE_DIRS} )

set(CMAKE_C_STANDARD 99)

set(SOURCE_FILES main.c util/math/mathUtils.c util/math/mathUtils.h struct/vertex.c struct/vertex.h struct/settings.c struct/settings.h)
add_executable(viewing3d ${SOURCE_FILES})


target_link_libraries(viewing3d m ${OPENGL_LIBRARIES} ${GLUT_LIBRARY} )