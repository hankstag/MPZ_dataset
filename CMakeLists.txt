cmake_minimum_required(VERSION 3.16)
project(ffield_loader)

list(PREPEND CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/cmake)

# Libigl
option(LIBIGL_WITH_OPENGL            "Use OpenGL"         ON)
option(LIBIGL_WITH_OPENGL_GLFW       "Use GLFW"           ON)
include(libigl)

# set up HighFive
option(HIGHFIVE_USE_EIGEN ON)
find_package(HDF5 REQUIRED)
add_library(highfive INTERFACE)
target_include_directories(highfive SYSTEM INTERFACE ext/HighFive/include/ ${HDF5_INCLUDE_DIRS})
target_link_libraries(highfive INTERFACE ${HDF5_LIBRARIES})

add_library(field
  src/field/field.cpp 
  src/field/util.cpp 
  src/field/tangent_frame.cpp 
  src/field/Claussen.cpp
  src/field/util.cpp
  src/field/turning_number.cpp
  src/field/angles.cpp
  src/edge_map.cpp
)

target_link_libraries(field igl::core highfive)
target_include_directories(field PUBLIC src)

# Add your project files
add_executable(${PROJECT_NAME} src/main.cpp)
target_link_libraries(${PROJECT_NAME} PUBLIC igl::core field)
