cmake_minimum_required(VERSION 3.8)
add_compile_options(-std=c++11)

project( catalogue  )
include_directories(./)
file(GLOB SOURCE "src/main.cpp")

add_executable( catalogue ${SOURCE} )