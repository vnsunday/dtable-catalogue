cmake_minimum_required(VERSION 3.8)
add_compile_options(-std=c++17)

project( book_examples  )
include_directories(./)

add_executable( ex4_4 ex4_3.cpp )
add_executable( ex4_5 ex4_5.cpp )
add_executable( exII_RandomVariable_4_4 exII_RandomVariable_4_4.cpp )
add_executable( define_III_MultiVariables define_III_MultiVariables.cpp )