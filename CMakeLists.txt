cmake_minimum_required(VERSION 3.27)
project(My_CPP_Libs)

set(CMAKE_CXX_STANDARD 17)

add_executable(My_CPP_Libs main.cpp
        Lib
        Lib/algebra.h
        Lib/algebra.tpp

)
