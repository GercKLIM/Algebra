cmake_minimum_required(VERSION 3.27)
project(My_CPP_Libs)

set(CMAKE_CXX_STANDARD 17)

add_executable(My_CPP_Libs
        main.cpp

        Lib
        Lib/algebra.h

        Lib/MathFiles
        Lib/MathFiles/MathFiles.h
        Lib/MathFiles/MathFiles.tpp

        Lib/vectors
        Lib/vectors/vectors.h
        Lib/vectors/vectors.tpp

        Lib/matrixes
        Lib/matrixes/matrixes.h
        Lib/matrixes/matrixes.tpp

        Lib/SLAE
        Lib/SLAE/SLAE.h
        Lib/SLAE/SLAE.tpp

        Lib/maths
        Lib/maths/maths.h
        Lib/maths/maths.tpp
)
