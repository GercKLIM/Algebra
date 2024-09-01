cmake_minimum_required(VERSION 3.27)
project(My_CPP_Libs)

# Установка стандартов C++
set(CMAKE_CXX_STANDARD 17)
#set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Основная цель
add_executable(My_CPP_Libs
        main.cpp
        # Включите другие исходные файлы по мере необходимости
        Lib/Algebra.h
        Lib/MathFiles/MathFiles.h
        Lib/MathFiles/MathFiles.tpp
        Lib/Vectors/Vector.h
        Lib/Vectors/Vector.tpp
        Lib/Matrix/Matrix.h
        Lib/Matrix/Matrix.tpp
        Lib/SLAE/SLAE.h
        Lib/SLAE/SLAE.tpp
        Lib/Math/Math.h
        Lib/Math/Math.tpp
        Lib/Equation/Equation.h
        Lib/Equation/equation.tpp
        Lib/Equation/Equation.tpp

        # Файлы тестов -> не включать!
        #Lib/Vectors/test_vectors.cpp
        #Lib/SLAE/test_SLAE.cpp
        #Lib/Matrix/test_matrixes.cpp
        #Lib/Math/test_maths.cpp
        #Lib/MathFiles/test_MathFiles.cpp
)


# ### ТЕСТЫ ###

# Включение Google Test
include(FetchContent)
FetchContent_Declare(
        googletest
        URL https://github.com/google/googletest/archive/release-1.11.0.zip
)

FetchContent_MakeAvailable(googletest)
enable_testing()



### Добавление тестов ###

# Тесты Vectors
add_executable(test_Vector
        TESTS/test_Vector.cpp

        # Файлы библиотеки
        Lib/Algebra.h
        Lib/MathFiles/MathFiles.h
        Lib/MathFiles/MathFiles.tpp
        Lib/Vectors/Vector.h
        Lib/Vectors/Vector.tpp
        Lib/Matrix/Matrix.h
        Lib/Matrix/Matrix.tpp
        Lib/SLAE/SLAE.h
        Lib/SLAE/SLAE.tpp
        Lib/Math/Math.h
        Lib/Math/Math.tpp
        Lib/Equation/Equation.h
        Lib/Equation/Equation.tpp
)

# Тесты Matrix
add_executable(test_Matrix
        TESTS/test_Matrix.cpp

        # Файлы библиотеки
        Lib/Algebra.h
        Lib/MathFiles/MathFiles.h
        Lib/MathFiles/MathFiles.tpp
        Lib/Vectors/Vector.h
        Lib/Vectors/Vector.tpp
        Lib/Matrix/Matrix.h
        Lib/Matrix/Matrix.tpp
        Lib/SLAE/SLAE.h
        Lib/SLAE/SLAE.tpp
        Lib/Math/Math.h
        Lib/Math/Math.tpp
        Lib/Equation/Equation.h
        Lib/Equation/Equation.tpp
)

# Тесты SLAE
add_executable(test_SLAE
        TESTS/test_SLAE.cpp

        # Файлы библиотеки
        Lib/Algebra.h
        Lib/MathFiles/MathFiles.h
        Lib/MathFiles/MathFiles.tpp
        Lib/Vectors/Vector.h
        Lib/Vectors/Vector.tpp
        Lib/Matrix/Matrix.h
        Lib/Matrix/Matrix.tpp
        Lib/SLAE/SLAE.h
        Lib/SLAE/SLAE.tpp
        Lib/Math/Math.h
        Lib/Math/Math.tpp
        Lib/Equation/Equation.h
        Lib/Equation/Equation.tpp
)

# Тесты Equation
add_executable(test_Equation
        TESTS/test_Equation.cpp

        # Файлы библиотеки
        Lib/Algebra.h
        Lib/MathFiles/MathFiles.h
        Lib/MathFiles/MathFiles.tpp
        Lib/Vectors/Vector.h
        Lib/Vectors/Vector.tpp
        Lib/Matrix/Matrix.h
        Lib/Matrix/Matrix.tpp
        Lib/SLAE/SLAE.h
        Lib/SLAE/SLAE.tpp
        Lib/Math/Math.h
        Lib/Math/Math.tpp
        Lib/Equation/Equation.h
        Lib/Equation/Equation.tpp
)

# Тесты Math
add_executable(test_Math
        TESTS/test_Math.cpp

        # Файлы библиотеки
        Lib/Algebra.h
        Lib/MathFiles/MathFiles.h
        Lib/MathFiles/MathFiles.tpp
        Lib/Vectors/Vector.h
        Lib/Vectors/Vector.tpp
        Lib/Matrix/Matrix.h
        Lib/Matrix/Matrix.tpp
        Lib/SLAE/SLAE.h
        Lib/SLAE/SLAE.tpp
        Lib/Math/Math.h
        Lib/Math/Math.tpp
        Lib/Equation/Equation.h
        Lib/Equation/Equation.tpp
)

# Тесты MathFiles
add_executable(test_MathFiles
        TESTS/test_MathFiles.cpp

        # Файлы библиотеки
        Lib/Algebra.h
        Lib/MathFiles/MathFiles.h
        Lib/MathFiles/MathFiles.tpp
        Lib/Vectors/Vector.h
        Lib/Vectors/Vector.tpp
        Lib/Matrix/Matrix.h
        Lib/Matrix/Matrix.tpp
        Lib/SLAE/SLAE.h
        Lib/SLAE/SLAE.tpp
        Lib/Math/Math.h
        Lib/Math/Math.tpp
        Lib/Equation/Equation.h
        Lib/Equation/Equation.tpp
)


# Конфигурация запуска всех тестов
add_executable(test_ALL
        TESTS/test_Vector.cpp
        TESTS/test_Matrix.cpp
        TESTS/test_SLAE.cpp
        TESTS/test_Math.cpp
        TESTS/test_MathFiles.cpp
        TESTS/test_Equation.cpp

        # Файлы библиотеки
        Lib/Algebra.h
        Lib/MathFiles/MathFiles.h
        Lib/MathFiles/MathFiles.tpp
        Lib/Vectors/Vector.h
        Lib/Vectors/Vector.tpp
        Lib/Matrix/Matrix.h
        Lib/Matrix/Matrix.tpp
        Lib/SLAE/SLAE.h
        Lib/SLAE/SLAE.tpp
        Lib/Math/Math.h
        Lib/Math/Math.tpp
        Lib/Equation/Equation.h
        Lib/Equation/Equation.tpp
)


# Добавление библиотеки Google Test в цели
target_link_libraries(test_Vector gtest gtest_main)
target_link_libraries(test_Matrix gtest gtest_main)
target_link_libraries(test_SLAE gtest gtest_main)
target_link_libraries(test_Math gtest gtest_main)
target_link_libraries(test_Equation gtest gtest_main)
target_link_libraries(test_MathFiles gtest gtest_main)
target_link_libraries(test_ALL gtest gtest_main)


# Обнаружение и регистрация тестов
include(GoogleTest)
gtest_discover_tests(test_Vector)
gtest_discover_tests(test_Matrix)
gtest_discover_tests(test_SLAE)
gtest_discover_tests(test_Math)
gtest_discover_tests(test_MathFiles)
gtest_discover_tests(test_ALL)

