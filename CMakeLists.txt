cmake_minimum_required(VERSION 3.27)
project(My_CPP_Libs)

# Установка стандартов C++
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Основная цель
add_executable(My_CPP_Libs
        main.cpp
        # Включите другие исходные файлы по мере необходимости
        Lib/algebra.h
        Lib/MathFiles/MathFiles.h
        Lib/MathFiles/MathFiles.tpp
        Lib/vectors/vectors.h
        Lib/vectors/vectors.tpp
        Lib/matrixes/matrixes.h
        Lib/matrixes/matrixes.tpp
        Lib/SLAE/SLAE.h
        Lib/SLAE/SLAE.tpp
        Lib/maths/maths.h
        Lib/maths/maths.tpp

        # Файлы тестов -> не включать!
        #Lib/vectors/test_vectors.cpp
        #Lib/SLAE/test_SLAE.cpp
        #Lib/matrixes/test_matrixes.cpp
        #Lib/maths/test_maths.cpp
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
add_executable(test_vectors
        TESTS/test_vectors.cpp

        # Файлы библиотеки
        Lib/algebra.h
        Lib/MathFiles/MathFiles.h
        Lib/MathFiles/MathFiles.tpp
        Lib/vectors/vectors.h
        Lib/vectors/vectors.tpp
        Lib/matrixes/matrixes.h
        Lib/matrixes/matrixes.tpp
        Lib/SLAE/SLAE.h
        Lib/SLAE/SLAE.tpp
        Lib/maths/maths.h
        Lib/maths/maths.tpp
)

# Тесты Matrix
add_executable(test_matrix
        TESTS/test_matrixes.cpp

        # Файлы библиотеки
        Lib/algebra.h
        Lib/MathFiles/MathFiles.h
        Lib/MathFiles/MathFiles.tpp
        Lib/vectors/vectors.h
        Lib/vectors/vectors.tpp
        Lib/matrixes/matrixes.h
        Lib/matrixes/matrixes.tpp
        Lib/SLAE/SLAE.h
        Lib/SLAE/SLAE.tpp
        Lib/maths/maths.h
        Lib/maths/maths.tpp
)

# Тесты SLAE
add_executable(test_SLAE
        TESTS/test_SLAE.cpp

        # Файлы библиотеки
        Lib/algebra.h
        Lib/MathFiles/MathFiles.h
        Lib/MathFiles/MathFiles.tpp
        Lib/vectors/vectors.h
        Lib/vectors/vectors.tpp
        Lib/matrixes/matrixes.h
        Lib/matrixes/matrixes.tpp
        Lib/SLAE/SLAE.h
        Lib/SLAE/SLAE.tpp
        Lib/maths/maths.h
        Lib/maths/maths.tpp
)

# Тесты Maths
add_executable(test_maths
        TESTS/test_maths.cpp

        # Файлы библиотеки
        Lib/algebra.h
        Lib/MathFiles/MathFiles.h
        Lib/MathFiles/MathFiles.tpp
        Lib/vectors/vectors.h
        Lib/vectors/vectors.tpp
        Lib/matrixes/matrixes.h
        Lib/matrixes/matrixes.tpp
        Lib/SLAE/SLAE.h
        Lib/SLAE/SLAE.tpp
        Lib/maths/maths.h
        Lib/maths/maths.tpp
)

# Тесты MathFiles
add_executable(test_MathFiles
        TESTS/test_MathFiles.cpp

        # Файлы библиотеки
        Lib/algebra.h
        Lib/MathFiles/MathFiles.h
        Lib/MathFiles/MathFiles.tpp
        Lib/vectors/vectors.h
        Lib/vectors/vectors.tpp
        Lib/matrixes/matrixes.h
        Lib/matrixes/matrixes.tpp
        Lib/SLAE/SLAE.h
        Lib/SLAE/SLAE.tpp
        Lib/maths/maths.h
        Lib/maths/maths.tpp
)


# Конфигурация запуска всех тестов
add_executable(test_ALL
        TESTS/test_vectors.cpp
        TESTS/test_matrixes.cpp
        TESTS/test_SLAE.cpp
        TESTS/test_maths.cpp
        TESTS/test_MathFiles.cpp

        # Файлы библиотеки
        Lib/algebra.h
        Lib/MathFiles/MathFiles.h
        Lib/MathFiles/MathFiles.tpp
        Lib/vectors/vectors.h
        Lib/vectors/vectors.tpp
        Lib/matrixes/matrixes.h
        Lib/matrixes/matrixes.tpp
        Lib/SLAE/SLAE.h
        Lib/SLAE/SLAE.tpp
        Lib/maths/maths.h
        Lib/maths/maths.tpp
)


# Добавление библиотеки Google Test в цели
target_link_libraries(test_vectors gtest gtest_main)
target_link_libraries(test_matrix gtest gtest_main)
target_link_libraries(test_SLAE gtest gtest_main)
target_link_libraries(test_maths gtest gtest_main)
target_link_libraries(test_MathFiles gtest gtest_main)
target_link_libraries(test_ALL gtest gtest_main)


# Обнаружение и регистрация тестов
include(GoogleTest)
gtest_discover_tests(test_vectors)
gtest_discover_tests(test_matrix)
gtest_discover_tests(test_SLAE)
gtest_discover_tests(test_maths)
gtest_discover_tests(test_MathFiles)
gtest_discover_tests(test_ALL)

