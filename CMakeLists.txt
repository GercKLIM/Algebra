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
        Lib/vectors/test_vectors.cpp
        Lib/SLAE/test_SLAE.cpp
        Lib/matrixes/test_matrixes.cpp
        Lib/maths/test_maths.cpp
        Lib/MathFiles/test_MathFiles.cpp
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
        Lib/vectors/test_vectors.cpp
)

# Тесты Matrix
add_executable(test_matrix
        Lib/matrixes/test_matrixes.cpp
)

# Тесты SLAE
add_executable(test_SLAE
        Lib/SLAE/test_SLAE.cpp
)

# Тесты Maths
add_executable(test_maths
        Lib/maths/test_maths.cpp
)

# Тесты MathFiles
add_executable(test_MathFiles
        Lib/MathFiles/test_MathFiles.cpp
)

add_executable(test_ALL
        Lib/vectors/test_vectors.cpp
        Lib/matrixes/test_matrixes.cpp
        Lib/SLAE/test_SLAE.cpp
        Lib/maths/test_maths.cpp
        Lib/MathFiles/test_MathFiles.cpp
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

