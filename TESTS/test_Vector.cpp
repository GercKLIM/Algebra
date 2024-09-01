/*
 *  ТЕСТИРОВАНИЕ АЛГЕБРЫ ВЕКТОРОВ
 *
 * */

#include <gtest/gtest.h>
#include "../Lib/algebra.h"



/* Тесты: Создание объекта вектора */


TEST(Vectors, Create_Vector_from_size) {

    Vector<double> vec(3);
    vec[0] = 1;
    vec[1] = 2;
    vec[2] = 3;

    std::vector<double> std_vec = {1, 2, 3};
    EXPECT_EQ(vec.to_std(), std_vec);
}

TEST(Vectors, Create_Vector_from_size_val) {

    Vector<double> vec(3, 1);
    std::vector<double> std_vec = {1, 1, 1};

    EXPECT_EQ(vec.to_std(), std_vec);
}

TEST(Vectors, Create_Vector_from_massive) {

    Vector<double> vec = {1, 2, 3};
    std::vector<double> std_vec = {1, 2, 3};

    EXPECT_EQ(vec.to_std(), std_vec);
}

TEST(Vectors, Create_Vector_from_std_vector) {

    std::vector<double> std_vec = {1, 2, 3};
    Vector<double> vec(std_vec);

    EXPECT_EQ(vec.to_std(), std_vec);
}

TEST(Vectors, size) {

    Vector<double> vec = {1, 2, 3};

    EXPECT_EQ(vec.size(), 3);
}

TEST(Vectors, operator_position) {

    Vector<double> vec = {1, 2, 3};

    EXPECT_EQ(((vec[0] == 1) and (vec[1] == 2) and (vec[2] == 3)), true);
}

TEST(Vectors, operator_plus) {

    Vector<int> vec1 = {1, 2, 3};
    Vector<int> vec2 = {1, 2, 3};
    Vector<int> vec3 = vec1 + vec2;
    std::vector<int> std_vec3 = {2, 4, 6};

    EXPECT_EQ(vec3.to_std(), std_vec3);
}

TEST(Vectors, operator_minus) {

    Vector<int> vec1 = {1, 2, 3};
    Vector<int> vec2 = {1, 2, 3};
    Vector<int> vec3 = vec1 - vec2;
    std::vector<int> std_vec3 = {0, 0, 0};

    EXPECT_EQ(vec3.to_std(), std_vec3);
}

TEST(Vectors, operator_multiplication_vectors) {

    Vector<int> vec1 = {1, 2, 3};
    Vector<int> vec2 = {1, 2, 3};
    Vector<int> vec3 = vec1 * vec2;
    std::vector<int> std_vec3 = {1, 4, 9};

    EXPECT_EQ(vec3.to_std(), std_vec3);
}


TEST(Vectors, operator_multiplication_scalar) {

    Vector<int> vec1 = {1, 2, 3};
    Vector<int> vec3 = 2 * vec1;
    std::vector<int> std_vec3 = {2, 4, 6};

    EXPECT_EQ(vec3.to_std(), std_vec3);
}


TEST(Vectors, operator_div_vector) {

    Vector<int> vec1 = {1, 2, 3};
    Vector<int> vec2 = {1, 2, 3};
    Vector<int> vec3 = vec1 / vec2;
    std::vector<int> std_vec3 = {1, 1, 1};

    EXPECT_EQ(vec3.to_std(), std_vec3);
}

TEST(Vectors, operator_div_scalar) {

    Vector<int> vec1 = {2, 4, 6};
    Vector<int> vec3 = vec1 / 2;
    std::vector<int> std_vec3 = {1, 2, 3};

    EXPECT_EQ(vec3.to_std(), std_vec3);
}

/* Тесты: функций вывода вектора -> без тестирования */


/* Тесты: Функций других операций  */

TEST(Vectors, dot) {
    Vector<double> vec1 = {1, 2, 3};
    Vector<double> vec2 = {1, 2, 3};

    EXPECT_EQ(vec1.dot(vec2), 14);

}

TEST(Vectors, norm) {
    Vector<double> vec = {1, 2, 2};

    EXPECT_EQ(vec.norm(2), 3);
}

TEST(Vectors, abs) {
    Vector<double> vec = {-1, -2, -3};
    std::vector<double> std_vec = {1, 2, 3};

    EXPECT_EQ((vec.abs()).to_std(), std_vec);
}

TEST(Vectors, sum) {
    Vector<double> vec = {1, 2, 3};

    EXPECT_EQ(vec.sum(), 6);
}

TEST(Vectors, sort) {
    Vector<double> vec = {3, 2, 1};
    std::vector<double> std_vec = {1, 2, 3};

    EXPECT_EQ((vec.sort()).to_std(), std_vec);

}

TEST(Vectors, max) {
    Vector<double> vec = {1, 3, 2};

    EXPECT_EQ(vec.max(), 3);
}

TEST(Vectors, shift) {
    Vector<double> vec = {1, 2, 3};
    std::vector<double> std_vec = {3, 1, 2};

    EXPECT_EQ(vec.shift(1).to_std(), std_vec);
}




