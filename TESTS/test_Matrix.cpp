/*
 *  ТЕСТИРОВАНИЕ АЛГЕБРЫ МАТРИЦ
 *
 * */

#include <gtest/gtest.h>
#include "../Lib/algebra.h"



/* Тесты: Создание объекта вектора */


TEST(Matrix, Create_Matrix) {

    Matrix<double> A;
    std::vector<std::vector<double>> std_A(1, std::vector<double>(1, 0));

    EXPECT_EQ(A.to_stds(), std_A);
}

TEST(Matrix, Create_Matrix_from_size) {

    Matrix<double> A(3);
    std::vector<std::vector<double>> std_A(3, std::vector<double>(1, 0));

    EXPECT_EQ(A.to_stds(), std_A);
}

TEST(Matrix, Create_Matrix_from_rows_colls) {

    Matrix<double> A(3, 3);
    std::vector<std::vector<double>> std_A(3, std::vector<double>(3, 0));

    EXPECT_EQ(A.to_stds(), std_A);
}

TEST(Matrix, Create_Matrix_from_rows_colls_value) {

    Matrix<double> A(3, 3, 1);
    std::vector<std::vector<double>> std_A(3, std::vector<double>(3, 1));

    EXPECT_EQ(A.to_stds(), std_A);
}

TEST(Matrix, Create_Matrix_from_stdvectorstdvector) {


    std::vector<std::vector<double>> std_A(3, std::vector<double>(3, 1));
    Matrix<double> A(std_A);

    EXPECT_EQ(A.to_stds(), std_A);
}

TEST(Matrix, Create_Matrix_from_stdvectorVector) {


    std::vector<Vector<double>> std_A(3, Vector<double>(3, 1));
    Matrix<double> A(std_A);

    EXPECT_EQ(A[0].to_std(), std_A[0].to_std());
}


TEST(Matrix, Create_Matrix_from_initlist) {


    std::vector<Vector<double>> std_A = {{1, 0, 0}, {0, 1, 0}, {0, 0, 1}};
    Matrix<double> A = {{1, 0, 0}, {0, 1, 0}, {0, 0, 1}};

    EXPECT_EQ(A.to_std(), std_A);
}



// TODO: Создание с аллокаторами
//TEST(Matrix, Create_Matrix_from_size_stdvector) {
//
//
//    std::vector<Vector<double>> std_A(3, Vector<double>(3, 1));
//    Matrix<double> A(3, std::vector<double>(3, 1));
//
//    EXPECT_EQ(A[0].to_std(), std_A[0].to_std());
//}

//TEST(Matrix, Create_Matrix_size_Vector) {
//
//
//    std::vector<Vector<double>> std_A(3, Vector<double>(3, 1));
//    Matrix<double> A(3, Vector<double>(3, 1));
//
//    EXPECT_EQ(A[0].to_std(), std_A[0].to_std());
//}


/* Тесты: операций с матрицами */


TEST(Matrix, Operator_position) {


    std::vector<Vector<double>> std_A(3, Vector<double>(3, 1));
    Matrix<double> A(std_A);

    EXPECT_EQ(A[0][0], std_A[0][0]);
}


TEST(Matrix, Operator_plus) {

    Matrix<double> A1(3, 3, 1);
    Matrix<double> A2(3, 3, 2);
    Matrix<double> A3(3, 3, 3);

    EXPECT_EQ((A1 + A2) == A3, true);
}

TEST(Matrix, Operator_minus) {

    Matrix<double> A1(3, 3, 3);
    Matrix<double> A2(3, 3, 1);
    Matrix<double> A3(3, 3, 2);

    EXPECT_EQ((A1 - A2) == A3, true);
}

TEST(Matrix, Operator_multiplication) {

    Matrix<double> A1(3, 3, 1);
    Matrix<double> A2(3, 3, 2);
    Matrix<double> A3(3, 3, 6);

    EXPECT_EQ((A1 * A2) == A3, true);
}

TEST(Matrix, Operator_multiplication_Vector_right) {

    Matrix<double> A1(3, 3, 2);
    Vector<double> vec(3, 3);
    Vector<double> vec2(3, 18);

    EXPECT_EQ((A1 * vec) == vec2, true);
}

TEST(Matrix, Operator_multiplication_Vector_left) {

    Matrix<double> A1(3, 3, 2);
    Vector<double> vec(3, 3);
    Vector<double> vec2(3, 18);

    EXPECT_EQ((vec * A1) == vec2, true);
}


TEST(Matrix, Operator_multiplication_scalar_left) {
    Matrix<double> A1(3, 3, 2);
    double scal = 3;
    Matrix<double> A2(3, 3, 6);

    EXPECT_EQ(((scal * A1)) == A2, true);
}

TEST(Matrix, Operator_multiplication_scalar_right) {
    Matrix<double> A1(3, 3, 2);
    double scal = 3;
    Matrix<double> A2(3, 3, 6);

    EXPECT_EQ(((A1 * scal) == A2), true);
}

TEST(Matrix, Operator_dis) {
    Matrix<double> A1(3, 3, 2);
    Matrix<double> A2(3, 3, -2);

    EXPECT_EQ(-A1, A2);
}

TEST(Matrix, Transpose) {
    Matrix<double> A1 = {{0, 1}, {0, 1}};
    Matrix<double> A2 = {{0, 0}, {1, 1}};

    EXPECT_EQ(A1.transpose(), A2);
}

TEST(Matrix, Create_identity) {
    Matrix<double> A1 = {{1, 0, 0}, {0, 1, 0},{0, 0, 1}};
    Matrix<double> A2;
    A2.create_identity_matrix(3);

    EXPECT_EQ(A1, A2);
}

TEST(Matrix, Multyply) {
    Matrix<double> A1 = {{2, 2}, {2, 2}};
    Matrix<double> A2 = {{2, 2}, {2, 2}};
    Matrix<double> A3 = {{8, 8}, {8, 8}};

    EXPECT_EQ(A1 * A2, A3);
}

TEST(Matrix, Round) {
    Matrix<double> A1 = {{2.01, 2.02}, {2.03, 2.04}};
    Matrix<double> A2 = {{2, 2}, {2, 2}};
    A1.round(0.1);

    EXPECT_EQ(A1, A2);
}

TEST(Matrix, Norm) {
    Matrix<double> A1 = {{2, 2}, {2, 2}};

    EXPECT_EQ(A1.norm(2) , 4);
}

TEST(Matrix, RotateLeft) {

    Matrix<double> A1 = {{1, 1}, {2, 2}};
    Matrix<double> A2 = {{2, 1}, {2, 1}};

    EXPECT_EQ(A1.rotateLeft(), A2);
}

TEST(Matrix, RotateRight) {

    Matrix<double> A2 = {{1, 1}, {2, 2}};
    Matrix<double> A1 = {{2, 1}, {2, 1}};

    EXPECT_EQ(A1.rotateRight(), A2);
}

TEST(Matrix, Inverse) {
    Matrix<double> A1 = {{1, 2}, {2, 2}};
    Matrix<double> A2 = {{-1, 1}, {1, -0.5}};

    EXPECT_EQ(A1.inverse(), A2);

}

TEST(Matrix, crop) {
    Matrix<double> A1 = {{1, 1, 1}, {2, 2, 2}, {3, 3, 3}};
    Matrix<double> A2 = {{1, 1}, {2, 2}};

    EXPECT_EQ(A1.crop(2) == A2, true);
}










