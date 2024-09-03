/*
 *  ТЕСТИРОВАНИЕ АЛГЕБРЫ СЛАУ
 *
 * */

#include <gtest/gtest.h>
#include "../Lib/algebra.h"



/* Тесты: Создание объекта СЛАУ*/


TEST(SLAE, Create_from_Matrix_Vector) {
    Matrix<double> A = {{1, 2}, {3, 4}};
    Vector<double> b = {5, 6};
    SLAE<double> slae(A, b);
    EXPECT_EQ(slae.size(), 2);
    EXPECT_EQ(slae.to_Matrix(), A);
    EXPECT_EQ(slae.to_Vector(), b);
}

TEST(SLAE, Create_from_stdvecvec_stdvec) {
    std::vector<std::vector<double>> A = {{1, 2}, {3, 4}};
    std::vector<double> b = {5, 6};
    SLAE<double> slae(A, b);
    EXPECT_EQ(slae.size(), 2);
    EXPECT_EQ(slae.to_Matrix(), Matrix<double>(A));
    EXPECT_EQ(slae.to_Vector(), Vector<double>(b));
}

TEST(SLAE, Create_from_Matrix_stdvec) {
    Matrix<double> A = {{1, 2}, {3, 4}};
    std::vector<double> b = {5, 6};
    SLAE<double> slae(A, b);
    EXPECT_EQ(slae.size(), 2);
    EXPECT_EQ(slae.to_Matrix(), A);
    EXPECT_EQ(slae.to_Vector(), Vector<double>(b));
}

TEST(SLAE, Create_from_stdvecvec_Vector) {
    std::vector<std::vector<double>> A = {{1, 2}, {3, 4}};
    Vector<double> b = {5, 6};
    SLAE<double> slae(A, b);
    EXPECT_EQ(slae.size(), 2);
    EXPECT_EQ(slae.to_Matrix(), Matrix<double>(A));
    EXPECT_EQ(slae.to_Vector(), b);
}



/* ТЕСТЫ: операций со СЛАУ */



TEST(SLAE, operator_comparison) {
    Matrix<double> A = {{1, 2}, {3, 4}};
    Vector<double> b = {5, 6};
    SLAE<double> slae1(A, b);
    SLAE<double> slae2(A, b);
    EXPECT_TRUE(slae1 == slae2);
}

TEST(SLAE, residual) {
    Matrix<double> A = {{1, 2}, {3, 4}};
    Vector<double> b = {5, 6};
    Vector<double> sol = {-4, 4.5}; // Пример решения, надо найти правильное решение
    SLAE<double> slae(A, b);
    Vector<double> res = slae.residual(sol);
    EXPECT_NEAR(res[0], 0.0, 1e-9);
    EXPECT_NEAR(res[1], 0.0, 1e-9);
}

TEST(SLAE, residual_norm) {
    Matrix<double> A = {{1, 2}, {3, 4}};
    Vector<double> b = {5, 6};
    Vector<double> sol = {-4, 4.5}; // Пример решения, надо найти правильное решение
    SLAE<double> slae(A, b);
    double res_norm = slae.residual(sol, 2);
    EXPECT_NEAR(res_norm, 0.0, 1e-9);
}

TEST(SLAE, is_sol) {
    Matrix<double> A = {{1, 2}, {3, 4}};
    Vector<double> b = {5, 6};
    Vector<double> sol = {-4, 4.5}; // Пример решения, надо найти правильное решение
    SLAE<double> slae(A, b);
    EXPECT_TRUE(slae.is_sol(sol));
}

//TEST(SLAE, QR_decomposition) {
//    Matrix<double> A = {{1, 2}, {3, 4}};
//    SLAE<double> slae(A, {0, 0});
//    auto [Q, R] = slae.QR();
//    Matrix<double> expectedQ = ... // Ожидаемое значение Q
//    Matrix<double> expectedR = ... // Ожидаемое значение R
//    EXPECT_EQ(Q, expectedQ);
//    EXPECT_EQ(R, expectedR);
//}

//TEST(SLAE, LU_decomposition) {
//    Matrix<double> A = {{1, 2}, {3, 4}};
//    SLAE<double> slae(A, {0, 0});
//    auto [L, U] = slae.LU();
//    Matrix<double> expectedL = ... // Ожидаемое значение L
//    Matrix<double> expectedU = ... // Ожидаемое значение U
//    EXPECT_EQ(L, expectedL);
//    EXPECT_EQ(U, expectedU);
//}

TEST(SLAE, forward_Gauss) {
    Matrix<double> A = {{2, -1, 0}, {-1, 2, -1}, {0, -1, 2}};
    Vector<double> b = {1, 0, 1};
    SLAE<double> slae(A, b);
    slae.forward_Gauss(1e-10);
    Matrix<double> expectedA = {{1, -0.5, 0},
                                {0, 1, -2.0/3.0},
                                {0, 0, 1}};
    EXPECT_EQ(slae.to_Matrix(), expectedA);
}



TEST(SLAE, back_Gauss) {
    Matrix<double> A = {{2, -1, 0}, {-1, 2, -1}, {0, -1, 2}};
    Vector<double> b = {1, 0, 1};
    SLAE<double> slae(A, b);
    slae.forward_Gauss(1e-9);
    Vector<double> sol = slae.back_Gauss();
    Vector<double> expected_sol = {1, 1, 1}; // Ожидаемое решение
    //EXPECT_EQ(sol, expected_sol);

    for (int i = 0; i < sol.size(); ++i) {
        EXPECT_NEAR(sol[i], expected_sol[i], 1e-6);  // Учитываем погрешности округления
    }

}




/* ТЕСТЫ: решения СЛАУ */


TEST(SLAE, solve) {
    Matrix<double> A = {{1, 2}, {3, 4}};
    Vector<double> b = {5, 6};
    SLAE<double> slae(A, b);
    Vector<double> sol = slae.solve(1e-9);
    Vector<double> expected_sol = {-4, 4.5}; // Ожидаемое решение
    EXPECT_EQ(sol, expected_sol);
}

TEST(SLAE, solve_Gauss_test1) {
    Matrix<double> A = {{2, -1, 0}, {-1, 2, -1}, {0, -1, 2}};
    Vector<double> b = {1, 0, 1};
    SLAE<double> slae(A, b);
    Vector<double> sol = slae.solve_Gauss_2(1e-9);
    Vector<double> expected_sol = {1., 1., 1.}; // Ожидаемое решение
    //EXPECT_EQ(sol, expected_sol);
    for (int i = 0; i < sol.size(); ++i) {
        EXPECT_NEAR(sol[i], expected_sol[i], 1e-6);  // Учитываем погрешности округления
    }
}

// Тесты для других решений
TEST(SLAE, solve_Cramer_test1) {
    Matrix<double> A = {{1, 2}, {3, 4}};
    Vector<double> b = {5, 6};
    SLAE<double> slae(A, b);
    Vector<double> sol = slae.solve_Cramer(1e-9);
    Vector<double> expected_sol = {-4, 4.5}; // Ожидаемое решение
    EXPECT_EQ(sol, expected_sol);
}
