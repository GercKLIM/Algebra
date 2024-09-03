/*
 *  ТЕСТИРОВАНИЕ АЛГЕБРЫ УРАВНЕНИЙ (ПРОСТЫХ)
 * */

#include <gtest/gtest.h>
#include "../Lib/algebra.h"



/* ТЕСТЫ: ОПЕРАЦИИ С УРАВНЕНИЯМИ */

/* Тест для оператора сложения уравнений */
TEST(EquationTest, AdditionOfEquations) {
    Equation<double> eq1([](double x) { return x * x; });
    Equation<double> eq2([](double x) { return x + 2; });

    auto sum_eq = eq1 + eq2;

    EXPECT_DOUBLE_EQ(sum_eq(2.0), 8.0);  // 2*2 + 2 + 2 = 8
    EXPECT_DOUBLE_EQ(sum_eq(3.0), 14.0); // 3*3 + 3 + 2 = 14
}

/* Тест для оператора сложения уравнения и числа */
TEST(EquationTest, AdditionOfEquationAndNumber) {
    Equation<double> eq1([](double x) { return x * x; });

    auto sum_eq = eq1 + 3.0;

    EXPECT_DOUBLE_EQ(sum_eq(2.0), 7.0);  // 2*2 + 3 = 7
    EXPECT_DOUBLE_EQ(sum_eq(3.0), 12.0); // 3*3 + 3 = 12
}

/* Тест для оператора вычитания уравнений */
TEST(EquationTest, SubtractionOfEquations) {
    Equation<double> eq1([](double x) { return x * x; });
    Equation<double> eq2([](double x) { return x + 2; });

    auto sub_eq = eq1 - eq2;

    EXPECT_DOUBLE_EQ(sub_eq(2.0), 0.0);   // 2*2 - (2 + 2) = 0
    EXPECT_DOUBLE_EQ(sub_eq(3.0), 4.0);   // 3*3 - (3 + 2) = 4
}

/* Тест для оператора вычитания уравнения и числа */
TEST(EquationTest, SubtractionOfEquationAndNumber) {
    Equation<double> eq1([](double x) { return x * x; });

    auto sub_eq = eq1 - 4.0;

    EXPECT_DOUBLE_EQ(sub_eq(2.0), 0.0);   // 2*2 - 4 = 0
    EXPECT_DOUBLE_EQ(sub_eq(3.0), 5.0);   // 3*3 - 4 = 5
}

/* Тест для оператора отрицания уравнения */
TEST(EquationTest, NegationOfEquation) {
    Equation<double> eq1([](double x) { return x * x; });

    auto neg_eq = -eq1;

    EXPECT_DOUBLE_EQ(neg_eq(2.0), -4.0);  // -(2*2) = -4
    EXPECT_DOUBLE_EQ(neg_eq(3.0), -9.0);  // -(3*3) = -9
}

/* Тест для оператора умножения уравнений */
TEST(EquationTest, MultiplicationOfEquations) {
    Equation<double> eq1([](double x) { return x * 2; });
    Equation<double> eq2([](double x) { return x + 3; });

    auto mul_eq = eq1 * eq2;

    EXPECT_DOUBLE_EQ(mul_eq(2.0), 20.0);  // 2*2 * (2 + 3) = 20
    EXPECT_DOUBLE_EQ(mul_eq(3.0), 36.0);  // 3*2 * (3 + 3) = 36
}

/* Тест для оператора умножения уравнения и числа */
TEST(EquationTest, MultiplicationOfEquationAndNumber) {
    Equation<double> eq1([](double x) { return x + 2; });

    auto mul_eq = eq1 * 3.0;

    EXPECT_DOUBLE_EQ(mul_eq(2.0), 12.0);  // (2 + 2) * 3 = 12
    EXPECT_DOUBLE_EQ(mul_eq(3.0), 15.0);  // (3 + 2) * 3 = 15
}

/* Тест для оператора деления уравнения на число */
TEST(EquationTest, DivisionOfEquationAndNumber) {
    Equation<double> eq1([](double x) { return x + 4; });

    auto div_eq = eq1 / 2.0;

    EXPECT_DOUBLE_EQ(div_eq(2.0), 3.0);   // (2 + 4) / 2 = 3
    EXPECT_DOUBLE_EQ(div_eq(4.0), 4.0);   // (4 + 4) / 2 = 4
}



/* ТЕСТЫ: ДРУГИЕ ОПЕРАЦИИ С УРАВНЕНИЯМИ ### */



/* Тестирование метода evaluate */
TEST(EquationTest, operator_evaluate) {
    Equation<double> eq([](double x) { return x * x - 4; });
    EXPECT_DOUBLE_EQ(eq(2.0), 0.0);
    EXPECT_DOUBLE_EQ(eq(3.0), 5.0);
    EXPECT_DOUBLE_EQ(eq(-2.0), 0.0);
}

/* Тестирование метода is_sol */
TEST(EquationTest, is_sol) {
    Equation<double> eq([](double x) { return x * x - 4; });
    EXPECT_TRUE(eq.is_sol(2.0));
    EXPECT_TRUE(eq.is_sol(-2.0));
    EXPECT_FALSE(eq.is_sol(0.0));
}

/* Тестирование метода differential (первой производной) */
TEST(EquationTest, Differential) {
    Equation<double> eq([](double x) { return x * x; });
    double diff_at_2 = eq.differential(2.0, 1e-6);
    EXPECT_NEAR(diff_at_2, 4.0, 1e-4);
}

/* Тестирование метода differential для n-й производной */
TEST(EquationTest, NDifferential) {
    Equation<double> eq([](double x) { return std::pow(x, 3); });
    double second_diff_at_2 = eq.differential(2.0, 2, 1e-6);
    EXPECT_NEAR(second_diff_at_2, 12.0, 1e-4);
}

/* Тестирование метода integrate */
TEST(EquationTest, Integrate) {
    Equation<double> eq([](double x) { return x; });
    double integral = eq.integrate(0.0, 0.0, 2.0, 1e-6);
    EXPECT_NEAR(integral, 2.0, 1e-4);
}