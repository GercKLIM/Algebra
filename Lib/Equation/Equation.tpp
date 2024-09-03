/*
 *  РЕАЛИЗАЦИЯ АЛГЕБРЫ УРАВНЕНИЙ (ПРОСТЫХ)
 *
 * */

#include "../Algebra.h"


/* ### ОПЕРАЦИИ С УРАВНЕНИЯМИ ### */



// Операция сложения уравнений
template<typename Y>
Equation<Y> operator+(const Equation<Y>& eq1, const Equation<Y>& eq2) {
    return Equation<Y>([=](Y x) { return eq1(x) + eq2(x); });
}

// Операция сложения уравнения и числа
template<typename Y>
Equation<Y> operator+(const Equation<Y>& eq, const Y& value) {
    return Equation<Y>([=](Y x) { return eq(x) + value; });
}

// Операция вычитания уравнений
template<typename Y>
Equation<Y> operator-(const Equation<Y>& eq1, const Equation<Y>& eq2) {
    return Equation<Y>([=](Y x) { return eq1(x) - eq2(x); });
}

// Операция вычитания уравнения и числа
template<typename Y>
Equation<Y> operator-(const Equation<Y>& eq, const Y& value) {
    return Equation<Y>([=](Y x) { return eq(x) - value; });
}

// Операция отрицания уравнения
template<typename Y>
Equation<Y> operator-(const Equation<Y>& eq) {
    return Equation<Y>([=](Y x) { return -eq(x); });
}

// Операция умножения уравнений
template<typename Y>
Equation<Y> operator*(const Equation<Y>& eq1, const Equation<Y>& eq2) {
    return Equation<Y>([=](Y x) { return eq1(x) * eq2(x); });
}

// Операция умножения уравнения на число
template<typename Y>
Equation<Y> operator*(const Equation<Y>& eq, const Y& value) {
    return Equation<Y>([=](Y x) { return eq(x) * value; });
}

// Операция деления уравнения на число
template<typename Y>
Equation<Y> operator/(const Equation<Y>& eq, const Y& value) {
    return Equation<Y>([=](Y x) { return eq(x) / value; });
}



/* ### ДРУГИЕ ОПЕРАЦИИ С УРАВНЕНИЯМИ ### */


template<typename T>
bool Equation<T>::is_sol(const T& sol) {
    return std::abs(equation(sol)) < std::numeric_limits<T>::epsilon();
}


/* Функция численной производной уравнения в точке */
template<typename T>
T Equation<T>::differential(const T& x, const T& eps) {
    return (equation(x + eps) - equation(x)) / eps;
}


/* Функция численной производной уравнения в точке с максимальной точностью */
template<typename T>
T Equation<T>::differential(const T& x) {
    T eps = std::sqrt(std::numeric_limits<T>::epsilon());
    return differential(x, eps);
}


/* Функция n-й численной производной уравнения в точке */
template<typename T>
T Equation<T>::differential(const T& x, const int& n, const T& eps) {
    if (n == 0) return equation(x);
    T res = equation(x);
    for (int i = 0; i < n; ++i) {
        res = differential(x, eps);
    }
    return res;
}


/* Функция n-й численной производной уравнения в точке с максимальной точностью */
template<typename T>
T Equation<T>::differential(const T& x, const int& n) {
    T eps = std::sqrt(std::numeric_limits<T>::epsilon());
    return differential(x, n, eps);
}


/* Функция численного интегрирования уравнения */
template<typename T>
T Equation<T>::integrate(const T& x, const T& a, const T& b, const T& eps) {
    T n = (b - a) / eps;
    T sum = 0;
    for (T i = 0; i < n; ++i) {
        T xi = a + i * eps;
        sum += equation(xi) * eps;
    }
    return sum;
}


/* Функция численного интегрирования уравнения с максимальной точностью */
template<typename T>
T Equation<T>::integrate(const T& x, const T& a, const T& b) {
    T eps = std::sqrt(std::numeric_limits<T>::epsilon());
    return integrate(x, a, b, eps);
}



/* ### РЕШЕНИЕ УРАВНЕНИЙ ### */