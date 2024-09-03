/*
 *  ОБЪЯВЛЕНИЕ АЛГЕБРЫ УРАВНЕНИЙ (ПРОСТЫХ)
 *
 * */

#pragma once
#include "../Algebra.h"


template <typename T>
class Equation {
private:
    std::function<T(T)> equation;
    T sol;
public:



    /* ### СОЗДАНИЕ УРАВНЕНИЯ ### */



    /* Создание объекта через Lambda-функцию */
    Equation(std::function<T(T)> eq) : equation(eq) {}


    // Метод для установки нового уравнения
    void setEquation(std::function<T(T)> eq) { equation = eq; }



    /* ### ОПЕРАЦИИ С УРАВНЕНИЯМИ ### */


    /* Операция вычисления уравнения */
    T operator()(const T& x) const {
        return equation(x);
    }

    /* Операция вычисления уравнения */
    T operator()(T& x) const {
        return equation(x);
    }


    /* Операция сравнения уравнений */
    //template<typename Y>
    //friend bool operator==(const Equation<Y>& slae1, const Equation<Y>& slae2);


    /* Операция сложения уравнений */
    template<typename Y>
    friend Equation<Y> operator+(const Equation<Y>& eq1, const Equation<Y>& eq2);


    /* Операция сложения уравнения и числа */
    template<typename Y>
    friend Equation<Y> operator+(const Equation<Y>& eq, const Y& value);


    /* Операция вычитания уравнений */
    template<typename Y>
    friend Equation<Y> operator-(const Equation<Y>& eq1, const Equation<Y>& eq2);


    /* Операция вычитания уравнения и числа */
    template<typename Y>
    friend Equation<Y> operator-(const Equation<Y>& eq, const Y& value);


    /* Операция отрицания уравнения */
    template<typename Y>
    friend Equation<Y> operator-(const Equation<Y>& eq);


    /* Операция умножения уравнений */
    template<typename Y>
    friend Equation<Y> operator*(const Equation<Y>& eq1, const Equation<Y>& eq2);


    /* Операция умножения на число */
    template<typename Y>
    friend Equation<Y> operator*(const Equation<Y>& eq, const Y& value);


    /* Операция деления уравнения на число */
    template<typename Y>
    friend Equation<Y> operator/(const Equation<Y>& e1, const Y& value);



    /* ### ДРУГИЕ ОПЕРАЦИИ С УРАВНЕНИЯМИ ### */



    /* Функция проверки решения уравнения */
    bool is_sol(const T& sol);


    /* Функция численной производной уравнения в точке */
    T differential(const T& x, const T& eps);


    /* Функция численной производной уравнения в точке c максимальной точностью */
    T differential(const T& x);


    /* Функция n-й численной производной уравнения в точке */
    T differential(const T& x, const int& n, const T& eps);


    /* Функция n-й численной производной уравнения в точке c максимальной точностью */
    T differential(const T& x, const int& n);


    /* Функция численного интегрирования уравнения в точке */
    T integrate(const T& x, const T& a, const T& b, const T& eps);


    /* Функция численного интегрирования уравнения в точке c максимальной точностью*/
    T integrate(const T& x, const T& a, const T& b);



    /* ### РЕШЕНИЕ УРАВНЕНИЙ ### */



};

/* ### ВЫЗОВ РЕАЛИЗАЦИИ ### */
#include "Equation.tpp"
