/*
 *  ОБЪЯВЛЕНИЕ ДРУГИХ ФУНКЦИЙ АЛГЕБРЫ (в разработке)
 *
 * */

#pragma once
#include "../Algebra.h"




/* Переопределение потока вывода для vector  */
template<typename T>
std::ostream& operator<<(std::ostream& os, const std::vector<T>& vec) {
    //os << "[";
    for (int i = 0; i < vec.size(); ++i) {
        os << std::setprecision(16) << vec[i];
        if (i != vec.size() - 1) {
            os << " ";
            //os << ", ";
        }
    }
    //os << "]";
    os << " ";
    return os;
}



/* ### ФУНКЦИИ ВЫЧИСЛЕНИЯ ПРОЗВОДНЫХ ### */


/* Функция, численно вычисляющая произвоную в точке point по i переменной */
//double Differential(std::vector<double> (*F)(const std::vector<double>&), const std::vector<double>& point, const int& i, const double& eps);


/* Функция, вычисляющая градиент функции в точке point */
//std::vector<double> Gradient(std::vector<double> (*F)(const std::vector<double>&), const std::vector<double>& point, const double& eps);


/* ### ФУНКЦИИ ВЫЧИСЛЕНИЯ ИНТЕГРАЛОВ ### */


/* ### ФУНКЦИИ РЕШЕНИЯ НЕЛИНЕЙНЫХ УРАВНЕНИЙ ### */


/* ### ВЫЗОВ РЕАЛИЗАЦИИ ### */
#include "Math.tpp"
