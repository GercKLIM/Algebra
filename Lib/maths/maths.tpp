/*
 *  РЕАЛИЗАЦИЯ ДРУГИХ ФУНКЦИЙ АЛГЕБРЫ
 *
 * */

#include "../algebra.h"


/* ### ФУНКЦИИ ВЫЧИСЛЕНИЯ ПРОЗВОДНЫХ ### */


/* Функция, численно вычисляющая произвоную в точке point по i переменной */
double Differential(std::vector<double> (*F)(const std::vector<double>&), const std::vector<double>& point, const int& i, const double& eps){

    std::vector<double> left_point(point);
    left_point[i] -= eps;
    std::vector<double> right_point(point);
    right_point[i] += eps;

    return (F(right_point)[i] - F(left_point)[i]) / (2 * eps);
}


/* Функция, вычисляющая градиент функции в точке point */
std::vector<double> Gradient(std::vector<double> (*F)(const std::vector<double>&), const std::vector<double>& point, const double& eps){

    int N = point.size();
    std::vector<double> grad(N, 0);

    for (int i = 0; i < N; i++){
        grad[i] = Differential(F, point, i, eps);
    }
    return grad;
}


/* ### ФУНКЦИИ ВЫЧИСЛЕНИЯ ИНТЕГРАЛОВ ### */


/* ### ФУНКЦИИ РЕШЕНИЯ НЕЛИНЕЙНЫХ УРАВНЕНИЙ ### */
