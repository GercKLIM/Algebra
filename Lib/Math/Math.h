/*
 *  ОБЪЯВЛЕНИЕ ДРУГИХ ФУНКЦИЙ АЛГЕБРЫ
 *
 * */

#pragma once
#include "../Algebra.h"


/* ### ФУНКЦИИ ВЫЧИСЛЕНИЯ ПРОЗВОДНЫХ ### */


/* Функция, численно вычисляющая произвоную в точке point по i переменной */
//double Differential(std::vector<double> (*F)(const std::vector<double>&), const std::vector<double>& point, const int& i, const double& eps);


/* Функция, вычисляющая градиент функции в точке point */
//std::vector<double> Gradient(std::vector<double> (*F)(const std::vector<double>&), const std::vector<double>& point, const double& eps);


/* ### ФУНКЦИИ ВЫЧИСЛЕНИЯ ИНТЕГРАЛОВ ### */


/* ### ФУНКЦИИ РЕШЕНИЯ НЕЛИНЕЙНЫХ УРАВНЕНИЙ ### */



/* ### ВЫЗОВ РЕАЛИЗАЦИИ ### */
#include "Math.tpp"
