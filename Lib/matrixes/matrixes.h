/*
 *  ОБЪЯВЛЕНИЕ АЛГЕБРЫ МАТРИЦ
 *
 * */

#pragma once
#include "../algebra.h"


/* ### ПЕРЕГРУЗКИ ОПЕРАТОРОВ ДЛЯ std::vector<vector<T>> ### */


/* Матричное умножение */
template <typename T>
std::vector<std::vector<T>> operator*(const std::vector<std::vector<T>>& A, const std::vector<std::vector<T>>& B);


/* Функция поэлементного сложения матриц */
template <typename T>
std::vector<std::vector<T>> operator+(const std::vector<std::vector<T>>& A, const std::vector<std::vector<T>>& B);


/* Функция поэлементного вычитания матриц */
template <typename T>
std::vector<std::vector<T>> operator-(const std::vector<std::vector<T>>& A, const std::vector<std::vector<T>>& B);


/* Функция для умножения матрицы на вектор */
template <typename T>
std::vector<T> operator*(const std::vector<std::vector<T>>& matrix, const std::vector<T>& vec);


/* ### ФУНКЦИИ ВЫВОДА МАТРИЦ НА ЭКРАН ### */

/* Функция вывода матрицы на экран */
template <typename T>
void print(const std::vector<std::vector<T>>& matrix);



/* ### ФУНКЦИИ ДРУГИХ ОПЕРАЦИЙ С МАТРИЦАМИ ### */

/* Функция для транспонирования матрицы */
template <typename T>
std::vector<std::vector<T>> transpose(const std::vector<std::vector<T>>& A);


/* Функция для создания единичной матрицы размера n x n */
template <typename T>
std::vector<std::vector<T>> create_identity_matrix(const int& n);


/* Функция для поэлементного умножения матриц */
template <typename T>
std::vector<std::vector<T>> Multyply(const std::vector<std::vector<T>>& A, const std::vector<std::vector<T>>& B);


/* Функция округления чисел в матрицах */
template <typename T>
std::vector<std::vector<T>> Matrix_round(const std::vector<std::vector<T>>& A, const double& eps);


/* Функция для вычисления нормы матрицы */
template <typename T>
T norm(const std::vector<std::vector<T>>& matrix, const int& p = 2);


/* Функция для вычисления числа обусловленности матрицы c нормой 1*/
template <typename T>
T cond_1(const std::vector<std::vector<T>>& matrix);


/* Функция для вычисления числа обусловленности матрицы c нормой 2*/
template <typename T>
T cond_2(const std::vector<std::vector<T>>& matrix);


/* Функция для вычисления числа обусловленности матрицы c нормой oo*/
template <typename T>
T cond_oo(const std::vector<std::vector<T>>& matrix);


/* Функция поворота матрицы вправо */
template <typename T>
std::vector<std::vector<T>> RotateRight(const std::vector<std::vector<T>>& A);


/* Функция поворота матрицы влево */
template <typename T>
std::vector<std::vector<T>> RotateLeft(const std::vector<std::vector<T>>& A);


// Функция для обратной матрицы с проверкой на вырожденность c определенной точностью
template <typename T>
std::vector<std::vector<T>> inverseMatrix(const std::vector<std::vector<T>>& A, const T& eps);


// Функция для обратной матрицы с проверкой на вырожденность c определенной точностью
template <typename T>
std::vector<std::vector<T>> inverseMatrix(const std::vector<std::vector<T>>& A);


// Функция обрезки матрицы снизу и справа
template <typename T>
std::vector<std::vector<T>> crop_matrix(const std::vector<std::vector<T>>& A, const int& k);


/* Функция, вычисляющая определитель матрицы 4х4 */
template <typename T>
double det(const std::vector<std::vector<T>>& matrix);


/* ### ВЫЗОВ РЕАЛИЗАЦИИ ### */
#include "matrixes.tpp"