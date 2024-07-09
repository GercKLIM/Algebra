/*
 *  ОБЪЯВЛЕНИЕ АЛГЕБРЫ ВЕКТОРОВ
 *
 * */

#pragma once
#include "../algebra.h"



/* ### ПЕРЕГРУЗКА ОПЕРАТОРОВ ДЛЯ std::vector<T> ### */



/* Операция cложения векторов */
template <typename T>
std::vector<T> operator+(const std::vector<T>& vec1, const  std::vector<T>& vec2);


/* Операция вычитания векторов */
template <typename T>
std::vector<T> operator-(const std::vector<T>& vec1, const std::vector<T>& vec2);


/* Операция почленного умножения векторов */
template <typename T>
std::vector<T> operator*(const std::vector<T>& vec1, const std::vector<T>& vec2);


/* Операция умножения вектора на число */
template <typename T>
std::vector<T> operator*(const T& c, const std::vector<T>& vec2);


template <typename T>
std::vector<T> operator*(const std::vector<T>& vec2, const T& c);

/* Операция деления вектора на число */
template<typename T>
std::vector<T> operator/(const std::vector<T>& vec, const T& c);


/* Операция почленного деления векторов */
template <typename T>
std::vector<T> operator/(const std::vector<T>& vec1, const std::vector<T>& vec2);


/* Определение оператора отрицания для матрицы */
template <typename T>
std::vector<std::vector<T>> operator-(const std::vector<std::vector<T>>& matrix);



/* ### ФУНКЦИИ ВЫВОДА ВЕКТОРА ### */



/* Функция вывода вектора на экран */
template <typename T>
void print(const std::vector<T>& vec);


/* Функция вывода обрезанного вектора на экран */
template <typename T>
void print_short(const std::vector<T>& vec, const int& n);


/* Функция, которая красиво выводит вектор*/
template<typename T>
void print_vec(const std::vector<T>& vec);



/* ### ФУНКЦИИ ДРУГИХ ОПЕРАЦИЙ С ВЕКТОРАМИ ### */



/* Функция для скалярного умножения векторов */
template <typename T>
T dot(const std::vector<T>& vec1, const std::vector<T>& vec2);


/* Функция для нормы вектора */
template <typename T>
T norm(const std::vector<T>& vec, const int& p = 2);


/* Функция, которая возращает матрицу комбинаций элементов вектора */
template<typename T>
std::vector<std::vector<T>> generateCombinations(const std::vector<T>& vec);


/* Функция, возвращает вектор модулей */
template<typename T>
std::vector<T> vec_abs(const std::vector<T>& vec);


/* Функция, возращающая сумму элементов вектора */
template<typename T>
T sum(const std::vector<T>& vec);


/* Функция, сортирующая вектор */
template< typename T>
std::vector<T> sorted(const std::vector<T>& vec_not_sort);


/* Функция, возращающая максимальный элемент вектора */
template<typename T>
T vec_max(const std::vector<T>& vec);


/* Функция, вычисляющая норму разности векторов */
double sqr(std::vector<double> vec1, std::vector<double> vec2);


/* Функция для сдвига вектора на n элементов */
template<typename T>
std::vector<T> shift(const std::vector<T>& vec, int n);


/* ### ВЫЗОВ РЕАЛИЗАЦИИ ### */
#include "vectors.tpp"

