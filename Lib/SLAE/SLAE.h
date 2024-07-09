/*
 *  ОБЪЯВЛЕНИЕ АЛГЕБРЫ CЛАУ
 *
 * */

#pragma once
#include "../algebra.h"


/* ### ФУНКЦИИ КОНВЕРТАЦИЙ в СЛАУ ### */


/* Функция для получения матрицы из СЛАУ */
template <typename T>
std::vector<std::vector<T>> SLAU_to_matrix(const std::vector<std::vector<T>>& SLAU);

/* Функция для получения векторая из СЛАУ */
template <typename T>
std::vector<T> SLAU_to_vec(const std::vector<std::vector<T>>& SLAU);


/* ### ФУНКЦИИ РЕШЕНИЯ СЛАУ ### */


/* ### ВЫЗОВ РЕАЛИЗАЦИИ ### */
#include "SLAE.tpp"