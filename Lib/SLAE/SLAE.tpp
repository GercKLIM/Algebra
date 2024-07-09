/*
 *  РЕАЛИЗАЦИЯ АЛГЕБРЫ CЛАУ
 *
 * */

#include "../algebra.h"


/* ### ФУНКЦИИ КОНВЕРТАЦИЙ в СЛАУ ### */


/* Функция для получения матрицы из СЛАУ */
template <typename T>
std::vector<std::vector<T>> SLAU_to_matrix(const std::vector<std::vector<T>>& SLAU){
    std::vector<std::vector<T>> matrix;
    matrix.resize(SLAU.size(), std::vector<T>(SLAU.size()));

    for (int i = 0; i < SLAU.size(); i++) {
        for (int j = 0; j < SLAU.size(); j++) {
            matrix[i][j] = SLAU[i][j];
        }
    }
    return matrix;
}


/* Функция для получения вектора из СЛАУ */
template <typename T>
std::vector<T> SLAU_to_vec(const std::vector<std::vector<T>>& SLAU){
    int s = SLAU.size();
    std::vector<T> vec(s);

    for (int i = 0; i < SLAU.size(); i++) {
        vec[i] = SLAU[i][s];
    }
    return vec;
}


/* ### ФУНКЦИИ РЕШЕНИЯ СЛАУ ### */
