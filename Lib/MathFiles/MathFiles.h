/*
 *  РЕАЛИЗАЦИЯ ФУНКЦИЙ ИМПОРТА И ЭКСПОРТА
 *  ВЕКТОРОВ И МАТРИЦ ИЗ *.txt ФАЙЛОВ
 *
 * */

#pragma once
#include "../Algebra.h"


/* ### ФУНКЦИИ ИМПОРТА ### */


/* Функция импорта чисел из файла */
std::vector<double> ImportData(const std::string& filename);


/* Функция импорта матрицы из текстового файла*/
template <typename T>
std::vector<std::vector<T>> importSLAU(const std::string& filename);


/* Переопределение потока вывода для std::vector */
template<typename T>
std::ostream& operator<<(std::ostream& os, const std::vector<T>& vec);


/* ### ВЫЗОВ РЕАЛИЗАЦИИ ### */
#include "MathFiles.tpp"
