/*
 *  РЕАЛИЗАЦИЯ ФУНКЦИЙ ИМПОРТА И ЭКСПОРТА
 *  ВЕКТОРОВ И МАТРИЦ ИЗ *.txt ФАЙЛОВ
 *
 * */

#pragma once
#include "../Algebra.h"


/* ### ФУНКЦИИ ИМПОРТА ### */


template <typename T>
class MathFile {
private:
    std::fstream file;
    std::string filename;

public:

    /* Создание объекта класса через путь к файлу */
    MathFile(const std::string& fname) : filename(fname) {}


    /* Функция записи матрицы в файл */
    bool write(const Matrix<T>& A);


    /* Функция записи вектора в файл */
    bool write(const Vector<T>& vec);


    /* Функция получения матрицы из файла */
    Matrix<T> getMatrix();


    /* Функция получение вектора из файла */
    Vector<T> getVector();


};

/* ### ВЫЗОВ РЕАЛИЗАЦИИ ### */
#include "MathFiles.tpp"
