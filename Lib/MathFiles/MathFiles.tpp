/*
 *  РЕАЛИЗАЦИЯ ФУНКЦИЙ ИМПОРТА И ЭКСПОРТА
 *  ВЕКТОРОВ И МАТРИЦ ИЗ *.txt ФАЙЛОВ
 *
 * */


#include "../Algebra.h"


/* ### ФУНКЦИИ ИМПОРТА ### */



/* Функция записи матрицы в файл */
template <typename T>
bool MathFile<T>::write(const Matrix<T>& A) {
    file.open(filename, std::ios::out | std::ios::trunc);  // Открываем файл на запись
    if (!file.is_open()) {
        return false;  // Если файл не открылся, возвращаем false
    }

    // Записываем размер матрицы
    file << A.rows() << " " << A.cols() << std::endl;

    // Записываем элементы матрицы
    for (size_t i = 0; i < A.rows(); ++i) {
        for (size_t j = 0; j < A.cols(); ++j) {
            file << A(i, j) << " ";
        }
        file << std::endl;
    }

    file.close();  // Закрываем файл
    return true;
}

/* Функция записи вектора в файл */
template <typename T>
bool MathFile<T>::write(const Vector<T>& vec) {
    file.open(filename, std::ios::out | std::ios::trunc);  // Открываем файл на запись
    if (!file.is_open()) {
        return false;  // Если файл не открылся, возвращаем false
    }

    // Записываем размер вектора
    file << vec.size() << std::endl;

    // Записываем элементы вектора
    for (size_t i = 0; i < vec.size(); ++i) {
        file << vec[i] << " ";
    }

    file << std::endl;
    file.close();  // Закрываем файл
    return true;
}

/* Функция получения матрицы из файла */
template <typename T>
Matrix<T> MathFile<T>::getMatrix() {
    file.open(filename, std::ios::in);  // Открываем файл на чтение
    if (!file.is_open()) {
        throw std::runtime_error("Unable to open file");  // Выбрасываем исключение при неудаче
    }

    size_t rows, cols;
    file >> rows >> cols;  // Читаем размер матрицы

    Matrix<T> A(rows, cols);

    // Читаем элементы матрицы
    for (size_t i = 0; i < rows; ++i) {
        for (size_t j = 0; j < cols; ++j) {
            file >> A(i, j);
        }
    }

    file.close();  // Закрываем файл
    return A;
}

/* Функция получения вектора из файла */
template <typename T>
Vector<T> MathFile<T>::getVector() {
    file.open(filename, std::ios::in);  // Открываем файл на чтение
    if (!file.is_open()) {
        throw std::runtime_error("Unable to open file");  // Выбрасываем исключение при неудаче
    }

    size_t size;
    file >> size;  // Читаем размер вектора

    Vector<T> vec(size);

    // Читаем элементы вектора
    for (size_t i = 0; i < size; ++i) {
        file >> vec[i];
    }

    file.close();  // Закрываем файл
    return vec;
}
