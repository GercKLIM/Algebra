/*
 *  ОБЪЯВЛЕНИЕ АЛГЕБРЫ МАТРИЦ
 *
 * */

#pragma once
#include "../algebra.h"




template<typename T>
class Matrix {
private:
    std::vector<Vector<T>> data;



public:

    /* ### СОЗДАНИЕ ОБЪЕКТА КЛАССА ### */


    /* Создание начального объекта */
    Matrix() : data(1, Vector<T>(1, 0)) {};


    /* Создание объекта по исходной размерности и значению */
    Matrix(int rows, int colls, T value) : data(rows, Vector<T>(colls, value)) {};


    /* Создание объекта по исходному std::vector<Vector<T>>*/
    Matrix(std::vector<Vector<T>> matrix) {data = matrix;};


    /* Создание объекта по исходному std::vector<std::vector<T>>*/
    Matrix(const std::vector<std::vector<T>>& matrix) {
        int rows = matrix.size();
        int cols = rows > 0 ? matrix[0].size() : 0;
        data.resize(rows, Vector<T>(cols));
        for (int i = 0; i < rows; ++i) {
            data[i] = Vector<T>(matrix[i]);
        }
    }


    Matrix(std::vector<std::vector<T>>& matrix) {
        int rows = matrix.size();
        int cols = rows > 0 ? matrix[0].size() : 0;
        data.resize(rows, Vector<T>(cols));
        for (int i = 0; i < rows; ++i) {
            data[i] = Vector<T>(matrix[i]);
        }
    }


    /* Создание объекта по исходной длине и  std::vector<T>*/
    Matrix(int size, std::vector<T> vec) {
        data(size, Vector<T>(vec.size(), 0));
        for (int i = 0; i < size; i++){
            Vector<T> Vec(vec);
            data[i] = Vec;
        }
    };


    /* Создание объекта по исходной длине и Vector<T>*/
    Matrix(int size, Vector<T> vec) {
        data(size, Vector<T>(vec.size(), 0));
        for (int i = 0; i < size; i++){
            data[i] = vec;
        }
    };


    /* Возвращение значения размерности */
    int size() const {
        return data.size();
    }



    /* ### ПЕРЕГРУЗКИ ОПЕРАТОРОВ ДЛЯ Matrix ### */



    /* Оператор индексации для неконстантного объекта */
    Vector<T>& operator[](int index) {
        return data[index];
    }

    /* Оператор индексации для константного объекта */
    const Vector<T>& operator[](int index) const {
        return data[index];
    }


    /* Операция поэлементного сложения матриц */
    template <typename Y>
    friend Matrix<Y> operator+(Matrix<Y>& A, Matrix<Y>& B);


    /* Операция поэлементного сложения const матриц */
    template <typename Y>
    friend Matrix<Y> operator+(const Matrix<Y>& A, const Matrix<Y>& B);


    /* Функция поэлементного вычитания матриц */
    template <typename Y>
    friend Matrix<Y> operator-(Matrix<Y>& A, Matrix<Y>& B);


    /* Функция поэлементного вычитания const матриц */
    template <typename Y>
    friend Matrix<Y> operator-(const Matrix<Y>& A, const Matrix<Y>& B);


    /* Матричное умножение */
    template <typename Y>
    friend Matrix<Y> operator*(Matrix<Y>& A, Matrix<Y>& B);


    /* Матричное умножение const */
    template <typename Y>
    friend Matrix<Y> operator*(const Matrix<Y>& A, const Matrix<Y>& B);


    /* Функция для умножения матрицы на вектор */
    template <typename Y>
    friend Matrix<Y> operator*(Matrix<Y>& matrix, Matrix<Y>& vec);


    /* Функция для умножения const матрицы на const вектор */
    template <typename Y>
    friend Matrix<Y> operator*(const Matrix<Y>& matrix, const Matrix<Y>& vec);


    /* Определение оператора отрицания для матрицы */
    template <typename Y>
    friend Matrix<Y> operator-(Matrix<Y>& matrix);


    /* Определение оператора отрицания для const матрицы */
    template <typename Y>
    friend Matrix<Y> operator-(const Matrix<Y>& matrix);


    /* Определение потока вывода для const матрицы */
    template <typename Y>
    friend std::ostream& operator<<(std::ostream& os, const Matrix<Y>& matrix);

    /* ### ФУНКЦИИ ВЫВОДА МАТРИЦ НА ЭКРАН ### */



    /* Функция вывода матрицы на экран */
    void print();

    /* Функция крисового вывода матрицы на экран */
    void print_matr();



    /* ### ФУНКЦИИ ДРУГИХ ОПЕРАЦИЙ С МАТРИЦАМИ ### */

    /* Функция для транспонирования матрицы */
    Matrix<T> transpose();


    /* Функция для создания единичной матрицы размера n x n */
    Matrix<T> create_identity_matrix(const int& n);


    /* Функция для поэлементного умножения матриц */
    Matrix<T> Multyply(const std::vector<std::vector<T>>& B);


    /* Функция округления чисел в матрицах */
    Matrix<T> round(const double& eps);


    /* Функция для вычисления нормы матрицы */
    T norm(const int& p);


    /* Функция поворота матрицы вправо */
    Matrix<T> RotateRight();


    /* Функция поворота матрицы влево */
    Matrix<T> RotateLeft();


    /* Функция для обратной матрицы с проверкой на вырожденность c определенной точностью */
    std::vector<std::vector<T>> inverseMatrix(const T& eps);


    /* Функция для обратной матрицы с проверкой на вырожденность */
    std::vector<std::vector<T>> inverseMatrix();


    /* Функция обрезки матрицы снизу и справа */
    std::vector<std::vector<T>> crop_matrix(const int& k);


};



/* ### ВЫЗОВ РЕАЛИЗАЦИИ ### */
#include "matrixes.tpp"