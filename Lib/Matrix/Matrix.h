/*
 *  ОБЪЯВЛЕНИЕ АЛГЕБРЫ МАТРИЦ
 *
 * */

#pragma once
#include "../Algebra.h"




template<typename T>
class Matrix {
private:
    std::vector<Vector<T>> data;



public:

    /* ### СОЗДАНИЕ ОБЪЕКТА КЛАССА ### */


    /* Создание начального объекта */
    Matrix() : data(1, Vector<T>(1, 0)) {};

    /* Конструктор с инициализатором списка */
    Matrix(std::initializer_list<std::initializer_list<T>> initList) {
        for (const auto& row : initList) {
            data.emplace_back(row);
        }
    }

    /* Создание объекта по исходной длине */
    Matrix(int size) : data(size, Vector<T>(1, 0)) {};


    /* Создание объекта по исходной размерности */
    Matrix(int rows, int colls) : data(rows, Vector<T>(colls, 0)) {};


    /* Создание объекта по исходной размерности и значению */
    Matrix(int rows, int colls, T value) : data(rows, Vector<T>(colls, value)) {};


    /* Создание объекта по исходному std::vector<Vector<T>>*/
    Matrix(std::vector<Vector<T>> matrix) {data = matrix;};


    /* Создание объекта по исходному const std::vector<std::vector<T>> */
    Matrix(const std::vector<std::vector<T>>& matrix) {
        int rows = matrix.size();
        int cols = rows > 0 ? matrix[0].size() : 0;
        data.resize(rows, Vector<T>(cols));
        for (int i = 0; i < rows; ++i) {
            data[i] = Vector<T>(matrix[i]);
        }
    }

    /* Создание объекта по исходному std::vector<std::vector<T>> */
    Matrix(std::vector<std::vector<T>>& matrix) {
        int rows = matrix.size();
        int cols = rows > 0 ? matrix[0].size() : 0;
        data.resize(rows, Vector<T>(cols));
        for (int i = 0; i < rows; ++i) {
            data[i] = Vector<T>(matrix[i]);
        }
    }


    /* Создание объекта по исходной длине и std::vector<T>*/
    Matrix(int size, const std::vector<T>& vec) {
        data.resize(size);

        for (int i = 0; i < size; ++i) {
            data[i] = Vector<T>(vec);
        }
    }


    /* Создание объекта по исходной длине и Vector<T> */
    Matrix(int size, const Vector<T>& vec) {
        data.resize(size, Vector<T>(vec.size(), T(0)));
        for (int i = 0; i < size; ++i) {
            data[i] = vec;
        }
    }


    /* Возвращение значения размерности */
    int size() const {
        return data.size();
    }


    /* Возвращение значения количества строк */
    int rows() const {
        return data.size();
    }


    /* Возвращение значения количества столбцов */
    int cols() const {
        return data[0].size();
    }


    /* Возвращение столбца */
    Vector<T> col(const int& i) const {

        Vector<T> coll(data.size(), 0);
        for (int j = 0; j < data.size(); j++){
            coll[j] = data[i];
        }
        return coll;
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


    /* Оператор индексации для константного объекта */
    const T& operator()(int i, int j) const {
        return data[i][j];
    }

    /* Оператор индексации для неконстантного объекта */
    T& operator()(int i, int j) {
        return data[i][j];
    }

    /* Операция сравнения */
    template<typename Y>
    friend bool operator==(const Matrix<Y> A1, const Matrix<Y> A2);

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


    /* Функция для умножения матрицы на вектор справа*/
    template<typename Y>
    friend Vector<Y> operator*(Vector<Y>& vec, Matrix<Y>& matrix);

    /* Операция умножения матрицы на вектор слева */
    template<typename Y>
    friend Vector<Y> operator*(Vector<Y>& vec, Matrix<Y>& matrix);

    /* Функция для умножения const матрицы на const вектор */
    template <typename Y>
    friend Matrix<Y> operator*(const Matrix<Y>& matrix, const Matrix<Y>& vec);


    /* Определение оператора отрицания для матрицы */
    template <typename Y>
    friend Matrix<Y> operator-(Matrix<Y>& A);


    /* Определение оператора отрицания для const матрицы */
    template <typename Y>
    friend Matrix<Y> operator-(const Matrix<Y>& matrix);


    /* Операция для умножения матрицы на число слева*/
    template <typename Y>
    friend Matrix<Y> operator*(const Matrix<Y>& A, const Y& scalar);


    /* Операция для умножения матрицы на число справа*/
    template <typename Y>
    friend Matrix<Y> operator*(const Y& scalar, const Matrix<Y>& A);


    /* Определение потока вывода для const матрицы */
    template <typename Y>
    friend std::ostream& operator<<(std::ostream& os, const Matrix<Y>& matrix);


    /* ### ФУНКЦИИ ВЫВОДА МАТРИЦ НА ЭКРАН ### */



    /* Функция вывода матрицы на экран */
    void print();

    /* Функция крисового вывода матрицы на экран */
    void print_matr();



    /* ### ФУНКЦИИ ДРУГИХ ОПЕРАЦИЙ С МАТРИЦАМИ ### */



    /* Функция преобразование матрицы в std::vector<Vector<T>> */
    std::vector<Vector<T>> to_std();


    /* Функция преобразование матрицы в std::vector<std::vector<T>> */
    std::vector<std::vector<T>> to_stds();



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
    Matrix<T> rotateRight();


    /* Функция поворота матрицы влево */
    Matrix<T> rotateLeft();


    /* Функция для обратной матрицы с проверкой на вырожденность c определенной точностью */
    Matrix<T> inverse(const T& eps);


    /* Функция для обратной матрицы с проверкой на вырожденность */
    Matrix<T> inverse();


    /* Функция обрезки матрицы снизу и справа */
    Matrix<T> crop(const int& k);


    T det() const;
};



/* ### ВЫЗОВ РЕАЛИЗАЦИИ ### */
#include "Matrix.tpp"