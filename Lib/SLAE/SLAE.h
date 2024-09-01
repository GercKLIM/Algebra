/*
 *  ОБЪЯВЛЕНИЕ АЛГЕБРЫ CЛАУ
 *
 * */

#pragma once
#include "../Algebra.h"


template <typename T>
class SLAE {
private:
    Matrix<T> A_data;
    Vector<T> b_data;
    Vector<T> sol_data;

public:



    /* ### СОЗДАНИЕ СЛАУ ### */



    SLAE(Matrix<T> A, Vector<T> b) {
        A_data(A);
        b_data(b);
    }


    SLAE(std::vector<std::vector<T>> A, std::vector<T> b) {
        A_data(A);
        b_data(b);
    }


    SLAE(Matrix<T> A, std::vector<T> b) {
        A_data(A);
        b_data(b);
    }


    SLAE(std::vector<std::vector<T>> A, Vector<T> b) {
        A_data(A);
        b_data(b);
    }



    /* ### ФУНКЦИИ КОНВЕРТАЦИЙ в СЛАУ ### */



    /* Функция для получения матрицы из СЛАУ */
    Matrix<T> to_Matrix();


    /* Функция для получения векторая из СЛАУ */
    Vector<T> to_Vector();



    /* ### ОПЕРАЦИИ СО СЛАУ ### */


    /* Операция сравнения */
    template<typename Y>
    friend bool operator==(const SLAE<Y>& slae1, const SLAE<Y>& slae2);


    /* Функция размерности СЛАУ */
    int size() {
        return A_data.size();
    }


    /* Функция для вычисления невязки найденного решения */
//    T residual();


    /* Функция для вычисления невязки данного решения */
    Vector<T> residual(const Vector<T> sol);


    /* Функция для вычисления нормы невязки данного решения */
    T residual(const Vector<T> sol, const int& norm_type);


    /* Функция для проверки решения СЛАУ*/
    bool is_sol(const Vector<T> sol);


    /* Функция QR-разложение */
    std::pair<Matrix<T>, Matrix<T>> QR();


    /* Функция LU-разложение */
    std::pair<Matrix<T>, Matrix<T>> LU();


    /* Функция прямого хода метода Гаусса */
    void forward_Gauss(const T& eps);



    /* Функция обратного хода метода Гаусса */
    Vector<T> back_Gauss();



    /* ### РЕШЕНИЕ СЛАУ ### */



    /* Вывод вычисленного решения */
    Vector<T> sol();


    /* Решение СЛАУ c автоматическим выбором метода */
    Vector<T> solve(T eps);


    /* Решение СЛАУ методом Гаусса */
    Vector<T> solve_Gauss(T eps);


    /* Решение СЛАУ методом Гаусса */
    Vector<T> solve_Gauss_2(T eps);


    /* Решение СЛАУ методом QR-разложения */
    Vector<T> solve_QR(T eps);


    /* Решение СЛАУ методом Крамера */
    Vector<T> solve_Cramer(T eps);


};


/* ### ВЫЗОВ РЕАЛИЗАЦИИ ### */
#include "SLAE.tpp"