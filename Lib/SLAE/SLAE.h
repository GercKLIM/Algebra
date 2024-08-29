/*
 *  ОБЪЯВЛЕНИЕ АЛГЕБРЫ CЛАУ
 *
 * */

#pragma once
#include "../algebra.h"


template <typename T>
class SLAE {
private:
    Matrix<T> A_data;
    Vector<T> b_data;
    Matrix<T> sol_data;

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



    /* ### ОПЕРАЦИИ СО СЛАУ ### */



    /* Функция для получения матрицы из СЛАУ */
    Matrix<T> to_matrix();


    /* Функция для получения векторая из СЛАУ */
    Vector<T> to_Vector();


    /* Функция для вычисления невязки найденного решения */
    T residual();


    /* Функция для вычисления невязки данного решения */
    T residual(const Vector<T> sol);


    /* Функция для проверки решения СЛАУ*/
    bool is_sol(const Vector<T> sol);


    /* Функция QR-разложение */
    std::pair<Matrix<T>, Matrix<T>> QR;


    /* Функция LU-разложение */
    std::pair<Matrix<T>, Matrix<T>> LU;




    /* ### РЕШЕНИЕ СЛАУ ### */

    /* Вывод вычисленного решения */
    Vector<T> sol();

    /* Решение СЛАУ c автоматическим выбором метода */
    Vector<T> sol(T eps);

    /* Решение СЛАУ методом Гаусса */
    Vector<T> sol_Gauss(T eps);

    /* Решение СЛАУ методом QR-разложения */
    Vector<T> sol_QR(T eps);

    /* Решение СЛАУ методом Крамера */
    Vector<T> sol_Cramer(T eps);


};


/* ### ВЫЗОВ РЕАЛИЗАЦИИ ### */
#include "SLAE.tpp"