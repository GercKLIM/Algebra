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



    /* ### ОПЕРАЦИИ СЛАУ ### */



    /* ### ДРУГИЕ ОПЕРАЦИИ СЛАУ ### */



    /* Функция для получения матрицы из СЛАУ */
    Matrix<T> SLAU_to_matrix();


    /* Функция для получения векторая из СЛАУ */
    Vector<T> SLAU_to_vec();


    /* Проверка решения СЛАУ через невязку */
    template <typename Y>
    friend T is_sol(const Vector<T> sol);



    /* ### РЕШЕНИЕ СЛАУ ### */



    /* Решение СЛАУ c автоматическим выбором метода */
    Vector<T> sol(T eps);

    /* Решение СЛАУ методом Гаусса */
    Vector<T> sol_Gauss(T eps);

    /* Решение СЛАУ методом QR-разложения */
    Vector<T> sol_QR(T eps);


};


/* ### ВЫЗОВ РЕАЛИЗАЦИИ ### */
#include "SLAE.tpp"