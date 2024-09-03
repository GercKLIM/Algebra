/*
 *  ОБЪЯВЛЕНИЕ АЛГЕБРЫ ВЕКТОРОВ
 *
 * */

#pragma once
#include "../Algebra.h"



template <typename T>
class Vector {
private:
    std::vector<T> data;


public:



    /* ### СОЗДАНИЕ ОБЪЕКТА КЛАССА ### */



    /* Конструктор по умолчанию */
    Vector() : data(1, 0){};


    /* Создание объекта по исходному std::vector */
    Vector(std::vector<T> vec) {
        data = vec;
    };


    /* Создание объекта по исходному вектору */
    Vector(const Vector<T>& vec) {
        int size = vec.size();
        data.resize(size);
        for (int i = 0; i < size; i++) {
            data[i] = vec[i];
        }
    }


    /* Создание объекта по исходной длине и значению */
    Vector(int size, T value = T()) : data(size, value) {};


    /* Создание объекта по массиву */
    Vector(std::initializer_list<T> initList) : data(initList) {}


    /* Возвращение значения длины вектора */
    int size() const {
        return data.size();
    }



    /* ### ПЕРЕГРУЗКА ОПЕРАТОРОВ ДЛЯ std::vector<T> ### */



    /* Доступ к элементам вектора  */
    T& operator[](int index) {
        return data[index];
    }


    /* Доступ к элементам const-вектора */
    const T& operator[](int index) const {
        return data[index];
    }


    /* Операция сравнения */
    template<typename Y>
    friend bool operator==(const Vector<Y>& v1, const Vector<Y>& v2);


    /* Операция сложения векторов */
    template<typename Y>
    friend Vector<Y> operator+(Vector<Y>& v1, Vector<Y>& v2);


    /* Операция cложения const векторов */
    template <typename Y>
    friend Vector<Y> operator+(const Vector<Y>& vec1, const  Vector<Y>& vec2);


    /* Операция вычитания векторов */
    template <typename Y>
    friend Vector<Y> operator-(Vector<Y>& vec1, Vector<Y>& vec2);


    /* Операция вычитания const векторов */
    template <typename Y>
    friend Vector<Y> operator-(const Vector<Y>& vec1, const Vector<Y>& vec2);


    /* Операция почленного умножения векторов */
    template <typename Y>
    friend Vector<Y> operator*(Vector<Y>& vec1, Vector<Y>& vec2);


    /* Операция почленного умножения const векторов */
    template <typename Y>
    friend Vector<Y> operator*(const Vector<Y>& vec1, const Vector<Y>& vec2);

    /* Операция умножения вектора на число */
    template <typename Y>
    friend Vector<Y> operator*(Y& c, Vector<Y>& vec2);


    /* Операция умножения const вектора на const число */
    template <typename Y>
    friend Vector<Y> operator*(const Y& c, const Vector<Y>& vec2);


    /* Операция умножения const вектора на  число */
    template <typename Y>
    friend Vector<Y> operator*(Y& c, const Vector<Y>& vec2);


    /* Операция умножения const вектора на  число */
    template <typename Y>
    friend Vector<Y> operator*(const Y& c, Vector<Y>& vec2);


    /* Операция деления вектора на число */
    template <typename Y>
    friend Vector<Y> operator/(Vector<Y>& vec, Y& c);


    /* Операция деления const вектора на const число */
    template <typename Y>
    friend Vector<Y> operator/(const Vector<Y>& vec, const Y& c);


    /* Операция деления вектора на const число */
    template <typename Y>
    friend Vector<Y> operator/(Vector<Y>& vec, const Y& c);


    /* Операция деления const вектора на число */
    template <typename Y>
    friend Vector<Y> operator/(const Vector<Y>& vec, Y& c);


    /* Операция почленного деления векторов */
    template <typename Y>
    friend Vector<Y> operator/(Vector<Y>& vec1, Vector<Y>& vec2);


    /* Операция почленного деления const векторов */
    template <typename Y>
    friend Vector<Y> operator/(const Vector<Y>& vec1, const Vector<Y>& vec2);


    /* Операция потокового вывода */
    template<typename Y>
    friend std::ostream& operator<<(std::ostream& os, Vector<T>& v);


    /* Операция потокового вывода для const */
    template<typename Y>
    friend std::ostream& operator<<(std::ostream& os, const Vector<Y>& v);



    /* ### ФУНКЦИИ ВЫВОДА ВЕКТОРА ### */



    /* Функция вывода вектора на экран */
    void print();


    /* Функция вывода обрезанного вектора на экран */
    void print(const int& n);


    /* Функция, которая красиво выводит вектор*/
    void print_vec();


    /* Функция, которая красиво обрезанный выводит вектор*/
    void print_vec(const int& n);



    /* ### ФУНКЦИИ ДРУГИХ ОПЕРАЦИЙ С ВЕКТОРАМИ ### */

    /* Функция преобразовывающая Vector -> std::vector*/
    std::vector<T> to_std();

    /* Функция для скалярного умножения векторов */
    T dot(const Vector<T>& vec2);


    /* Функция для скалярного умножения векторов */
    template<typename Y>
    friend Y dot(const Vector<Y>& vec1, const Vector<Y>& vec2);


    /* Функция для нормы вектора */
    T norm(const int& p = 2);



    /* Функция, возвращает вектор модулей */
    Vector<T> abs();


    /* Функция, возращающая сумму элементов вектора */
    T sum();


    /* Функция, сортирующая вектор */
    Vector<T> sort();


    /* Функция, возращающая максимальный элемент вектора */
    T max();


    /* Функция для сдвига вектора на n элементов */
    Vector<T> shift(int n);

};







/* ### ВЫЗОВ РЕАЛИЗАЦИИ ### */
#include "Vector.tpp"

