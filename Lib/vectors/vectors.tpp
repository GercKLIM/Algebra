/*
 *  РЕАЛИЗАЦИЯ АЛГЕБРЫ ВЕКТОРОВ
 *
 * */

#include "../algebra.h"


/* ### ПЕРЕГРУЗКА ОПЕРАТОРОВ ДЛЯ std::vector<T> ### */

/* Функция для сложения векторов */
template <typename T>
std::vector<T> operator+(const std::vector<T>& vec1, const std::vector<T>& vec2){
    std::vector<T> pert_vec = vec1;
    for (int i = 0; i < vec1.size(); i++) {
        pert_vec[i] += vec2[i];
    }
    return pert_vec;
}


/* Функция вычитания векторов */
template <typename T>
std::vector<T> operator-(const std::vector<T>& a, const std::vector<T>& b){
    // Проверка на возможность умножения
    if (a.size() != b.size()) {
        std::cout << "Error: size a != size b in substraction vectors." << std::endl;
        exit(1);
    }
    // Создание результирующего вектора
    std:: vector<T> result(a.size(), 0);

    // Умножение матрицы на вектор
    for (int i = 0; i < a.size(); ++i) {
        result[i] += a[i] - b[i];
    }
    return result;

}


/* Операция почленного умножения векторов */
template <typename T>
std::vector<T> operator*(const std::vector<T>& vec1, const std::vector<T>& vec2){
    if (vec1.size() != vec2.size()) {
        std::cout << "Error: vector1 size != vector2 size in operator*." << std::endl;
        exit(1);
    }
    std::vector<T> result(vec1.size(), 0);
    for (int i = 0; i < vec1.size(); i++){
        result[i] = vec1[i] * vec2[i];
    }
    return result;
}


/* Операция умножения вектора на число */
template <typename T>
std::vector<T> operator*(const T& c, const std::vector<T>& vec){
    std::vector<T> result(vec.size(), 0);
    for (int i = 0; i < vec.size(); i++){
        result[i] = vec[i] * c;
    }
    return result;
}


template <typename T>
std::vector<T> operator*(const std::vector<T>& vec, const T& c){
    std::vector<T> result(vec.size(), 0);
    for (int i = 0; i < vec.size(); i++){
        result[i] = vec[i] * c;
    }
    return result;
}

/* Операция деления вектора на число */
template<typename T>
std::vector<T> operator/(const std::vector<T>& vec, const T& c) {
    std::vector<T> result(vec.size(), 0);
    for (int i = 0; i < vec.size(); i++) {
        result[i] = vec[i] / c;
    }
    return result;
}


/* Операция почленного деления векторов */
template <typename T>
std::vector<T> operator/(const std::vector<T>& vec1, const std::vector<T>& vec2){
    if (vec1.size() != vec2.size()) {
        std::cout << "Error: vector1 size != vector2 size in operator*." << std::endl;
        exit(1);
    }
    std::vector<T> result(vec1.size(), 0);
    for (int i = 0; i < vec1.size(); i++){
        result[i] = vec1[i] / vec2[i];
    }
    return result;
}



/* ### ФУНКЦИИ ВЫВОДА ВЕКТОРА ### */



/* Функция вывода вектора на экран */
template <typename T>
void print(const std::vector<T>& vec) {
    for (T value : vec) {
        std::cout << value << ' ';
    }
    std::cout << std::endl;
}


/* Функция вывода обрезанного вектора на экран */
template <typename T>
void print_short(const std::vector<T>& vec, const int& n){

    for (int i = 0; i < n; ++i){
        std::cout << vec[i] << ' ';
    }
    std::cout << "..." << std::endl;
}


/* Функция, которая красиво выводит вектор*/
template<typename T>
void print_vec(const std::vector<T>& vec){
    std::cout << "(" << vec[0];
    for (int i = 1; i < vec.size(); i++){
        std::cout << ", " << vec[i];
    }
    std::cout << ")" << std::endl;
}



/* ### ФУНКЦИИ ДРУГИХ ОПЕРАЦИЙ С ВЕКТОРАМИ ### */



/* Функция для скалярного умножения векторов */
template <typename T>
T dot(const std::vector<T>& vec1, const std::vector<T>& vec2){
    if (vec1.size() != vec2.size()) {
        std::cout << "Error: vector1 size != vector2 size in operator*." << std::endl;
        exit(1);
    }
    T result;
    for (int i = 0; i < vec1.size(); i++){
        result += vec1[i] * vec2[i];
    }
    return result;
}


/* Функция для нормы вектора */
template <typename T>
T norm(const std::vector<T>& vec, const int& p){
    if (vec.empty()) {
        std::cerr << "Error: Empty vector in norm() \n";
        exit(1);
    }

    T result = 0.0;

    // Вычисление нормы
    if (p == 0) {
        // Норма oo
        for (const auto& element : vec) {
            T absElement = abs(element);
            if (absElement > result) {
                result = absElement;
            }
        }
    } else {
        // Общий случай для норм L1, L2 и т.д.
        for (const auto& element : vec) {
            result += pow(abs(element), p);
        }

        result = pow(result, 1.0 / p);
    }

    return result;
}


/* Функция, которая возращает матрицу комбинаций элементов вектора */
template<typename T>
std::vector<std::vector<T>> generateCombinations(const std::vector<T>& vec) {
    int n = vec.size();

    // Вектор для хранения всех комбинаций
    std::vector<std::vector<T>> combinations;

    // Внешний цикл по всем возможным комбинациям
    for (int i = 0; i < (1 << n); ++i) {
        std::vector<T> current(n);

        // Внутренний цикл для каждой позиции вектора
        for (int j = 0; j < n; ++j) {
            current[j] = (i & (1 << j)) ? vec[j] : -vec[j];
        }

        // Добавить текущую комбинацию в вектор
        combinations.push_back(current);
    }

    return combinations;
}


/* Функция, возвращает вектор модулей */
template<typename T>
std::vector<T> vec_abs(const std::vector<T>& vec){
    for (int i = 0; i < vec.size(); i++){
        vec[i] = fabs(vec[i]);
    }
    return vec;
}


/* Функция, возращающая сумму элементов вектора */
template<typename T>
T sum(const std::vector<T>& vec){
    T sum = 0;
    for (int i = 0; i < vec.size(); i++){
        sum += vec[i];
    }
    return sum;
}


/* Функция, сортирующая вектор */
template< typename T>
std::vector<T> sorted(const std::vector<T>& vec_not_sort) {
    std::vector<T> vec(vec_not_sort);
    int n = vec.size();
    for (int i = 0; i < n - 1; ++i) {
        for (int j = 0; j < n - i - 1; ++j) {
            if (vec[j] > vec[j + 1]) {
                // Обмен элементов, если они не упорядочены
                T temp = vec[j];
                vec[j] = vec[j + 1];
                vec[j + 1] = temp;
            }
        }
    }
    return vec;
}


/* Функция, возращающая максимальный по модулю элемент вектора */
template<typename T>
T vec_max(const std::vector<T>& vec){
    int n = vec.size();
    T max = 0;
    for (int i = 0; i < n; i++) {
        if (abs(vec[i]) > max)
            max = abs(vec[i]);
    }
    return max;
}


/* Функция, вычисляющая норму разности векторов */
double sqr(std::vector<double> vec1, std::vector<double> vec2) {
    int m = vec1.size();
    double sum;
    for (int i = 0; i < m; i++) {
        sum = (vec1[i] - vec2[i]) * (vec1[i] - vec2[i]);
    }
    return sum;
}

/* Функция для сдвига вектора на n элементов */
template<typename T>
std::vector<T> shift(const std::vector<T>& vec, int n) {
    std::vector<T> shiftedVec(vec.size()); // Создаем вектор той же длины
    int size = vec.size();

    // Если сдвиг больше длины вектора, находим остаток от деления
    n = n % size;

    // Перемещаем элементы вправо
    if (n >= 0) {
        for (int i = 0; i < size; ++i) {
            shiftedVec[(i + n) % size] = vec[i];
        }
    }
        // Перемещаем элементы влево
    else {
        for (int i = 0; i < size; ++i) {
            shiftedVec[(size + i + n) % size] = vec[i];
        }
    }

    return shiftedVec;
}
