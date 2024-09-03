/*
 *  РЕАЛИЗАЦИЯ АЛГЕБРЫ ВЕКТОРОВ
 *
 * */

#include "../Algebra.h"



/* ### ПЕРЕГРУЗКА ОПЕРАТОРОВ ДЛЯ std::vector<T> ### */


template<typename Y>
bool operator==(const Vector<Y>& v1, const Vector<Y>& v2) {
    // Сначала сравниваем размеры
    if (v1.size() != v2.size()) {
        return false;
    }

    // Затем сравниваем элементы
    for (int i = 0; i < v1.size(); i++){
        if (v1[i] != v2[i]){
            return false;
        }
    }
    return true;
}


/* Операция сложения векторов */
template<typename T>
Vector<T> operator+(Vector<T>& v1, Vector<T>& v2) {
    if (v1.size() != v2.size()) {
        std::cout << "Error with operator+ Vector: different size" << std::endl;
    }

    Vector<T> result(v1.size());
    for (size_t i = 0; i < v1.size(); ++i) {
        result[i] = v1[i] + v2[i];
    }

    return result;
}


template<typename T>
Vector<T> operator+(const Vector<T>& v1, const Vector<T>& v2) {
    if (v1.size() != v2.size()) {
        std::cout << "Error with operator+ Vector: different size" << std::endl;
    }

    Vector<T> result(v1.size());
    for (size_t i = 0; i < v1.size(); ++i) {
        result[i] = v1[i] + v2[i];
    }

    return result;
}


/* Функция вычитания векторов */
template<typename T>
Vector<T> operator-(Vector<T>& v1, Vector<T>& v2) {
    if (v1.size() != v2.size()) {
        std::cout << "Error with operator- Vector: different size" << std::endl;
    }

    Vector<T> result(v1.size());
    for (size_t i = 0; i < v1.size(); ++i) {
        result[i] = v1[i] - v2[i];
    }

    return result;
}


/* Функция вычитания const векторов */
template<typename T>
Vector<T> operator-(const Vector<T>& v1, const Vector<T>& v2) {
    if (v1.size() != v2.size()) {
        std::cout << "Error with operator- Vector: different size" << std::endl;
    }

    Vector<T> result(v1.size());
    for (size_t i = 0; i < v1.size(); ++i) {
        result[i] = v1[i] - v2[i];
    }

    return result;
}


/* Операция почленного умножения векторов */
template<typename T>
Vector<T> operator*(Vector<T>& v1, Vector<T>& v2) {
    if (v1.size() != v2.size()) {
        std::cout << "Error with operator- Vector: different size" << std::endl;
    }

    Vector<T> result(v1.size());
    for (size_t i = 0; i < v1.size(); ++i) {
        result[i] = v1[i] * v2[i];
    }

    return result;
}


/* Операция почленного умножения const векторов */
template<typename T>
Vector<T> operator*(const Vector<T>& v1, const Vector<T>& v2) {
    if (v1.size() != v2.size()) {
        std::cout << "Error with operator- Vector: different size" << std::endl;
    }

    Vector<T> result(v1.size());
    for (size_t i = 0; i < v1.size(); ++i) {
        result[i] = v1[i] * v2[i];
    }

    return result;
}


/* Операция умножения вектора на число */
template <typename T>
Vector<T> operator*(T& c, Vector<T>& vec){
    Vector<T> result(vec.size(), 0);
    for (int i = 0; i < vec.size(); i++){
        result[i] = vec[i] * c;
    }
    return result;
}


/* Операция умножения const вектора на const число */
template <typename T>
Vector<T> operator*(const T& c, const Vector<T>& vec){
    Vector<T> result(vec.size(), 0);
    for (int i = 0; i < vec.size(); i++){
        result[i] = vec[i] * c;
    }
    return result;
}


/* Операция умножения const вектора на число */
template <typename T>
Vector<T> operator*(T& c, const Vector<T>& vec){
    Vector<T> result(vec.size(), 0);
    for (int i = 0; i < vec.size(); i++){
        result[i] = vec[i] * c;
    }
    return result;
}


/* Операция умножения вектора на const число */
template <typename T>
Vector<T> operator*(const T& c, Vector<T>& vec){
    Vector<T> result(vec.size(), 0);
    for (int i = 0; i < vec.size(); i++){
        result[i] = vec[i] * c;
    }
    return result;
}


/* Операция деления вектора на число */
template <typename T>
Vector<T> operator/(Vector<T>& vec, T& c){
    Vector<T> result(vec.size(), 0);
    for (int i = 0; i < vec.size(); i++){
        result[i] = vec[i] / c;
    }
    return result;
}


/* Операция деления сonst вектора на const число */
template <typename T>
Vector<T> operator/(const Vector<T>& vec, const T& c){
    Vector<T> result(vec.size(), 0);
    for (int i = 0; i < vec.size(); i++){
        result[i] = vec[i] / c;
    }
    return result;
}


/* Операция деления вектора на const число */
template <typename T>
Vector<T> operator/(Vector<T>& vec, const T& c){
    Vector<T> result(vec.size(), 0);
    for (int i = 0; i < vec.size(); i++){
        result[i] = vec[i] / c;
    }
    return result;
}


/* Операция деления const вектора на число */
template <typename T>
Vector<T> operator/(const Vector<T>& vec, T& c){
    Vector<T> result(vec.size(), 0);
    for (int i = 0; i < vec.size(); i++){
        result[i] = vec[i] / c;
    }
    return result;
}


/* Операция почленного деления векторов */
template<typename T>
Vector<T> operator/(Vector<T>& v1, Vector<T>& v2) {
    if (v1.size() != v2.size()) {
        std::cout << "Error with operator- Vector: different size" << std::endl;
    }

    Vector<T> result(v1.size());
    for (size_t i = 0; i < v1.size(); ++i) {
        result[i] = v1[i] / v2[i];
    }

    return result;
}


/* Операция почленного деления векторов */
template<typename T>
Vector<T> operator/(const Vector<T>& v1, const Vector<T>& v2) {
    if (v1.size() != v2.size()) {
        std::cout << "Error with operator- Vector: different size" << std::endl;
    }

    Vector<T> result(v1.size());
    for (size_t i = 0; i < v1.size(); ++i) {
        result[i] = v1[i] / v2[i];
    }

    return result;
}


/* Операция потокового вывода */
template<typename T>
std::ostream& operator<<(std::ostream& os, Vector<T>& v) {
    for (size_t i = 0; i < v.size(); ++i) {
        os << v[i];
        if (i < v.size() - 1) {
            os << " ";
        }
    };
    return os;
}

/* Операция потокового вывода для const */
template<typename T>
std::ostream& operator<<(std::ostream& os, const Vector<T>& v) {
    for (size_t i = 0; i < v.size(); ++i) {
        os << v[i];
        if (i < v.size() - 1) {
            os << " ";
        }
    }
    return os;
}



/* ### ФУНКЦИИ ВЫВОДА ВЕКТОРА ### */



/* Функция вывода вектора на экран */
template <typename T>
void Vector<T>::print() {
    std::cout << data << std::endl;
}


/* Функция вывода обрезанного вектора на экран */
template <typename T>
void Vector<T>::print(const int& n){

    for (int i = 0; i < n; ++i){
        std::cout << data[i] << ' ';
    }
    std::cout << /*"..." <<*/ std::endl;
}


/* Функция, которая красиво выводит вектор*/
template<typename T>
void Vector<T>::print_vec(){
    std::cout << "(" << data[0];
    for (int i = 1; i < data.size(); i++){
        std::cout << ", " << data[i];
    }
    std::cout << ")" << std::endl;
}


/* Функция, которая красиво выводит вектор*/
template<typename T>
void Vector<T>::print_vec(const int& n){
    std::cout << "(" << data[0];
    for (int i = 1; i < n; i++){
        std::cout << ", " << data[i];
    }
    std::cout << ")" << std::endl;
}


/* ### ФУНКЦИИ ДРУГИХ ОПЕРАЦИЙ С ВЕКТОРАМИ ### */

template <typename T>
std::vector<T> Vector<T>::to_std(){
    return data;
}

/* Функция для скалярного умножения векторов */
template <typename T>
T Vector<T>::dot(const Vector<T>& vec2){
    if (data.size() != vec2.size()) {
        std::cout << "Error: vector1 size != vector2 size in dot." << std::endl;
        exit(1);
    }
    T result;
    for (int i = 0; i < vec2.size(); i++){
        result += data[i] * vec2[i];
    }
    return result;
}


/* Функция для скалярного умножения векторов */
template<typename T>
T dot(const Vector<T>& vec1, const Vector<T>& vec2){
    if (vec1.size() != vec2.size()) {
        std::cout << "Error: vector1 size != vector2 size in dot." << std::endl;
        exit(1);
    }
    T result;
    for (int i = 0; i < vec2.size(); i++){
        result += vec1[i] * vec2[i];
    }
    return result;
}


/* Функция для нормы вектора */
template <typename T>
T Vector<T>::norm(const int& p){
    if (data.empty()) {
        std::cout << "Error: Empty vector in norm() \n";
        //exit(1);
    }

    T result = 0.0;

    // Вычисление нормы
    if (p == 0) {
        // Норма oo
        for (const auto& element : data) {
            T absElement = fabs(element);
            if (absElement > result) {
                result = absElement;
            }
        }
    } else {
        // Общий случай для норм L1, L2 и т.д.
        for (const auto& element : data) {
            result += pow(fabs(element), p);
        }

        result = pow(result, 1.0 / p);
    }

    return result;
}


/* Функция, которая возращает матрицу комбинаций элементов вектора */
//template<typename T>
//std::vector<std::vector<T>> generateCombinations(const std::vector<T>& vec) {
//    int n = vec.size();
//
//    // Вектор для хранения всех комбинаций
//    std::vector<std::vector<T>> combinations;
//
//    // Внешний цикл по всем возможным комбинациям
//    for (int i = 0; i < (1 << n); ++i) {
//        std::vector<T> current(n);
//
//        // Внутренний цикл для каждой позиции вектора
//        for (int j = 0; j < n; ++j) {
//            current[j] = (i & (1 << j)) ? vec[j] : -vec[j];
//        }
//
//        // Добавить текущую комбинацию в вектор
//        combinations.push_back(current);
//    }
//
//    return combinations;
//}


/* Функция, возвращает вектор модулей */
template<typename T>
Vector<T> Vector<T>::abs(){
    for (int i = 0; i < data.size(); i++){
        data[i] = fabs(data[i]);
    }
    return data;
}


/* Функция, возращающая сумму элементов вектора */
template<typename T>
T Vector<T>::sum(){
    T sum = 0;
    for (int i = 0; i < data.size(); i++){
        sum += data[i];
    }
    return sum;
}


/* Функция, сортирующая вектор */
template< typename T>
Vector<T> Vector<T>::sort() {
    Vector<T> vec(data);
    int n = data.size();
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
T Vector<T>::max(){
    int n = data.size();
    T max = 0;
    for (int i = 0; i < n; i++) {
        if (std::abs(data[i]) > max)
            max = std::abs(data[i]);
    }
    return max;
}


/* Функция для сдвига вектора на n элементов */
template<typename T>
Vector<T> Vector<T>::shift(int n) {
    Vector<T> shiftedVec(data.size(), 0); // Создаем вектор той же длины
    int size = data.size();

    // Если сдвиг больше длины вектора, находим остаток от деления
    n = n % size;

    // Перемещаем элементы вправо
    if (n >= 0) {
        for (int i = 0; i < size; ++i) {
            shiftedVec[(i + n) % size] = data[i];
        }
    }
        // Перемещаем элементы влево
    else {
        for (int i = 0; i < size; ++i) {
            shiftedVec[(size + i + n) % size] = data[i];
        }
    }

    return shiftedVec;
}
