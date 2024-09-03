/*
 *  РЕАЛИЗАЦИЯ АЛГЕБРЫ МАТРИЦ
 *
 * */

#include "../Algebra.h"



/* ### ФУНКЦИИ ВЫВОДА МАТРИЦ НА ЭКРАН ### */



/* Функция вывода матрицы на экран */
template <typename T>
void Matrix<T>::print() {
    for (int i = 0; i < data.size(); i++)
        std::cout << data[i] << std::endl;
}



/* Функция вывода матрицы на экран */
template <typename T>
void Matrix<T>::print_matr() {
    std::cout << "{" << "{" << data[0] << "}";
    for (int i = 1; i < data.size(); i++) {
        std::cout << ", {" << data[i] << "}" /*<< std::endl*/;
    }
    std::cout << "}";
}



/* ### ПЕРЕГРУЗКИ ОПЕРАТОРОВ ДЛЯ Matrix ### */

/* Операция сравнения */
template<typename T>
bool operator==(const Matrix<T> A1, const Matrix<T> A2){

    if (A1.size() != A2.size()){
        return false;
    }

    for (int i = 0; i < A1.size(); i++){
        if (not(A1[i] == A2[i])) {
            return false;
        }
    }

    return true;
}

/* Операция поэлементного сложения матриц */
template<typename T>
Matrix<T> operator+(Matrix<T>& A, Matrix<T>& B) {
    if (A.data.size() != B.data.size() || A.data[0].size() != B.data[0].size()) {
        throw std::invalid_argument("Matrix dimensions must match for addition");
    }

    Matrix<T> result(A.data.size(), A.data[0].size(), 0);

    for (size_t i = 0; i < A.data.size(); ++i) {
        for (size_t j = 0; j < A.data[i].size(); ++j) {
            result.data[i][j] = A.data[i][j] + B.data[i][j];
        }
    }

    return result;
}


/* Операция поэлементного сложения сonst матриц */
template<typename T>
Matrix<T> operator+(const Matrix<T>& A, const Matrix<T>& B) {
    if (A.data.size() != B.data.size() || A.data[0].size() != B.data[0].size()) {
        throw std::invalid_argument("Matrix dimensions must match for addition");
    }

    Matrix<T> result(A.data.size(), A.data[0].size(), 0);

    for (size_t i = 0; i < A.data.size(); ++i) {
        for (size_t j = 0; j < A.data[i].size(); ++j) {
            result.data[i][j] = A.data[i][j] + B.data[i][j];
        }
    }

    return result;
}


/* Операция поэлементного вычитания матриц */
template<typename T>
Matrix<T> operator-(Matrix<T>& A, Matrix<T>& B) {
    if (A.data.size() != B.data.size() || A.data[0].size() != B.data[0].size()) {
        throw std::invalid_argument("Matrix dimensions must match for addition");
    }

    Matrix<T> result(A.data.size(), A.data[0].size(), 0);

    for (size_t i = 0; i < A.data.size(); ++i) {
        for (size_t j = 0; j < A.data[i].size(); ++j) {
            result.data[i][j] = A.data[i][j] - B.data[i][j];
        }
    }

    return result;
}


/* Операция поэлементного вычитания сonst матриц */
template <typename T>
Matrix<T> operator-(const Matrix<T>& A, const Matrix<T>& B) {
    if (A.data.size() != B.data.size() || A.data[0].size() != B.data[0].size()) {
        throw std::invalid_argument("Matrix dimensions must match for addition");
    }

    Matrix<T> result(A.data.size(), A.data[0].size(), 0);

    for (size_t i = 0; i < A.data.size(); ++i) {
        for (size_t j = 0; j < A.data[i].size(); ++j) {
            result.data[i][j] = A.data[i][j] - B.data[i][j];
        }
    }

    return result;
}



/* Матричное умножение */
template<typename T>
Matrix<T> operator*(Matrix<T>& A, Matrix<T>& B) {
    size_t rowsA = A.size();
    size_t colsA = A[0].size();
    size_t rowsB = B.size();
    size_t colsB = B[0].size();

    // Проверка на возможность умножения матриц
    if (colsA != rowsB) {
        throw std::invalid_argument("Matrix dimensions do not match for multiplication.");
    }

    Matrix<T> result(rowsA, colsB, T(0));  // Создание матрицы-результата

    for (size_t i = 0; i < rowsA; ++i) {
        for (size_t j = 0; j < colsB; ++j) {
            for (size_t k = 0; k < colsA; ++k) {
                result[i][j] += A[i][k] * B[k][j];
            }
        }
    }

    return result;
}


/* Матричное умножение */
template <typename T>
Matrix<T> operator*(const Matrix<T>& A, const Matrix<T>& B) {
    if (A.data[0].size() != B.data.size()) {
        throw std::invalid_argument("Matrix dimensions must match for multiplication");
    }

    Matrix<T> result(A.data.size(), B.data[0].size(), 0);

    for (size_t i = 0; i < A.data.size(); ++i) {
        for (size_t j = 0; j < B.data[0].size(); ++j) {
            for (size_t k = 0; k < A.data[0].size(); ++k) {
                result.data[i][j] += A.data[i][k] * B.data[k][j];
            }
        }
    }

    return result;
}


/* Определение оператора отрицания для матрицы */
template <typename Y>
Matrix<Y> operator-(Matrix<Y>& A) {

    Matrix<Y> result(A.size(), Vector<Y>(A[0].size(), 0));

    for (int i = 0; i < A.size(); i++){
        for (int j = 0; j < A[0].size(); j++) {
            result[i][j] = -1 * A[i][j];
        }
    }
    return result;
}



/* Операция для умножения матрицы на число */
template <typename T>
Matrix<T> operator*(const Matrix<T>& A, const T& scalar){
    // Создание результирующей матрицы с теми же размерами
    Matrix<T> result(A.size(), Vector<T>(A[0].size(), 0));

    // Умножение каждого элемента матрицы на число
    for (size_t i = 0; i < A.size(); ++i) {
        for (size_t j = 0; j < A[0].size(); ++j) {
            result[i][j] = A[i][j] * scalar;
        }
    }

    return result;
}


/* Операция для умножения числа на матрицу */
template <typename T>
Matrix<T> operator*(const T& scalar, const Matrix<T>& A){
    // Создание результирующей матрицы с теми же размерами
    Matrix<T> result(A.size(), Vector<T>(A[0].size(), 0));

    // Умножение каждого элемента матрицы на число
    for (size_t i = 0; i < A.size(); ++i) {
        for (size_t j = 0; j < A[0].size(); ++j) {
            result[i][j] = A[i][j] * scalar;
        }
    }

    return result;
}


/* Операция умножения матрицы на вектор справа */
template<typename T>
Vector<T> operator*(Matrix<T>& matrix, Vector<T>& vec) {
    if (matrix[0].size() != vec.size()) {
        throw std::invalid_argument("Matrix and vector dimensions must match for multiplication");
    }

    Vector<T> result(matrix.size(), 0);

    for (size_t i = 0; i < matrix.size(); ++i) {
        for (size_t j = 0; j < matrix[0].size(); ++j) {
            result[i] += matrix[i][j] * vec[j];
        }
    }

    return result;
}

/* Операция умножения матрицы на вектор слева */
template<typename Y>
Vector<Y> operator*(Vector<Y>& vec, Matrix<Y>& matrix) {
    if (matrix.data[0].size() != vec.size()) {
        throw std::invalid_argument("Matrix and vector dimensions must match for multiplication");
    }

    Vector<Y> result(matrix.data.size(), 0);

    for (size_t i = 0; i < matrix.data.size(); ++i) {
        for (size_t j = 0; j < matrix.data[0].size(); ++j) {
            result[i] += matrix.data[i][j] * vec[j];
        }
    }

    return result;
}


// Определение оператора отрицания для матрицы
template <typename T>
Matrix<T> operator-(const Matrix<T>& matrix) {
    int rows = matrix.size();
    int cols = matrix[0].size();
    std::vector<std::vector<T>> result(rows, std::vector<T>(cols, 0));
    for (int i = 0; i < rows; ++i) {
        for (int j = 0; j < cols; ++j) {
            result[i][j] = -matrix[i][j];
        }
    }
    return result;
}


/* Реализация оператора потока вывода для const матрицы */
template <typename T>
std::ostream& operator<<(std::ostream& os, const Matrix<T>& matrix) {
    for (const auto& row : matrix.data) {
        os << row << std::endl;
    }
    return os;
}






/* ### ФУНКЦИИ ДРУГИХ ОПЕРАЦИЙ С МАТРИЦАМИ ### */


/* Функция преобразование матрицы в std::vector<Vector<T>> */
template <typename T>
std::vector<Vector<T>> Matrix<T>::to_std(){

    return data;
}


/* Функция преобразование матрицы в std::vector<std::vector<T>> */
template <typename T>
std::vector<std::vector<T>> Matrix<T>::to_stds(){

    std::vector<std::vector<T>> std_data(data.size(), std::vector<T>(data[0].size()));
    for (int i = 0; i < data.size(); i++){
        for (int j = 0; j < data[0].size(); j++){
            std_data[i][j] = data[i][j];
        }
    }

    return std_data;
}


/* Функция транспонирования матрицы */
template <typename T>
Matrix<T> Matrix<T>::transpose() {
    int rows = data.size();
    int cols = data.size();
    std::vector<std::vector<T>> result(cols, std::vector<T>(rows));

    for (int i = 0; i < rows; i++) {
        for (int j = 0; j < cols; j++) {
            result[j][i] = data[i][j];
        }
    }

    return result;
}


/* Функция для поэлементного умножения матриц */
template <typename T>
std::vector<std::vector<T>> Multyply(const std::vector<std::vector<T>>& A, const std::vector<std::vector<T>>& B){
    int m = A.size();    // Количество строк в матрице A
    int n = A[0].size(); // Количество столбцов в матрице A
    int p = B[0].size(); // Количество столбцов в матрице B

    if (n != B.size()) {
        printf("Error: impossible multiply matrix");
        exit(1);
    }

    std::vector<std::vector<T>> result(m, std::vector<T>(p, 0.0));

    for (int i = 0; i < m; i++) {
        for (int j = 0; j < p; j++) {
            result[i][j] = A[i][j] * B[i][j];
        }
    }
    return result;
}


/* Функция округления чисел в матрицах */
template <typename T>
Matrix<T> Matrix<T>::round(const double& eps) {


    int rows = data.size();
    int cols = data[0].size();

    for (int i = 0; i < rows; ++i) {
        for (int j = 0; j < cols; ++j) {

            // Округляем значение
            data[i][j] = std::round(data[i][j] / eps) * eps;
        }
    }

    Matrix<T> roundA(data);
    return roundA;
}


/* Функция для вычисления нормы матрицы */
template <typename T>
T Matrix<T>::norm(const int& p) {
    // Проверка на пустую матрицу
    if (data.size() == 0 || data[0].size() == 0) {
        std::cerr << "Error: Empty matrix in norm()\n";
        exit(1);
    }

    int rows = data.size();
    int cols = data[0].size();

    T result = 0.0;

    if (p == 0) {

        // Норма Чебышева (максимум суммы абсолютных значений в столбце)
        for (int j = 0; j < cols; ++j) {
            T colSum = 0.0;
            for (int i = 0; i < rows; ++i) {
                colSum += std::abs(data[i][j]);
            }
            result = std::max(result, colSum);
        }
    } else {

        // Общая p-норма (например, Фробениусова норма для p=2)
        T sum = 0.0;
        for (int i = 0; i < rows; ++i) {
            for (int j = 0; j < cols; ++j) {
                sum += std::pow(std::abs(data[i][j]), p);
            }
        }
        result = std::pow(sum, 1.0 / p);
    }

    return result;
}


/* Функция поворота матрицы вправо */
template <typename T>
Matrix<T> Matrix<T>::rotateRight(){

    Matrix<T> data_rotate(data.size(), Vector<T>(data.size(), 0));

    for (int i = 0; i < data.size(); ++i) {
        for (int j = 0; j < data.size(); ++j) {
            data_rotate[data.size() - 1 - j][i] = data[i][j];
        }
    }

    return data_rotate;

}


/* Функция поворота матрицы влево */
template <typename T>
Matrix<T> Matrix<T>::rotateLeft(){

    Matrix<T> data_rotate(data.size(), Vector<T>(data.size(), 0));

    for (int i = 0; i < data.size(); ++i) {
        for (int j = 0; j < data.size(); ++j) {
            data_rotate[j][data.size() - 1 - i] = data[i][j];
        }
    }

    return data_rotate;
}


// Функция для создания единичной матрицы размера n x n
template <typename T>
Matrix<T> Matrix<T>::create_identity_matrix(const int& size) {
    std::vector<Vector<T>> E(size, Vector<T>(size, 0));
    for (int i = 0; i < size; i++) {
        E[i][i] = 1;
    }
    data = E;
    return E;
}


template <typename T>
Vector<T> method_Gaussa(const Matrix<T>& A, const Vector<T>& b, const T& eps) {
    int n = A.size();
    Matrix<T> augmentedMatrix(n, Vector<T>(n + 1)); // Расширенная матрица

    // Формируем расширенную матрицу [A|b]
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            augmentedMatrix[i][j] = A[i][j];
        }
        augmentedMatrix[i][n] = b[i];
    }

    // Прямой ход метода Гаусса
    for (int i = 0; i < n; ++i) {
        // Поиск максимального элемента в столбце для выбора главного элемента
        int maxRow = i;
        for (int k = i + 1; k < n; ++k) {
            if (std::abs(augmentedMatrix[k][i]) > std::abs(augmentedMatrix[maxRow][i])) {
                maxRow = k;
            }
        }

        // Проверка на вырожденность
        if (std::abs(augmentedMatrix[maxRow][i]) < eps) {
            throw std::runtime_error("Matrix is singular or near-singular");
        }

        // Перестановка строк
        std::swap(augmentedMatrix[i], augmentedMatrix[maxRow]);

        // Приводим главный элемент в 1, делим всю строку на него
        T mainElement = augmentedMatrix[i][i];
        for (int j = i; j <= n; ++j) {
            augmentedMatrix[i][j] /= mainElement;
        }

        // Обнуление элементов ниже главного
        for (int k = i + 1; k < n; ++k) {
            T factor = augmentedMatrix[k][i];
            for (int j = i; j <= n; ++j) {
                augmentedMatrix[k][j] -= factor * augmentedMatrix[i][j];
            }
        }
    }

    // Обратный ход метода Гаусса
    Vector<T> x(n);
    for (int i = n - 1; i >= 0; --i) {
        x[i] = augmentedMatrix[i][n];
        for (int j = i + 1; j < n; ++j) {
            x[i] -= augmentedMatrix[i][j] * x[j];
        }
    }

    return x;
}


// Функция для обратной матрицы с проверкой на вырожденность
template <typename T>
Matrix<T> Matrix<T>::inverse(const T& eps) {
    Matrix<T> E = create_identity_matrix(data.size());
    Matrix<T> E_rotate = RotateLeft(E);
    Vector<T> e(data.size());
    Matrix<T> X(data.size(), std::vector<T>(data.size(), 0));
    Matrix<T> Data(data);

    for (int i = 0; i < data.size(); i++){
        e = E_rotate[i];
        X[i] = method_Gaussa(Data, e, eps);

    }
    Matrix<T> A_inv = RotateLeft(X);
    return A_inv;
}


// Функция для обратной матрицы с проверкой на вырожденность с максимальной точностью
//template <typename T>
//Matrix<T> Matrix<T>::inverse(){
//
//    // миниму разрядности
//    T eps = std::numeric_limits<T>::epsilon();
//
//    Matrix<T> E = create_identity_matrix(data.size());
//    Matrix<T> E_rotate = E.rotateLeft();
//    Vector<T> e(data.size());
//    Matrix<T> X(data.size(), std::vector<T>(data.size(), 0));
//    Matrix<T> Data(data);
//
//    for (int i = 0; i < data.size(); i++){
//        e = E_rotate[i];
//        X[i] = method_Gaussa(Data, e, eps);
//
//    }
//    Matrix<T> A_inv = X.rotateLeft();
//    return A_inv;
//}

// Функция для обратной матрицы с проверкой на вырожденность с максимальной точностью
template <typename T>
Matrix<T> Matrix<T>::inverse(){
    size_t n = data.size();
    if (n == 0 || data[0].size() != n) {
        throw std::invalid_argument("Matrix must be square to compute inverse.");
    }

    // Создаем копию матрицы и единичную матрицу
    Matrix<T> A_copy = *this;
    Matrix<T> I(n, std::vector<T>(n, 0));
    for (size_t i = 0; i < n; ++i) {
        I[i][i] = 1;
    }

    // Прямой ход алгоритма Гаусса
    for (size_t i = 0; i < n; ++i) {

        // Найти максимальный элемент в текущем столбце
        size_t max_row = i;
        for (size_t k = i + 1; k < n; ++k) {
            if (std::abs(A_copy[k][i]) > std::abs(A_copy[max_row][i])) {
                max_row = k;
            }
        }

        // Переставляем строки
        if (i != max_row) {
            std::swap(A_copy[i], A_copy[max_row]);
            std::swap(I[i], I[max_row]);
        }

        // Проверяем на вырожденность
        if (std::abs(A_copy[i][i]) < 1e-10) {
            throw std::runtime_error("Matrix is singular or nearly singular");
        }

        // Приведение к диагональному виду
        T diag_val = A_copy[i][i];
        for (size_t j = 0; j < n; ++j) {
            A_copy[i][j] /= diag_val;
            I[i][j] /= diag_val;
        }

        // Приведение остальных строк к нулю в текущем столбце
        for (size_t k = 0; k < n; ++k) {
            if (k != i) {
                T factor = A_copy[k][i];
                for (size_t j = 0; j < n; ++j) {
                    A_copy[k][j] -= factor * A_copy[i][j];
                    I[k][j] -= factor * I[i][j];
                }
            }
        }
    }
    return I;
}







// Функция обрезки матрицы снизу и справа
template <typename T>
Matrix<T> Matrix<T>::crop(const int& new_size){

    int n = data.size();
    if (n < new_size) {
        return data;
    }

    Matrix<T> data_crop(new_size, std::vector<T>(new_size, 0));
    for (int i = 0; i < (new_size); i++){
        for (int j = 0; j < (new_size); j++){
            data_crop[i][j] = data[i][j];
        }
    }

    return data_crop;
}


template<typename T>
T Matrix<T>::det() const {
    int n = rows();
    if (n != cols()) {
        throw std::invalid_argument("Matrix must be square to compute determinant.");
    }

    // Детераминант для матрицы 2x2
    if (n == 2) {
        return data[0][0] * data[1][1] - data[0][1] * data[1][0];
    }

    // Детераминант для матрицы 3x3
    if (n == 3) {
        return data[0][0] * (data[1][1] * data[2][2] - data[1][2] * data[2][1]) -
               data[0][1] * (data[1][0] * data[2][2] - data[1][2] * data[2][0]) +
               data[0][2] * (data[1][0] * data[2][1] - data[1][1] * data[2][0]);
    }

    // Для матриц размером больше 3x3 используем рекурсивный подход
    T det = 0;
    Matrix<T> minor(n - 1, n - 1);

    for (int i = 0; i < n; ++i) {
        // Создаем минор матрицы
        for (int j = 1; j < n; ++j) {
            int colIndex = 0;
            for (int k = 0; k < n; ++k) {
                if (k == i) continue;
                minor(j - 1, colIndex) = data[j][k];
                ++colIndex;
            }
        }

        // Вычисляем детерминант с учетом знака
        det += (i % 2 == 0 ? 1 : -1) * data[0][i] * minor.det();
    }

    return det;
}