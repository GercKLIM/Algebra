/*
 *  РЕАЛИЗАЦИЯ АЛГЕБРЫ МАТРИЦ
 *
 * */

#include "../algebra.h"



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


/* Операция для умножения  числа на матрицу */
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








/* Операция умножения матрицы на вектор */
template<typename T>
std::vector<T> operator*(const Matrix<T>& matrix, const Vector<T>& vec) {
    if (matrix.data[0].size() != vec.size()) {
        throw std::invalid_argument("Matrix and vector dimensions must match for multiplication");
    }

    Vector<T> result(matrix.data.size(), 0);

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
Matrix<T> Matrix<T>::round(const double& eps){
    Matrix<T> roundA(data);
    int size = data.size();

    for (int i = 0; i < size; ++i) {
        for (int j = 0; j < size; ++j) {
            roundA[i][j] = (round(data[i][j]) >= 0)? round(abs(data[i][j]) * (1 / eps)) / (1 / eps): -1 * round(abs(data[i][j]) * (1 / eps)) / (1 / eps);
        }
    }
    return roundA;
}


/* Функция для вычисления нормы матрицы */
template <typename T>
T Matrix<T>::norm(const int& p) {
    // Проверка на пустую матрицу
    if (data.empty() || data[0].empty()) {
        std::cout << "Error: Empty matrix in norm()\n" << std::endl;
        exit(1);
    }

    int rows = data.size();
    int cols = data[0].size();

    T result = 0.0;

    // Вычисление нормы матрицы
    if (p == 0) {
        // Норма матрицы Чебышева (максимальное значение по модулю в строке)
        for (int i = 0; i < rows; ++i) {
            T rowSum = 0.0;
            for (int j = 0; j < cols; ++j) {
                rowSum += abs(data[i][j]);
            }
            if (rowSum > result) {
                result = rowSum;
            }
        }
    } else {
        // Общий случай для норм матрицы (Фробениуса и др.)
        for (int j = 0; j < cols; ++j) {
            T colSum = 0.0;
            for (T i = 0; i < rows; ++i) {
                colSum += pow(abs(data[i][j]), p);
            }
            result += pow(colSum, 1.0 / p);
        }

        result = pow(result, 1.0 / p);
    }
    return result;
}


/* Функция поворота матрицы вправо */
template <typename T>
std::vector<std::vector<T>> RotateRight(const std::vector<std::vector<T>>& A){

    Matrix<T> A_rotate(A.size(), Vector<T>(A.size(), 0));

    for (int i = 0; i < A.size(); ++i) {
        for (int j = 0; j < A.size(); ++j) {
            A_rotate[A.size() - 1 - j][i] = A[i][j];
        }
    }

    return A_rotate;

}


/* Функция поворота матрицы влево */
template <typename T>
std::vector<std::vector<T>> RotateLeft(const std::vector<std::vector<T>>& A){

    Matrix<T> A_rotate(A.size(), Vector<T>(A.size(), 0));

    for (int i = 0; i < A.size(); ++i) {
        for (int j = 0; j < A.size(); ++j) {
            A_rotate[j][A.size() - 1 - i] = A[i][j];
        }
    }

    return A_rotate;
}


// Функция для создания единичной матрицы размера n x n
template <typename T>
Matrix<T> Matrix<T>::create_identity_matrix(const int& n) {
    Matrix<T> identity(n, Vector<T>(n, 0));
    for (int i = 0; i < n; i++) {
        identity[i][i] = 1;
    }
    return identity;
}


// Функция для обратной матрицы с проверкой на вырожденность
//template <typename T>
//std::vector<std::vector<T>> inverseMatrix(const std::vector<std::vector<T>>& A, const T& eps) {
//    std::vector<std::vector<T>> E = create_identity_matrix<T>(A.size());
//    std::vector<std::vector<T>> E_rotate = RotateLeft(E);
//    std::vector<T> e(A.size());
//    std::vector<std::vector<T>> X(A.size(), std::vector<T>(A.size(), 0));
//
//
//    for (int i = 0; i < A.size(); i++){
//        e = E_rotate[i];
//        X[i] = method_Gaussa(A, e, eps);
//
//    }
//    std::vector<std::vector<T>> A_inv = RotateLeft(X);
//    return A_inv;
//}


// Функция для обратной матрицы с проверкой на вырожденность с максимальной точностью
//template <typename T>
//std::vector<std::vector<T>> inverseMatrix(const std::vector<std::vector<T>>& A){
//    T eps = std::numeric_limits<T>::epsilon();
//    return inverseMatrix(A, eps);
//}



/* Функция транспонирования матрицы */
template <typename T>
std::vector<std::vector<T>> transpose(const std::vector<std::vector<T>>& A) {
    int rows = A.size();
    int cols = A[0].size();
    std::vector<std::vector<T>> result(cols, std::vector<T>(rows));

    for (int i = 0; i < rows; i++) {
        for (int j = 0; j < cols; j++) {
            result[j][i] = A[i][j];
        }
    }

    return result;
}


// Функция обрезки матрицы снизу и справа
template <typename T>
std::vector<std::vector<T>> crop_matrix(const std::vector<std::vector<T>>& A, const int& k){

    int n = A.size();
    std::vector<std::vector<T>> A_crop(n - k, std::vector<T>(n - k, 0));
    for (int i = 0; i < (n - k); i++){
        for (int j = 0; j < (n - k); j++){
            A_crop[i][j] = A[i][j];
        }
    }

    return A_crop;
}