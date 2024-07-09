/*
 *  РЕАЛИЗАЦИЯ АЛГЕБРЫ МАТРИЦ
 *
 * */

/* Функция вывода матрицы на экран */
template <typename T>
void print(const std::vector<std::vector<T>>& matrix) {
    for (std::vector<T> row : matrix) {
        for (T value: row) {
            std::cout << value << ' ';
        }
        std::cout << '\n';
    }
    std::cout << std::endl;
}

/* Операция для умножения матрицы на число */
//template <typename T>
//vector<vector<T>> operator*(const vector<vector<T>>& A, const T& scalar){
//    // Создание результирующей матрицы с теми же размерами
//    vector<vector<T>> result(A.size(), vector<T>(A[0].size(), 0));
//
//    // Умножение каждого элемента матрицы на число
//    for (size_t i = 0; i < A.size(); ++i) {
//        for (size_t j = 0; j < A[0].size(); ++j) {
//            result[i][j] = A[i][j] * scalar;
//        }
//    }
//
//    return result;
//}


/* Операция для умножения  числа на матрицу */
//template <typename T>
//vector<vector<T>> operator*(const T& scalar, const vector<vector<T>>& A){
//    // Создание результирующей матрицы с теми же размерами
//    vector<vector<T>> result(A.size(), vector<T>(A[0].size(), 0));
//
//    // Умножение каждого элемента матрицы на число
//    for (size_t i = 0; i < A.size(); ++i) {
//        for (size_t j = 0; j < A[0].size(); ++j) {
//            result[i][j] = A[i][j] * scalar;
//        }
//    }
//
//    return result;
//}


/* Операция поэлементного сложения матриц */
template <typename T>
std::vector<std::vector<T>> operator+(const std::vector<std::vector<T>>& A, const std::vector<std::vector<T>>& B){
    // Проверка на совпадение размеров матриц
    if (A.size() != B.size() || A[0].size() != B[0].size()) {
        std::cout << "Error: size A != size B in addition matrix." << std::endl;
        exit(1);
    }

    // Создание результирующей матрицы с теми же размерами
    std::vector<std::vector<T>> result(A.size(), std::vector<T>(A[0].size(), 0));

    // Поэлементное сложение
    for (size_t i = 0; i < A.size(); ++i) {
        for (size_t j = 0; j < A[0].size(); ++j) {
            result[i][j] = A[i][j] + B[i][j];
        }
    }

    return result;
}


/* Операция поэлементного вычитания матриц */
template <typename T>
std::vector<std::vector<T>> operator-(const std::vector<std::vector<T>>& A, const std::vector<std::vector<T>>& B){
    // Проверка на совпадение размеров матриц
    if (A.size() != B.size() || A[0].size() != B[0].size()) {
        std::cout << "Error: size A != size B in substraction matrix." << std::endl;
        exit(1);
    }

    // Создание результирующей матрицы с теми же размерами
    std::vector<std::vector<T>> result(A.size(), std::vector<T>(A[0].size(), 0));

    // Поэлементное сложение
    for (size_t i = 0; i < A.size(); ++i) {
        for (size_t j = 0; j < A[0].size(); ++j) {
            result[i][j] = A[i][j] - B[i][j];
        }
    }
    return result;
}


/* Операция умножения матрицы на вектор */
template <typename T>
std::vector<T> operator*(const std::vector<std::vector<T>>& matrix, const std::vector<T>& vec) {
    // Проверка на возможность умножения
    if (matrix[0].size() != vec.size()) {
        std::cout << "Error: size A != size b in multiply Matrix By Vector." << std::endl;
        exit(1);
    }
    // Создание результирующего вектора
    std::vector<T> result(matrix.size(), 0);

    // Умножение матрицы на вектор
    for (int i = 0; i < matrix.size(); ++i) {
        for (size_t j = 0; j < matrix[0].size(); ++j) {
            result[i] += matrix[i][j] * vec[j];
        }
    }
    return result;
}


/* Матричное умножение */
template <typename T>
std::vector<std::vector<T>> operator*(const std::vector<std::vector<T>>& A, const std::vector<std::vector<T>>& B){
    int m = A.size();    // Количество строк в матрице A
    int n = A[0].size(); // Количество столбцов в матрице A
    int p = B[0].size(); // Количество столбцов в матрице B

    if (n != B.size()) {
        std::cout << "Error: impossible multiply matrix" << std::endl;
        exit(1);
    }

    std::vector<std::vector<T>> result(m, std::vector<T>(p, 0.0));

    for (int i = 0; i < m; i++) {
        for (int j = 0; j < p; j++) {
            for (int k = 0; k < n; k++) {
                result[i][j] += A[i][k] * B[k][j];
            }
        }
    }
    return result;
}


// Определение оператора отрицания для матрицы
template <typename T>
std::vector<std::vector<T>> operator-(const std::vector<std::vector<T>>& matrix) {
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
std::vector<std::vector<T>> Matrix_round(const std::vector<std::vector<T>>& A, const double& eps){
    std::vector<std::vector<T>> roundA = A;
    int size = A.size();

    for (int i = 0; i < size; ++i) {
        for (int j = 0; j < size; ++j) {
            roundA[i][j] = (round(A[i][j]) >= 0)? round(abs(A[i][j]) * (1 / eps)) / (1 / eps): -1 * round(abs(A[i][j]) * (1 / eps)) / (1 / eps);
        }
    }
    return roundA;
}


/* Функция для вычисления нормы матрицы */
template <typename T>
T norm(const std::vector<std::vector<T>>& matrix, const int& p) {
    // Проверка на пустую матрицу
    if (matrix.empty() || matrix[0].empty()) {
        std::cout << "Error: Empty matrix in norm()\n" << std::endl;
        exit(1);
    }

    int rows = matrix.size();
    int cols = matrix[0].size();

    T result = 0.0;

    // Вычисление нормы матрицы
    if (p == 0) {
        // Норма матрицы Чебышева (максимальное значение по модулю в строке)
        for (int i = 0; i < rows; ++i) {
            T rowSum = 0.0;
            for (int j = 0; j < cols; ++j) {
                rowSum += abs(matrix[i][j]);
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
                colSum += pow(abs(matrix[i][j]), p);
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

    std::vector<std::vector<T>> A_rotate(A.size(), std::vector<T>(A.size(), 0));

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

    std::vector<std::vector<T>> A_rotate(A.size(), std::vector<T>(A.size(), 0));

    for (int i = 0; i < A.size(); ++i) {
        for (int j = 0; j < A.size(); ++j) {
            A_rotate[j][A.size() - 1 - i] = A[i][j];
        }
    }

    return A_rotate;
}


// Функция для создания единичной матрицы размера n x n
template <typename T>
std::vector<std::vector<T>> create_identity_matrix(const int& n) {
    std::vector<std::vector<T>> identity(n, std::vector<T>(n, 0));
    for (int i = 0; i < n; i++) {
        identity[i][i] = 1;
    }
    return identity;
}


// Функция для обратной матрицы с проверкой на вырожденность
template <typename T>
std::vector<std::vector<T>> inverseMatrix(const std::vector<std::vector<T>>& A, const T& eps) {
    std::vector<std::vector<T>> E = create_identity_matrix<T>(A.size());
    std::vector<std::vector<T>> E_rotate = RotateLeft(E);
    std::vector<T> e(A.size());
    std::vector<std::vector<T>> X(A.size(), std::vector<T>(A.size(), 0));


    for (int i = 0; i < A.size(); i++){
        e = E_rotate[i];
        X[i] = method_Gaussa(A, e, eps);

    }
    std::vector<std::vector<T>> A_inv = RotateLeft(X);
    return A_inv;
}


// Функция для обратной матрицы с проверкой на вырожденность с максимальной точностью
template <typename T>
std::vector<std::vector<T>> inverseMatrix(const std::vector<std::vector<T>>& A){
    T eps = std::numeric_limits<T>::epsilon();
    return inverseMatrix(A, eps);
}


/* Функция для вычисления числа обусловленности матрицы c нормой 1*/
template <typename T>
T cond_1(const std::vector<std::vector<T>>& matrix){
    T n_1 = norm_1(matrix);
    if (n_1 == 0) {
        printf("Error: Det(A) = 0  =>  cond_1(A) = oo");
        return std::numeric_limits<T>::infinity();
    }
    std::vector<std::vector<T>> inverse_matrix = inverseMatrix(matrix);
    T n_2 = norm_1(inverse_matrix);
    T cond = n_1 * n_2;
    return cond;
}


/* Функция для вычисления числа обусловленности матрицы c нормой 2*/
template <typename T>
T cond_2(const std::vector<std::vector<T>>& matrix){
    T n_1 = norm_2(matrix);
    if (n_1 == 0) {
        std::cout << "Error: Det(A) = 0  =>  cond_2(A) = oo" << std::endl;
        return std::numeric_limits<T>::infinity();
    }
    std::vector<std::vector<T>> inverse_matrix = inverseMatrix(matrix);
    T n_2 = norm_2(inverse_matrix);
    T cond = n_1 * n_2;
    return cond;
}


/* Функция для вычисления числа обусловленности матрицы с нормой oo*/
template <typename T>
T cond_oo(const std::vector<std::vector<T>>& matrix){
    T n_1 = norm_oo(matrix);
    if (n_1 == 0) {
        printf("Error: Det(A) = 0  =>  cond_oo(A) = oo");
        return std::numeric_limits<T>::infinity();
    }
    std::vector<std::vector<T>> inverse_matrix = inverseMatrix(matrix);
    T n_2 = norm_oo(inverse_matrix);
    T cond = n_1 * n_2;
    return cond;
}

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

/* Функция, вычисляющая определитель матрицы 4х4 */
template <typename T>
double det(const std::vector<std::vector<T>>& matrix) {
    return
            matrix[0][0] * (
                    matrix[1][1] * (matrix[2][2] * matrix[3][3] - matrix[2][3] * matrix[3][2]) -
                    matrix[1][2] * (matrix[2][1] * matrix[3][3] - matrix[2][3] * matrix[3][1]) +
                    matrix[1][3] * (matrix[2][1] * matrix[3][2] - matrix[2][2] * matrix[3][1])
            ) -
            matrix[0][1] * (
                    matrix[1][0] * (matrix[2][2] * matrix[3][3] - matrix[2][3] * matrix[3][2]) -
                    matrix[1][2] * (matrix[2][0] * matrix[3][3] - matrix[2][3] * matrix[3][0]) +
                    matrix[1][3] * (matrix[2][0] * matrix[3][2] - matrix[2][2] * matrix[3][0])
            ) +
            matrix[0][2] * (
                    matrix[1][0] * (matrix[2][1] * matrix[3][3] - matrix[2][3] * matrix[3][1]) -
                    matrix[1][1] * (matrix[2][0] * matrix[3][3] - matrix[2][3] * matrix[3][0]) +
                    matrix[1][3] * (matrix[2][0] * matrix[3][1] - matrix[2][1] * matrix[3][0])
            ) -
            matrix[0][3] * (
                    matrix[1][0] * (matrix[2][1] * matrix[3][2] - matrix[2][2] * matrix[3][1]) -
                    matrix[1][1] * (matrix[2][0] * matrix[3][2] - matrix[2][2] * matrix[3][0]) +
                    matrix[1][2] * (matrix[2][0] * matrix[3][1] - matrix[2][1] * matrix[3][0]));
}