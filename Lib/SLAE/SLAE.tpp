/*
 *  РЕАЛИЗАЦИЯ АЛГЕБРЫ CЛАУ
 *
 * */

#include "../Algebra.h"


/* ### ФУНКЦИИ КОНВЕРТАЦИЙ в СЛАУ ### */


/* Функция для получения матрицы из СЛАУ */
template <typename T>
Matrix<T> SLAE<T>::to_Matrix(){
    std::vector<std::vector<T>> matrix;
    matrix.resize(A_data.size(), std::vector<T>(A_data.size()));

    for (int i = 0; i < A_data.size(); i++) {
        for (int j = 0; j < A_data.size(); j++) {
            matrix[i][j] = A_data[i][j];
        }
    }
    return matrix;
}


/* Функция для получения вектора из СЛАУ */
template <typename T>
Vector<T> SLAE<T>::to_Vector(){
    int s = A_data.size();
    std::vector<T> vec(s);

    for (int i = 0; i < A_data.size(); i++) {
        vec[i] = A_data[i][s];
    }
    return vec;
}



/* ### ОПЕРАЦИИ СО СЛАУ ### */



/* Операция сравнения */
template<typename Y>
bool operator==(const SLAE<Y>& slae1, const SLAE<Y>& slae2){
    if ((slae1.to_Vector() == slae2.to_Vector()) and (slae1.to_Matrix() == slae2.to_Matrix())){
        return true;
    }
    return false;
}


template<typename T>
Vector<T> SLAE<T>::residual(const Vector<T> sol){
    int n = A_data.size();
    Vector<T> residual(n, 0);


    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            residual[i] += A_data[i][j] * sol[j];
        }
        residual[i] = b_data[i] - residual[i];
    }
    return residual;
}


//template<typename T>
//T SLAE<T>::residual(const Vector<T> sol, const int& norm_type){
//    int n = A_data.size();
//    Vector<T> residual(n, 0);
//
//
//    for (int i = 0; i < n; i++) {
//        for (int j = 0; j < n; j++) {
//            residual[i] += A_data[i][j] * sol[j];
//        }
//        residual[i] = b_data[i] - residual[i];
//    }
//    return residual.norm(norm_type);
//}


/* Функция для проверки решения СЛАУ*/
template<typename T>
bool SLAE<T>::is_sol(const Vector<T> sol){

    int n = A_data.size();
    Vector<T> residual(n, 0);


    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            residual[i] += A_data[i][j] * sol[j];
        }
        residual[i] = b_data[i] - residual[i];
    }
    return (residual.norm(2) == 0);
}


/* Функция QR-разложение */
template<typename T>
std::pair<Matrix<T>, Matrix<T>> SLAE<T>::QR(){
    int rows = A_data.rows();
    int cols = A_data.cols();

    Matrix<T> Q(rows, cols); // Ортогональная матрица
    Matrix<T> R(cols, cols); // Верхнетреугольная матрица

    for (int k = 0; k < cols; ++k) {
        // 1. Инициализируем q_k как k-й столбец матрицы A
        Vector<T> q_k(rows);
        for (int i = 0; i < rows; ++i) {
            q_k[i] = A_data(i, k);
        }

        // 2. Проекция q_k на все предыдущие столбцы q_j и вычитание проекций
        for (int j = 0; j < k; ++j) {
            T r_jk = 0;
            for (int i = 0; i < rows; ++i) {
                r_jk += Q(i, j) * A_data(i, k);
            }

            R(j, k) = r_jk;
            for (int i = 0; i < rows; ++i) {
                q_k[i] -= r_jk * Q(i, j);
            }
        }

        // 3. Норма вектора q_k и нормализация
        T norm_q_k = q_k.norm();
        R(k, k) = norm_q_k;

        for (int i = 0; i < rows; ++i) {
            Q(i, k) = q_k[i] / norm_q_k;
        }
    }

    return std::make_pair(Q, R);
}


/* Функция LU-разложение */
template<typename T>
std::pair<Matrix<T>, Matrix<T>> SLAE<T>::LU() {
    int n = A_data.rows();
    Matrix<T> L(n, n); // Нижнетреугольная матрица
    Matrix<T> U(n, n); // Верхнетреугольная матрица

    // Инициализация матриц L и U нулями
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            L(i, j) = 0;
            U(i, j) = 0;
        }
    }

    // LU-разложение
    for (int i = 0; i < n; ++i) {
        // Вычисление верхнетреугольной матрицы U
        for (int k = i; k < n; ++k) {
            T sum = 0;
            for (int j = 0; j < i; ++j) {
                sum += L(i, j) * U(j, k);
            }
            U(i, k) = A_data(i, k) - sum;
        }

        // Вычисление нижнетреугольной матрицы L
        for (int k = i; k < n; ++k) {
            if (i == k) {
                L(i, i) = 1; // Диагональные элементы L равны 1
            } else {
                T sum = 0;
                for (int j = 0; j < i; ++j) {
                    sum += L(k, j) * U(j, i);
                }
                L(k, i) = (A_data(k, i) - sum) / U(i, i);
            }
        }
    }

    return std::make_pair(L, U);
}



/* ### ФУНКЦИИ РЕШЕНИЯ СЛАУ ### */



/* Вывод вычисленного решения */
template<typename T>
Vector<T> SLAE<T>::sol(){
    return sol_data;
}


/* Решение СЛАУ c автоматическим выбором метода */
template<typename T>
Vector<T> SLAE<T>::solve(T eps) {
    // Размер матрицы
    int n = A_data.rows(); // Предполагаем, что метод rows() возвращает количество строк
    if (n != A_data.cols() || b_data.size() != n) {
        throw std::invalid_argument("Matrix and vector dimensions are inconsistent.");
    }

    // Проверка на размер матрицы
    if (n > 5) {
        // Используем QR-разложение для больших матриц
        return solve_QR(eps);
    } else if (n == A_data.cols()) {
        // Проверяем, применим ли метод Крамера
        T det_A = A_data.det();
        if (std::abs(det_A) > eps) {
            // Используем метод Крамера, если детерминант не близок к нулю
            return solve_Cramer(eps);
        } else {
            // Если детерминант близок к нулю, используем метод Гаусса
            return solve_Gauss(eps);
        }
    } else {
        // Если матрица не квадратная, использовать метод Гаусса
        return solve_Gauss(eps);
    }
}



/* Функция прямого хода метода Гаусса */
template<typename T>
void SLAE<T>::forward_Gauss(const T& eps) {
    int n = A_data.rows();

    for (int i = 0; i < n; ++i) {
        // Поиск главного элемента в текущем столбце
        T maxElem = std::abs(A_data(i, i));
        int maxRow = i;
        for (int k = i + 1; k < n; ++k) {
            if (std::abs(A_data(k, i)) > maxElem) {
                maxElem = std::abs(A_data(k, i));
                maxRow = k;
            }
        }

        // Если главный элемент слишком мал, то система может быть вырождена
        if (maxElem < eps) {
            throw std::runtime_error("The system is singular or nearly singular.");
        }

        // Перестановка строк
        for (int k = i; k < n; ++k) {
            std::swap(A_data(i, k), A_data(maxRow, k));
        }
        std::swap(b_data[i], b_data[maxRow]);

        // Нормализация строки с ведущим элементом
        T leadElem = A_data(i, i);
        for (int k = i; k < n; ++k) {
            A_data(i, k) /= leadElem;
        }
        b_data[i] /= leadElem;

        // Обнуление всех элементов ниже ведущего
        for (int k = i + 1; k < n; ++k) {
            T factor = A_data(k, i);
            for (int j = i; j < n; ++j) {
                A_data(k, j) -= factor * A_data(i, j);
            }
            b_data[k] -= factor * b_data[i];
        }
    }
}


/* Функция обратного хода метода Гаусса */
template<typename T>
Vector<T> SLAE<T>::back_Gauss() {
    int n = A_data.rows();
    Vector<T> x(n);

    for (int i = n - 1; i >= 0; --i) {
        x[i] = b_data[i];
        for (int j = i + 1; j < n; ++j) {
            x[i] -= A_data(i, j) * x[j];
        }
    }

    return x;
}



/* Решение СЛАУ методом Гаусса */
template<typename T>
Vector<T> SLAE<T>::solve_Gauss(T eps) {
    int n = A_data.rows();
    Matrix<T> A = A_data;  // Копия матрицы коэффициентов
    Vector<T> b = b_data;  // Копия вектора правой части

    // Прямой ход
    for (int i = 0; i < n; ++i) {
        // Поиск главного элемента для выбора ведущей строки
        T maxElem = std::abs(A(i, i));
        int maxRow = i;
        for (int k = i + 1; k < n; ++k) {
            if (std::abs(A(k, i)) > maxElem) {
                maxElem = std::abs(A(k, i));
                maxRow = k;
            }
        }

        // Если главный элемент слишком мал, то система может быть вырождена
        if (maxElem < eps) {
            throw std::runtime_error("The system is singular or nearly singular.");
        }

        // Перестановка строк
        for (int k = i; k < n; ++k) {
            std::swap(A(i, k), A(maxRow, k));
        }
        std::swap(b[i], b[maxRow]);

        // Нормализация строки с ведущим элементом
        T leadElem = A(i, i);
        for (int k = i; k < n; ++k) {
            A(i, k) /= leadElem;
        }
        b[i] /= leadElem;

        // Обнуление всех элементов ниже ведущего
        for (int k = i + 1; k < n; ++k) {
            T factor = A(k, i);
            for (int j = i; j < n; ++j) {
                A(k, j) -= factor * A(i, j);
            }
            b[k] -= factor * b[i];
        }
    }

    // Обратный ход
    Vector<T> x(n);
    for (int i = n - 1; i >= 0; --i) {
        x[i] = b[i];
        for (int j = i + 1; j < n; ++j) {
            x[i] -= A(i, j) * x[j];
        }
    }

    return x;
}


/* Решение СЛАУ методом Гаусса (вариант 2)*/
template<typename T>
Vector<T> SLAE<T>::solve_Gauss_2(T eps) {
    Matrix<T> A = A_data;  // Копия матрицы коэффициентов
    Vector<T> b = b_data;  // Копия вектора правой части

    // Прямой ход
    forward_Gauss(A, b, eps);

    // Обратный ход
    return back_Gauss(A, b);
}


/* Решение СЛАУ методом QR-разложения */
template<typename T>
Vector<T> SLAE<T>::solve_QR(T eps) {
    // 1. Получение QR-разложения матрицы A
    auto [Q, R] = QR(); // Вызов функции QR-разложения

    int n = R.rows();
    Vector<T> b_tilde(n);

    // 2. Вычисление Q^T * b
    // Q^T * b = Q транспонированная умноженная на b
    for (int i = 0; i < n; ++i) {
        b_tilde[i] = 0;
        for (int j = 0; j < n; ++j) {
            b_tilde[i] += Q(j, i) * b_data[j];
        }
    }

    // 3. Решение системы R * x = b_tilde методом Гаусса
    // Обратное замещение
    Vector<T> x(n);

    for (int i = n - 1; i >= 0; --i) {
        x[i] = b_tilde[i];
        for (int j = i + 1; j < n; ++j) {
            x[i] -= R(i, j) * x[j];
        }
    }

    return x;
}



/* Решение СЛАУ методом Крамера */
template<typename T>
Vector<T> SLAE<T>::solve_Cramer(T eps) {
    int n = A_data.rows();

    // Проверка, что матрица квадратная
    if (n != A_data.cols() || b_data.size() != n) {
        throw std::invalid_argument("Matrix A and vector b must have compatible dimensions.");
    }

    // Вычисление детерминанта матрицы A
    T det_A = A_data.det();
    if (std::abs(det_A) < eps) {
        throw std::runtime_error("The determinant of the coefficient matrix is zero or near zero. The system is singular or nearly singular.");
    }

    // Вектор решений
    Vector<T> x(n);

    // Решение для каждого неизвестного x_i
    for (int i = 0; i < n; ++i) {
        // Создание матрицы A_i (замена i-го столбца матрицы A на вектор b)
        Matrix<T> A_i = A_data;
        for (int j = 0; j < n; ++j) {
            A_i(j, i) = b_data[j];
        }

        // Вычисление детерминанта матрицы A_i
        T det_A_i = A_i.det();

        // Решение x_i
        x[i] = det_A_i / det_A;
    }

    return x;
}

