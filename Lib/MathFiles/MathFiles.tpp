/*
 *  РЕАЛИЗАЦИЯ ФУНКЦИЙ ИМПОРТА И ЭКСПОРТА
 *  ВЕКТОРОВ И МАТРИЦ ИЗ *.txt ФАЙЛОВ
 *
 * */


#include "../algebra.h"


/* ### ФУНКЦИИ ИМПОРТА ### */


/* Функция импорта чисел из файла */
std::vector<double> ImportData(const std::string& filename) {
    std::vector<double> numbers;

    // Открытие файла для чтения
    std::ifstream file(filename);
    if (!file.is_open()) {
        std::cerr << "Failed to open file: " << filename << std::endl;
        return numbers;
    }

    std::string line;
    while (std::getline(file, line)) {
        // Преобразование строки в число типа double и добавление его в вектор
        double num = std::atof(line.c_str());
        numbers.push_back(num);
    }

    // Закрытие файла
    file.close();

    return numbers;
}


/* Функция импорта матрицы из текстового файла */
template <typename T>
std::vector<std::vector<T>> importSLAU(const std::string& filename) {
    std::vector<std::vector<T>> matrix;
    std::vector<T> vec;
    std::ifstream file(filename);

    if (!file.is_open()) {
        std::cout << "Error: not open file \n" << std::endl;
        exit(1);
    }

    int size;
    file >> size;

    matrix.resize(size, std::vector<T>(size+1));

    for (int i = 0; i < size; ++i) {
        for (int j = 0; j < size; ++j) {
            T value;
            if (file >> value) {
                matrix[i][j] = value;
            }
        }
    }

    file.close();
    return matrix;
};


/* Переопределение потока вывода для vector  */
template<typename T>
std::ostream& operator<<(std::ostream& os, const std::vector<T>& vec) {
    //os << "[";
    for (int i = 0; i < vec.size(); ++i) {
        os << std::setprecision(16) << vec[i];
        if (i != vec.size() - 1) {
            os << " ";
            //os << ", ";
        }
    }
    //os << "]";
    os << " ";
    return os;
}
