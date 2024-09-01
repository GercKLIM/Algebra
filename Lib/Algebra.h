//
// Объявление функций и переопределений для std::vector
// для возможности абстракции в математические векторы, матрицы и другой доп. функционал
//


#pragma once

#include <iostream>
#include <fstream>
#include<istream>
#include <sstream>
#include <vector>
#include <string>
#include <cmath>
#include <limits>
#include <iomanip>
#include <initializer_list>
#include <algorithm>


/* ВЫЗОВ РЕАЛИЗАЦИЙ */

/* Функции для Импорта/Экспорта матриц и векторов их тектового файла */
#include "MathFiles/MathFiles.h"

/* Функции алгебры векторов */
#include "Vectors/Vector.h"

/* Функции алгебры матриц */
#include "Matrix/Matrix.h"

/* Функции алгебры СЛАУ */
#include "SLAE/SLAE.h"

/* Функции алгебры уравнений (простых) */
#include "Equation/Equation.h"

/* Некоторые другие функции математики */
#include "Math/Math.h"
