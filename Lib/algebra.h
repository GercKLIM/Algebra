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


/* ВЫЗОВ РЕАЛИЗАЦИЙ */

/* Функции для Импорта/Экспорта матриц и векторов их тектового файла */
#include "MathFiles/MathFiles.h"

/* Функции алгебры векторов */
#include "vectors/vectors.h"

/* Функции алгебры матриц */
#include "matrixes/matrixes.h"

/* Функции алгебры СЛАУ */
#include "SLAE/SLAE.h"

/* Некоторые другие функции математики */
#include "maths/maths.h"
