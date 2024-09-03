/*
 *  БИБЛИОТЕКА "АЛГЕБРА"
 *
 */

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
#include <functional>
#include <stdexcept>


/* ВЫЗОВ РЕАЛИЗАЦИЙ */

/* Функции алгебры векторов */
#include "Vector/Vector.h"

/* Функции алгебры матриц */
#include "Matrix/Matrix.h"

/* Функции алгебры СЛАУ */
#include "SLAE/SLAE.h"

/* Функции алгебры уравнений (простых) */
#include "Equation/Equation.h"

/* Функции для Импорта/Экспорта матриц и векторов их тектового файла */
#include "MathFiles/MathFiles.h"

/* Некоторые другие функции математики */
#include "Math/Math.h"
