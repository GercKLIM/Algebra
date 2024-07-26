#include <iostream>
#include "Lib/algebra.h"

int main() {

    //std::vector<Vector<double>> E = {{1, 0, 0}, {0, 1, 0}, {0, 0, 1}};
    std::vector<std::vector<double>> E = {{1, 0, 0}, {0, 1, 0}, {0, 0, 1}};
    Matrix<double> m1(E);
    std::cout << m1 - m1 << std::endl;


    std::cout << "Complete!" << std::endl;
    return 0;
}
