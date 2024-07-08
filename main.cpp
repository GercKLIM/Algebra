#include <iostream>
#include "Lib/algebra.h"

int main() {
    std::cout << "E = " << std::endl;
    std::vector<std::vector<double>> E = create_identity_matrix<double>(5);
    print(E);
    //std::cout << E << std::endl;
    return 0;
}
