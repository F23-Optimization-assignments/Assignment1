#include <iostream>
#include <vector>

#include "src/fraction.h"
#include "src/matrix.h"
#include "src/simplex.h"

int main() {
    // Input type of the problem - min (0) OR max(1)
    int objNum;
    std::cin >> objNum;
    auto obj = static_cast<Objective>(objNum);

    // Input number of ALL variables and constraints
    size_t vars, cons;
    std::cin >> vars >> cons;

    // Input function's coefficients
    std::vector<Fraction> coeffs(vars);
    for (auto &i: coeffs) {
        std::cin >> i;
    }

    // Input matrix A of constraints
    Matrix<Fraction> m(cons, vars);
    std::cin >> m;

    // Input vector of values of constraints
    std::vector<Fraction> b(cons);
    for (auto &i: b) {
        std::cin >> i;
    }

    // Input indices of initial basic vars in corresponding constraints
    std::vector<size_t> basic(cons);
    for (auto &i: basic) {
        std::cin >> i;
    }

    // Compose simplex, find and print solution
    Simplex<Fraction> simplex(obj, coeffs, m, b, basic);
    std::cout << simplex.find_solution();
}