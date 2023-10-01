#include <iostream>
#include <vector>

#include "src/fraction.h"
#include "src/matrix.h"
#include "src/simplex.h"

int main() {
    // Input number of ALL variables and constraints
    size_t vars, cons;
    std::cin >> vars >> cons;

    // Input function's coefficients
    std::vector<Fraction> coeffs(vars);
    for (auto& i : coeffs) {
        std::cin >> i;
    }

    // Input matrix A of constraints
    Matrix<Fraction> m(cons, vars);
    std::cin >> m;

    // Input vector of values of constraints
    std::vector<Fraction> b(cons);
    for (auto& i : b) {
        std::cin >> i;
    }

    // Input indices of initial basic vars in corresponding constraints
    std::vector<size_t> basic(cons);
    for (auto& i : basic) {
        std::cin >> i;
    }

    // Compose simplex, find and print solution
    Simplex<Fraction> simplex(coeffs, m, b, basic);
    std::cout << simplex.find_solution();
}

// Ex. Test1: Incorrect number of basic variables
/**
6 3

-2 3 -6 -1 0 0

2 1 -2 0 0 0
1 2 4 0 1 0
1 -1 2 0 0 1

24 22 10

3 4 5
*/

// Ex. Test2: Incorrect index of the base variable
/**
6 3

-2 3 -6 -1 0 0

2 1 -2 2 0 0
1 0 4 0 1 0
1 0 2 0 0 1

24 22 10

3 4 5
*/

// Ex. 1 : Zlata's Lab Sample
// Answer = -64
/**
6 3

-2 3 -6 -1 0 0

2 1 -2 1 0 0
1 2 4 0 1 0
1 -1 2 0 0 1

24 22 10

3 4 5
*/

// Ex.2 : Assignment task
// Answer = -17
/**
6 2
-6 -8 -5 -9 0 0

2 1 1 3 1 0
1 3 1 2 0 1

5 3
4 5
*/

// Ex. 3 : Assignment task
// Answer = -2
/**
6 4
-2 -1 0 0 0 0

2 1 1 0 0 0
2 3 0 1 0 0
4 1 0 0 1 0
1 5 0 0 0 1

4 3 5 1
2 3 4 5
*/

// Ex. 4 : Assignment task
// Answer = -21/2
/**
6 3
-2 -3 -4 0 0 0

0 -2 -3 -1 0 0
1 1 2 0 1 0
1 2 3 0 0 1

-5 4 7
3 4 5
*/

// Ex. 5 : Tutorial example
// Answer = Unbounded
/**
4 2
-2 -1 0 0

1 -1 1 0
2 0 0 1

10 40

2 3
*/