#include <iostream>
#include <vector>

#include "src/fraction.h"
#include "src/matrix.h"
#include "src/simplex.h"

int main() {
    size_t vars, cons;
    std::cin >> vars >> cons;

    std::vector<Fraction> coeffs(vars);
    for (auto& i : coeffs) {
        std::cin >> i;
    }

    Matrix<Fraction> m(cons, vars);
    std::cin >> m;

    std::vector<Fraction> b(cons);
    for (auto& i : b) {
        std::cin >> i;
    }

    std::vector<size_t> basic(cons);
    for (auto& i : basic) {
        std::cin >> i;
    }
    Simplex<Fraction> simplex(coeffs, m, b, basic);
    simplex.find_solution();
}

// Zlata's Lab Sample
/**
6 3

-2 3 -6 -1 0 0

2 1 -2 1 0 0
1 2 4 0 1 0
1 -1 2 0 0 1

24 22 10

3 4 5
 */

//Other example, answer = -17
/**
6 2
-6 -8 -5 -9 0 0

2 1 1 3 1 0
1 3 1 2 0 1

5 3
4 5
 */

//Ans = -2
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

/**
ans = -21/2
6 3
-2 -3 -4 0 0 0

0 -2 -3 -1 0 0
1 1 2 0 1 0
1 2 3 0 0 1

-5 4 7
3 4 5
*/