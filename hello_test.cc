#include <gtest/gtest.h>
#include <sstream>
#include "src/simplex.h"
#include "src/fraction.h"


TEST(SimplexOptimumTest, Example1) {
    // Define the input data for the Simplex function
    std::vector<Fraction> coeffs = { -2, 3, -6, -1, 0, 0 };
    Matrix<Fraction> m(3, 6);
    std::istringstream input_matrix("2 1 -2 1 0 0\n1 2 4 0 1 0\n1 -1 2 0 0 1");
    input_matrix >> m;
    std::vector<Fraction> b = { 24, 22, 10 };
    std::vector<size_t> basic_indices = { 3, 4, 5 };

    // Create a Simplex object and find the solution
    Simplex<Fraction> simplex(coeffs, m, b, basic_indices);
    Solution<Fraction> solution = simplex.find_solution();

    // Define the expected optimum value
    int expected_optimum = -64;

    // Check if the actual optimum matches the expected optimum
    EXPECT_EQ(solution.is_final(), true);
    EXPECT_EQ(solution.is_unbounded(), false);
    EXPECT_EQ(solution.get_optimum(), expected_optimum);
}