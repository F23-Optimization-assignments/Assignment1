#include <gtest/gtest.h>
#include <sstream>
#include "src/simplex.h"
#include "src/fraction.h"


TEST(SimplexOptimumTest, Example1) {
    // Define the input data for the Simplex function
    std::vector<Fraction> coeffs = {-2, 3, -6, -1, 0, 0};
    Matrix<Fraction> m(3, 6);
    std::istringstream input_matrix("2 1 -2 1 0 0\n1 2 4 0 1 0\n1 -1 2 0 0 1");
    input_matrix >> m;
    std::vector<Fraction> b = {24, 22, 10};
    std::vector<size_t> basic_indices = {3, 4, 5};

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

TEST(SimplexOptimumTest, Example2) {
    std::vector<Fraction> coeffs = {-6, -8, -5, -9, 0, 0};
    Matrix<Fraction> m(2, 6);
    std::istringstream input_matrix("2 1 1 3 1 0\n1 3 1 2 0 1");
    input_matrix >> m;
    std::vector<Fraction> b = {5, 3};
    std::vector<size_t> basic_indices = {4, 5};

    Simplex<Fraction> simplex(coeffs, m, b, basic_indices);
    Solution<Fraction> solution = simplex.find_solution();

    int expected_optimum = -17;

    EXPECT_EQ(solution.is_final(), true);
    EXPECT_EQ(solution.is_unbounded(), false);
    EXPECT_EQ(solution.get_optimum(), expected_optimum);
}

TEST(SimplexOptimumTest, Example3) {
    std::vector<Fraction> coeffs = {-2, -1, 0, 0, 0, 0};
    Matrix<Fraction> m(4, 6);
    std::istringstream input_matrix("2 1 1 0 0 0\n2 3 0 1 0 0\n4 1 0 0 1 0\n1 5 0 0 0 1");
    input_matrix >> m;
    std::vector<Fraction> b = {4, 3, 5, 1};
    std::vector<size_t> basic_indices = {2, 3, 4, 5};

    Simplex<Fraction> simplex(coeffs, m, b, basic_indices);
    Solution<Fraction> solution = simplex.find_solution();

    int expected_optimum = -2;

    EXPECT_EQ(solution.is_final(), true);
    EXPECT_EQ(solution.is_unbounded(), false);
    EXPECT_EQ(solution.get_optimum(), expected_optimum);
}

TEST(SimplexOptimumTest, Example4) {
    std::vector<Fraction> coeffs = {-2, -3, -4, 0, 0, 0};
    Matrix<Fraction> m(3, 6);
    std::istringstream input_matrix("0 -2 -3 -1 0 0\n1 1 2 0 1 0\n1 2 3 0 0 1");
    input_matrix >> m;
    std::vector<Fraction> b = {-5, 4, 7};
    std::vector<size_t> basic_indices = {3, 4, 5};

    Simplex<Fraction> simplex(coeffs, m, b, basic_indices);
    Solution<Fraction> solution = simplex.find_solution();

    Fraction expected_optimum = Fraction(-21, 2);

    EXPECT_EQ(solution.is_final(), true);
    EXPECT_EQ(solution.is_unbounded(), false);
    EXPECT_EQ(solution.get_optimum(), expected_optimum);
}

TEST(SimplexOptimumTest, Example5) {
    std::vector<Fraction> coeffs = {-2, -1, 0, 0};
    Matrix<Fraction> m(2, 4);
    std::istringstream input_matrix("1 -1 1 0\n2 0 0 1");
    input_matrix >> m;
    std::vector<Fraction> b = {10, 40};
    std::vector<size_t> basic_indices = {2, 3};

    Simplex<Fraction> simplex(coeffs, m, b, basic_indices);
    Solution<Fraction> solution = simplex.find_solution();

    EXPECT_EQ(solution.is_final(), true);
    EXPECT_EQ(solution.is_unbounded(), true);
}

TEST(SimplexOptimumTest, ExampleFromReport1) {
    std::vector<Fraction> coeffs = {-4, 12, 0, 0, 0};
    Matrix<Fraction> m(2, 5);
    std::istringstream input_matrix("2 2 4 1 0\n2 -2 3 0 1");
    input_matrix >> m;
    std::vector<Fraction> b = {2, 1};
    std::vector<size_t> basic_indices = {3, 4};

    Simplex<Fraction> simplex(coeffs, m, b, basic_indices);
    Solution<Fraction> solution = simplex.find_solution();

    Fraction expected_optimum = -2;

    EXPECT_EQ(solution.is_final(), true);
    EXPECT_EQ(solution.is_unbounded(), false);
    EXPECT_EQ(solution.get_optimum(), expected_optimum);
}

TEST(SimplexOptimumTest, ExampleFromReport2) {
    std::vector<Fraction> coeffs = {-4, -5, 0, 0, 0, 0};
    Matrix<Fraction> m(2, 6);
    std::istringstream input_matrix("2 3 -1 0 1 0\n8 6 0 -1 0 1");
    input_matrix >> m;
    std::vector<Fraction> b = {8, 9};
    std::vector<size_t> basic_indices = {4, 5};

    Simplex<Fraction> simplex(coeffs, m, b, basic_indices);
    Solution<Fraction> solution = simplex.find_solution();


    EXPECT_EQ(solution.is_final(), true);
    EXPECT_EQ(solution.is_unbounded(), true);
}

// Run all tests
int main(int argc, char **argv) {
    testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
