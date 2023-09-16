#include <iostream>
#include <vector>

#include "fraction.h"
#include "matrix.h"
#include "vector_ops.h"

class SimplexException : public std::exception {
private:
    const char* msg;
public:
    explicit SimplexException(const char* msg) : std::exception(), msg(msg) { }
    [[nodiscard]] const char * what() const noexcept override {
        return msg;
    }
};

template<typename T>
class Simplex {
private:
    std::vector<T> func;
    Matrix<T> A;
    std::vector<T> b;
    std::vector<size_t> basic_indices;
    std::vector<T> basic_coeffs;

    void validate() const {
        // TODO: implement
    }

    [[nodiscard]] std::vector<T> find_delta() const {
        std::vector<T> delta_row;
        for (size_t idx = 0; idx < func.size(); ++idx) {
            delta_row.push_back(dot_product(A.get_column(idx), basic_coeffs) - func[idx]);
        }
        return delta_row;
    }

    [[nodiscard]] size_t find_pivot_col(const std::vector<T>& v) const {
        size_t max_idx = 0;
        for (size_t idx = 1; idx < v.size(); ++idx) {
            if (v[idx] > v[max_idx]) {
                max_idx = idx;
            }
        }
        return max_idx;
    }

    [[nodiscard]] size_t find_pivot_row(size_t col) const {
        bool found = false;
        size_t row_idx = 0;
        T best_ratio;
        for (size_t row = 0; row < A.get_rows(); ++row) {
            auto ratio = b[row] / A[row][col];
            if (ratio <= 0) { // TODO: check!
                continue;
            }
            if (!found) {
                found = true;
                row_idx = row;
                best_ratio = ratio;
            } else if (ratio < best_ratio) {
                row_idx = row;
            }
        }
        if (!found) {
            throw SimplexException("did not found positive ratio pivot");
        }
        return row_idx;
    }

    void update_pivot_row(size_t pivot_row, size_t pivot_col)  {
        auto pivot = A[pivot_row][pivot_col];
        b[pivot_row] /= pivot;
        for (size_t idx = 0; idx < A.get_columns(); ++idx) {
            A[pivot_row][idx] /= pivot;
        }
        basic_indices[pivot_row] = pivot_col;
        basic_coeffs[pivot_row] = func[pivot_col];
    }

    void update_constraints_values(size_t pivot_row, size_t pivot_col) {
        for (size_t row = 0; row < A.get_rows(); ++row) {
            if (pivot_row == row) {
                continue;
            }
            b[row] -= A[row][pivot_col] * b[pivot_row];
        }
    }

    void update_other_rows(size_t pivot_row, size_t pivot_col) {
        for (size_t row = 0; row < A.get_rows(); ++row) {
            if (row == pivot_row) {
                continue;
            }
            for (size_t col = 0; col < A.get_columns(); ++col) {
                if (col == pivot_col) {
                    continue;
                }
                A[row][col] -= A[row][pivot_col] * A[pivot_row][col];
            }
        }
    }

    void update_pivot_column(size_t pivot_row, size_t pivot_col) {
        for (size_t row = 0; row < A.get_rows(); ++row) {
            if (row == pivot_row) {
                A[row][pivot_col] = 1;
            } else {
                A[row][pivot_col] = 0;
            }
        }
    }

    bool iterate()  {
        auto delta = find_delta();
        auto pivot_col = find_pivot_col(delta);
        if (delta[pivot_col] <= 0) {
            return false;
        }
        auto pivot_row = find_pivot_row(pivot_col);
        update_pivot_row(pivot_row, pivot_col);
        update_constraints_values(pivot_row, pivot_col);
        update_other_rows(pivot_row, pivot_col);
        update_pivot_column(pivot_row, pivot_col);
        return true;
    }

    [[nodiscard]] T function_value() const {
        return dot_product(b, basic_coeffs);
    }

    void print_solution() const {
        std::cout << "Optimum is " << function_value() << '\n';
        for (size_t idx = 0; idx < basic_indices.size(); ++idx) {
            std::cout << "Variable x_" << basic_indices[idx] << " is " << b[idx] << '\n';
        }
        std::cout << "Others are zeros.\n";
    }
public:
    Simplex(const std::vector<T>& coefficients,
            const Matrix<T>& A,
            const std::vector<T>& b,
            const std::vector<size_t>& basic_indices) : func(
            coefficients), A(A), b(b), basic_indices(
            basic_indices), basic_coeffs() {
        validate();
        for (const size_t& idx : basic_indices) {
            basic_coeffs.push_back(func[idx]);
        }
    }

    void find_solution() {
        while (iterate());
        print_solution();
    }
};

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
/*
6 2
-6 -8 -5 -9 0 0

2 1 1 3 1 0
1 3 1 2 0 1

5 3
4 5
 */

//Ans = -2
/*
6 4
-2 -1 0 0 0 0

2 1 1 0 0 0
2 3 0 1 0 0
4 1 0 0 1 0
1 5 0 0 0 1

4 3 5 1
2 3 4 5
 */