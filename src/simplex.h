#ifndef ASSIGNMENT1_SIMPLEX_H
#define ASSIGNMENT1_SIMPLEX_H

#include <optional>
#include "matrix.h"
#include "vector_ops.h"

template<typename T>
class Solution {
private:
    bool final;
    bool unbounded;
    T optimum;
    std::vector<T> xs;

    template<class U>
    friend std::ostream& operator<<(std::ostream&, const Solution<U>&);
public:
    Solution(bool isFinal,
             bool isUnbounded,
             T optimum,
             const std::vector<T>& xs) :
            final(isFinal), unbounded(isUnbounded), optimum(optimum), xs(xs) {}
    explicit Solution(const bool unbounded) : unbounded(unbounded), final(true) {}
    Solution() : final(false), unbounded(false) { };

    [[nodiscard]] bool is_final() const {
        return final;
    }

    [[nodiscard]] bool is_unbounded() const {
        return unbounded;
    }

    [[ nodiscard ]] std::vector<T> get_vector() const {
        return xs;
    }

    [[ nodiscard ]] T get_optimum() const {
        return optimum;
    }
};

template<class U>
std::ostream& operator<<(std::ostream& stream, const Solution<U>& solution) {
    if (solution.unbounded) {
        return (stream << "Solution is unbounded!\n");
    }
    stream << (solution.final ? "[ FINAL ]\n" : "[ NOT FINAL ]\n");
    stream << "Optimum is " << solution.optimum << '\n';
    for (size_t idx = 0; idx < solution.xs.size(); ++idx) {
        stream << "X_" << idx << " = " << solution.xs[idx] << '\n';
    }
    return stream;
}

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

    [[nodiscard]] std::optional<size_t> find_pivot_col(const std::vector<T>& v) const {
        std::optional<size_t> max_idx;
        for (size_t idx = 0; idx < v.size(); ++idx) {
            if (v[idx] > 0 && (!max_idx || v[max_idx.value()] < v[idx])) {
                max_idx = idx;
            }
        }
        return max_idx;
    }

    [[nodiscard]] std::optional<size_t> find_pivot_row(size_t col) const {
        std::optional<size_t> pivot_row;
        T best_ratio;
        for (size_t row = 0; row < A.get_rows(); ++row) {
            auto ratio = b[row] / A[row][col];
            if (ratio <= 0) {
                continue;
            }
            if (!pivot_row || ratio < best_ratio) {
                pivot_row = row;
                best_ratio = ratio;
            }
        }
        return pivot_row;
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

    [[nodiscard]] T function_value() const {
        return dot_product(b, basic_coeffs);
    }

    [[ nodiscard ]] std::vector<T> vector_value() const {
        std::vector<T> vars(func.size());
        for (size_t idx = 0; idx < basic_indices.size(); ++idx) {
            vars[basic_indices[idx]] = b[idx];
        }
        return vars;
    }

    Solution<T> iterate()  {
        // TODO: may be additional checks?
        auto delta = find_delta();
        auto col_opt = find_pivot_col(delta);
        if (!col_opt) {
            return Solution<T>(true, false, function_value(), vector_value());
        }
        size_t pivot_col = col_opt.value();
        auto row_opt = find_pivot_row(col_opt.value());
        if (!row_opt) {
            return Solution<T>(true);
        }
        size_t pivot_row = row_opt.value();
        update_pivot_row(pivot_row, pivot_col);
        update_constraints_values(pivot_row, pivot_col);
        update_other_rows(pivot_row, pivot_col);
        update_pivot_column(pivot_row, pivot_col);
        return Solution<T>(false, false, function_value(), vector_value());
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
        basic_coeffs.shrink_to_fit();
    }

    Solution<T> find_solution() {
        Solution<T> solution;
        do {
            solution = iterate();
        } while (!solution.is_final());
        return solution;
    }
};

#endif //ASSIGNMENT1_SIMPLEX_H
