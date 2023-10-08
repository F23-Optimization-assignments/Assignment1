#ifndef ASSIGNMENT1_SIMPLEX_H
#define ASSIGNMENT1_SIMPLEX_H

#include <optional>
#include <algorithm>
#include "matrix.h"
#include "vector_ops.h"

#define MALFORMED_INPUT 0
#define DEGENERACY 1
#define ALTERNATIVE_OPTIMA 2
#define UNBOUNDED 3

class ValidationReport {
private:
    bool OK;
    // 1 - Degeneracy; 2 - Alternative Optima; 3 - Unbounded Solution;
    // 0 - Incorrect Input (Mismatch number of base vars.)
    int special_case_id;
    std::vector<size_t> incorrect_vars;

    friend std::ostream &operator<<(std::ostream &, const ValidationReport &);

public:
    ValidationReport(bool isOK,
                     int special_case_id,
                     const std::vector<size_t> &incorrect_vars) :
            OK(isOK), special_case_id(special_case_id), incorrect_vars(incorrect_vars) {}

    [[nodiscard]] bool is_OK() const {
        return OK;
    }

    [[nodiscard]] int special_case() const {
        return special_case_id;
    }

    [[nodiscard]] std::vector<size_t> get_incorrect_vars() const {
        return incorrect_vars;
    }
};

std::ostream &operator<<(std::ostream &stream, const ValidationReport &report) {
    if (report.is_OK()) {
        return stream;
    }

    if (report.special_case_id == MALFORMED_INPUT) {
        stream << ("There are incorrect base variables. Recheck the columns of variables.\n");

        stream << "Found basic variables: ";
        for (size_t incorrect_var: report.incorrect_vars) {
            stream << "X_" << incorrect_var << "; ";
        }
        stream << std::endl;

        return stream;
    }

    stream << ("[ SPECIAL CASE ]\n");

    switch (report.special_case_id) {
        case DEGENERACY:
            stream << ("There is a tie for a minimum ratio, which can increase the infinite loop.\n");
            break;
        case ALTERNATIVE_OPTIMA:
            stream << ("One of the constraints is in parallel with the objective function. "
                       "There are a lot of solutions (points) to the given problem.\n");
            break;
        case UNBOUNDED:
            stream << ("Can be increased or decreased infinitely (without violating "
                       "any constraint. An unbounded objective function.\n");
            break;
        default:
            stream << ("Something went wrong in the output.\n");
    }


    stream << "Affected variables:\n";
    for (size_t incorrect_var: report.incorrect_vars) {
        stream << "X_" << incorrect_var << "; ";
    }
    stream << std::endl;

    return stream;
}

template<typename T>
class Solution {
private:
    bool final;
    bool unbounded;
    T optimum;
    std::vector<T> xs;

    template<class U>
    friend std::ostream &operator<<(std::ostream &, const Solution<U> &);

public:
    Solution(bool isFinal,
             bool isUnbounded,
             T optimum,
             const std::vector<T> &xs) :
            final(isFinal), unbounded(isUnbounded), optimum(optimum), xs(xs) {}

    explicit Solution(const bool unbounded) : final(true), unbounded(unbounded) {}

    Solution() : final(false), unbounded(false) {};

    [[nodiscard]] bool is_final() const {
        return final;
    }

    [[nodiscard]] bool is_unbounded() const {
        return unbounded;
    }

    [[nodiscard]] std::vector<T> get_vector() const {
        return xs;
    }

    [[nodiscard]] T get_optimum() const {
        return optimum;
    }
};

template<class U>
std::ostream &operator<<(std::ostream &stream, const Solution<U> &solution) {
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
    bool is_correct_input;
    std::vector<T> func;
    Matrix<T> A;
    std::vector<T> b;
    std::vector<size_t> basic_indices;
    std::vector<T> basic_coeffs;

    void validate_on_start() {
        std::vector<size_t> checked_base_vars;

        for (size_t idx = 0; idx < A.get_columns(); idx++) {
            bool is_one_appears = false;
            bool is_only_zeroes = true;
            for (const T &num: A.get_column(idx)) {
                if (num == 1 && !is_one_appears) {
                    is_one_appears = true;
                    continue;
                }

                if (num != 0) {
                    is_only_zeroes = false;
                    break;
                }
            }

            if (is_only_zeroes && is_one_appears) checked_base_vars.push_back(idx);
        }

        if (checked_base_vars.size() < basic_indices.size()) {
            is_correct_input = false;
            std::cout << ValidationReport(is_correct_input, MALFORMED_INPUT, checked_base_vars);
            print_basic_indices();

            return;
        }

        std::sort(basic_indices.begin(), basic_indices.end());

        std::vector<size_t> incorrectVars;
        for (size_t base_var: basic_indices) {
            bool is_there = false;
            for (size_t num: checked_base_vars) {
                if (base_var == num) {
                    is_there = true;
                }
            }

            if (!is_there) {
                is_correct_input = false;
                incorrectVars.push_back(base_var);
            }
        }

        if (!is_correct_input) {
            std::cout << ValidationReport(is_correct_input, MALFORMED_INPUT, checked_base_vars);
            print_basic_indices();
        }
    }

    [[nodiscard]] std::vector<T> find_delta() const {
        std::vector<T> delta_row;
        for (size_t idx = 0; idx < func.size(); ++idx) {
            delta_row.push_back(dot_product(A.get_column(idx), basic_coeffs) - func[idx]);
        }
        return delta_row;
    }

    [[nodiscard]] std::optional<size_t> find_pivot_col(const std::vector<T> &v) const {
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
            if (A[row][col] == 0) {
                continue;
            }
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

    void update_pivot_row(size_t pivot_row, size_t pivot_col) {
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

    [[nodiscard]] std::vector<T> vector_value() const {
        std::vector<T> vars(func.size());
        for (size_t idx = 0; idx < basic_indices.size(); ++idx) {
            vars[basic_indices[idx]] = b[idx];
        }
        return vars;
    }

    void print_basic_indices() {
        std::cout << "Your basic variables: ";
        for (size_t inx: basic_indices) {
            std::cout << "X_" << inx << "; ";
        }
        std::cout << std::endl;
    }

    Solution<T> iterate() {
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
    Simplex(const std::vector<T> &coefficients,
            const Matrix<T> &A,
            const std::vector<T> &b,
            const std::vector<size_t> &basic_indices) : func(
            coefficients), A(A), b(b), basic_indices(
            basic_indices), basic_coeffs() {
        is_correct_input = true;
        validate_on_start();
        for (const size_t &idx: basic_indices) {
            basic_coeffs.push_back(func[idx]);
        }
        basic_coeffs.shrink_to_fit();
    }

    Solution<T> find_solution() {
        Solution<T> solution;

        if (!is_correct_input) {
            return Solution<T>();
        }

        do {
            solution = iterate();
        } while (!solution.is_final());
        return solution;
    }
};

#endif //ASSIGNMENT1_SIMPLEX_H
