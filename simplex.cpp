#include <iostream>
#include <vector>
#include <cmath>

class Fraction {
private:
    int p, q;

    [[nodiscard]] int gcd(int a, int b) const {
        return (b == 0 ? a : gcd(b, a % b));
    }

    [[ nodiscard ]] int lcm(int a, int b) const {
        return a / gcd(a, b) * b;
    }

    void simplify() {
        bool sign = (p < 0) ^ (q < 0);
        p = std::abs(p);
        q = std::abs(q);
        int g = gcd(p, q);
        p = p / g;
        q = q / g;
        p *= (sign ? -1 : 1);
    }

    friend std::ostream& operator<<(std::ostream& stream, const Fraction& fraction);
    friend std::istream& operator>>(std::istream& stream, Fraction& fraction);

public:
    Fraction(int p, int q) : p(p), q(q) { simplify(); }
    explicit Fraction(int p) : Fraction(p, 1) { }
    Fraction() : Fraction(0) { }

    Fraction& operator*=(const Fraction& other) {
        p *= other.p;
        q *= other.q;
        simplify();
        return *this;
    }

    Fraction operator*(const Fraction& other) const {
        Fraction f(*this);
        return f *= other;
    }

    Fraction& operator/=(const Fraction& other) {
        p *= other.q;
        q *= other.p;
        simplify();
        return *this;
    }

    Fraction operator/(const Fraction& other) const {
        Fraction f(*this);
        return f /= other;
    }

    Fraction& operator+=(const Fraction& other) {
        int l = lcm(q, other.q);
        int mp = l / q;
        p *= mp;
        q *= mp;
        p += (l / other.q) * other.p;
        simplify();
        return *this;
    }

    Fraction operator+(const Fraction& other) const {
        Fraction f(*this);
        return f += other;
    }

    Fraction& operator-=(const Fraction& other) {
        int l = lcm(q, other.q);
        int mp = l / q;
        p *= mp;
        q *= mp;
        p -= (l / other.q) * other.p;
        simplify();
        return *this;
    }

    Fraction operator-(const Fraction& other) const {
        Fraction f(*this);
        return f -= other;
    }

    Fraction& operator=(const int& other) {
        return *this = Fraction(other);
    }

    bool operator<(const Fraction& other) const {
        int l = lcm(q, other.q);
        int lp = p * (l / q), rp = other.p * (l / other.q);
        return lp < rp;
    }

    bool operator>(const Fraction& other) const {
        return other < *this;
    }

    bool operator==(const Fraction& other) const {
        return !(*this < other || *this > other);
    }

    bool operator<=(const Fraction& other) const {
        return !(*this > other);
    }

    bool operator>=(const Fraction& other) const {
        return !(*this < other);
    }

    bool operator<(const int& other) const {
        return *this < Fraction(other);
    }

    bool operator<=(const int& other) const {
        return *this <= Fraction(other);
    }

    bool operator==(const int& other) const {
        return *this == Fraction(other);
    }

    bool operator>(const int& other) const {
        return *this > Fraction(other);
    }

    bool operator>=(const int& other) const {
        return *this >= Fraction(other);
    }
};

std::ostream& operator<<(std::ostream& stream, const Fraction& fraction) {
    if (fraction.q == 1) {
        return stream << fraction.p;
    }
    return stream << fraction.p << '/' << fraction.q;
}

std::istream& operator>>(std::istream& stream, Fraction& fraction) {
    stream >> fraction.p;
    fraction.q = 1;
    fraction.simplify();
    return stream;
}


class VectorException : public std::exception {
private:
    const char* msg;
public:
    explicit VectorException(const char* msg) : std::exception(), msg(msg) { }
    [[nodiscard]] const char * what() const noexcept override {
            return msg;
    }
};

template<typename T>
T dot_product(const std::vector<T>& a, const std::vector<T>& b) {
    if (a.size() != b.size()) {
        throw VectorException("incompatible vectors' sizes for dot product");
    }
    T res;
    for (size_t idx = 0; idx < a.size(); ++idx) {
        res += a[idx] * b[idx];
    }
    return res;
}


class SimplexException : public std::exception {
private:
    const char* msg;
public:
    explicit SimplexException(const char* msg) : std::exception(), msg(msg) { }
    [[nodiscard]] const char * what() const noexcept override {
            return msg;
    }
};


class MatrixException : public std::exception {
private:
    const char* msg;
public:
    explicit MatrixException(const char* msg) : std::exception(), msg(msg) { }
    [[nodiscard]] const char * what() const noexcept override {
            return msg;
    }
};

template<typename T>
class Matrix {
private:
    size_t rows, columns;
    std::vector<std::vector<T>> entries;

    void check_sum(const Matrix<T>& other) const noexcept {
        if (rows != other.rows || columns != other.columns) {
            throw MatrixException("incompatible matrices' dimensions to sum/subtract them");
        }
    }

    void check_product(const Matrix<T>& other) const {
        if (columns != other.rows) {
            throw MatrixException("incompatible matrices' dimensions to multiply them");
        }
    }

    template<typename U>
    friend std::istream& operator>>(std::istream& stream, Matrix<U>& matrix);

public:
    explicit Matrix(const size_t& n) : rows(n), columns(n), entries(n, std::vector<T>(n)) { }
    Matrix(const size_t& n, const size_t& m) : rows(n), columns(m), entries(n, std::vector<T>(m)) { }

    [[ nodiscard ]] size_t get_rows() const {
        return rows;
    }

    [[ nodiscard ]] size_t get_columns() const {
        return columns;
    }

    [[ nodiscard ]] std::vector<T> get_row(size_t i) const {
        return entries[i];
    }

    [[nodiscard ]] std::vector<T> get_column(size_t i) const {
        std::vector<T> column;
        for (size_t j = 0; j < rows; ++j) {
            column.push_back(entries[j][i]);
        }
        return column;
    }

    std::vector<T>& operator[](size_t i) {
        return entries[i];
    }

    const std::vector<T>& operator[] (size_t i) const {
        return entries[i];
    }

    Matrix<T>& operator+=(const Matrix<T>& other) {
        check_sum(other);
        for (size_t i = 0; i < rows; ++i) {
            for (size_t j = 0; j < columns; ++j) {
                entries[i][j] += other[i][j];
            }
        }
        return *this;
    }

    Matrix<T> operator+(const Matrix<T>& other) const {
        Matrix<T> sum = *this;
        sum += other;
        return sum;
    }

    Matrix<T>& operator-=(const Matrix<T>& other) {
        check_sum();
        for (size_t i = 0; i < rows; ++i) {
            for (size_t j = 0; j < columns; ++j) {
                entries[i][j] -= other[i][j];
            }
        }
        return *this;
    }

    Matrix<T> operator-(const Matrix<T>& other) const {
        Matrix<T> sum = *this;
        sum -= other;
        return sum;
    }

    Matrix<T> operator*(const Matrix<T>& other) const {
        check_product(other);
        Matrix<T> product(rows, other.columns);
        for (size_t i = 0; i < rows; ++i) {
            for (size_t j = 0; j < other.columns; ++j) {
                for (size_t k = 0; k < other.rows; ++k) {
                    product[i][j] += entries[i][k] * other[k][j];
                }
            }
        }
        return product;
    }

    Matrix<T>& operator*=(const Matrix<T>& other) {
        return *this = (*this * other);
    }

    bool operator==(const Matrix<T>& other) const {
        if (rows != other.rows || columns != other.columns) {
            return false;
        }
        for (size_t i = 0; i < rows; ++i) {
            for (size_t j = 0; j < columns; ++j) {
                if (entries[i][j] != other[i][j]) {
                    return false;
                }
            }
        }
        return true;
    }

    bool operator!=(const Matrix<T>& other) const {
        return *this != other;
    }
};

template<typename T>
std::istream& operator>>(std::istream& stream, Matrix<T>& matrix) {
    for (size_t i = 0; i < matrix.rows; ++i) {
        for (size_t j = 0; j < matrix.columns; ++j) {
            stream >> matrix[i][j];
        }
    }
    return stream;
}

template<typename T>
std::ostream& operator<<(std::ostream& stream, const Matrix<T>& matrix) {
    for (size_t i = 0; i < matrix.get_rows(); ++i) {
        for (size_t j = 0; j < matrix.get_columns(); ++j) {
            stream << matrix[i][j] << ' ';
        }
        stream << ' ';
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
        basic_indices[pivot_col] = pivot_col;
        basic_coeffs[pivot_col] = func[pivot_col];
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
            std::cout << "Variable #" << basic_indices[idx] << " is " << b[idx] << '\n';
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