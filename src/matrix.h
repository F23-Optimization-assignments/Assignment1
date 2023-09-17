#ifndef ASSIGNMENT1_MATRIX_H
#define ASSIGNMENT1_MATRIX_H

#include <iostream>
#include <vector>

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

#endif //ASSIGNMENT1_MATRIX_H
