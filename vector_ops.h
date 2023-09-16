//
// Created by Matvey Korinenko on 16.09.2023.
//

#ifndef ASSIGNMENT1_VECTOR_OPS_H
#define ASSIGNMENT1_VECTOR_OPS_H

#include <iostream>
#include <vector>

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

#endif //ASSIGNMENT1_VECTOR_OPS_H
