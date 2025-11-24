#pragma once
#include <iostream>
#include <vector>
#include "Vector.hpp"
#include <math.h>
#include <complex>

template <typename K>
struct Matrix {
	std::vector<std::vector<K>> data;
	std::pair<size_t, size_t> shape() const;
	void add(const Matrix<K>& m);
	void sub(const Matrix<K>& m);
	void scl(K val);
	bool square() const;
	K trace() const;
	size_t size() const;
	bool empty() const;
	Vector<K> mul_vec(const Vector<K>& vec);
	Matrix<K> mul_mat(const Matrix<K>& vec);
	Matrix<K> transpose() const;
	Matrix<K> row_echelon() const;
	Matrix<K> row_echelon_det(size_t *swap) const;
	Matrix<K> reduced_row_echelon() const;
	K determinant() const;
	Matrix<K> inverse() const;
	size_t rank() const;
	Matrix<K> submatrix(size_t r, size_t c) const;
	std::vector<K>& operator[](size_t i) { return data[i]; }
	const std::vector<K>& operator[](size_t i) const { return data[i]; }
	Matrix<K> operator*(K f) const;
	Matrix<K> operator*(const Vector<K>& v) const;
	Matrix<K> operator*(const Matrix<K>& m) const;
	Matrix() = default;
	Matrix(std::initializer_list<std::initializer_list<K>> list) : data(list.begin(), list.end()) {};
	friend std::ostream& operator<<(std::ostream& os, const Matrix<K>& m) {
		for (size_t i = 0; i < m.data.size(); ++i)
		{
			os << "{";

			for (size_t j = 0; j < m.data[i].size(); ++j)
			{
				os << m.data[i][j];
				if (j < m.data[i].size() - 1)
					os << ", ";
			}
			os << "} ";
		}
		return os;
	}
};

template <typename K>
Matrix<K> identity(size_t n);

template <typename K>
Matrix<K> projection(float fov, float ratio, float near, float far);

constexpr double EPS = 1e-16;

template <typename K>
bool isZero(const K& x) {
	return std::abs(x) < EPS;
}

#include "Matrix.tpp"
