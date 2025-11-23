#pragma once
#include <iostream>
#include <vector>
#include <math.h>
#include <complex>


template <typename K>
struct Vector {
	std::vector<K> data;
	void add(const Vector<K>& v);
	void sub(const Vector<K>& v);
	void scl(const K& val);
	float norm_1() const;
	float norm() const;
	float norm_inf() const;
	size_t size() const;
	bool empty() const;
	K dot(const Vector<K>& v) const;
	K& operator[](size_t i) { return data[i]; }
	const K& operator[](size_t i) const { return data[i]; }
	Vector() = default;
	//variadic template constructor for multiple args
	template <typename... Args>
	Vector(Args... args) : data{static_cast<K>(args)...} {};
	Vector(std::initializer_list<K> list) : data{list} {};
	friend std::ostream& operator<<(std::ostream& os, const Vector<K>& v) {
		os << "{";
		for (size_t i = 0; i < v.data.size(); ++i)
		{
			os << v.data[i];
			if (i < v.data.size() - 1)
				os << ", ";
		}
		os << "}";
		return os;
	}
};
template <typename K>
float angle_cos(const Vector<K>&u, const Vector<K>&v);

template <typename K>
Vector<K> cross_product(const Vector<K>&u, const Vector<K>&v);

template <typename K>
Vector<K> linear_combination(const std::vector<Vector<K>>& u, const std::vector<K>& coefs);

#include "Vector.tpp"