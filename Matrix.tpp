template <typename K>
void Matrix<K>::add(const Matrix<K>& m) {
	if (this->shape() != m.shape())
		throw std::invalid_argument("Matrice should have same dimensions");
	for(size_t i = 0; i < data.size(); ++i)
		for(size_t j = 0; j < data[i].size(); ++j)
			data[i][j] += m.data[i][j];
}

template <typename K>
void Matrix<K>::sub(const Matrix<K>& m) {
	if (this->shape() != m.shape())
		throw std::invalid_argument("Matrice should have same dimensions");
	for(size_t i = 0; i < data.size(); ++i)
		for(size_t j = 0; j < data[i].size(); ++j)
			data[i][j] -= m.data[i][j];
}

template <typename K>
void Matrix<K>::scl(K val) {
	for(size_t i = 0; i < data.size(); ++i)
		for(size_t j = 0; j < data[i].size(); ++j)
			data[i][j] *= val;
}

template <typename K>
std::pair<size_t, size_t> Matrix<K>::shape() const {
	size_t cols = data[0].size();
	size_t rows = data.size();
	for(size_t i = 0; i < rows; ++i)
	{
		if (data[i].size() != (size_t)cols)
			throw std::invalid_argument("All rows should have the same size");
	}
	return {rows, cols};
}

template <typename K>
size_t Matrix<K>::size() const {
	return data.size();
}

template <typename K>
bool Matrix<K>::empty() const {
	return data.empty();
}


template <typename K>
Vector<K> Matrix<K>::mul_vec(const Vector<K>& vec) {
	if (empty() || vec.size() != data[0].size())
		throw std::invalid_argument("Vector size does not match number of matrix columns");
	Vector<K> res;
	res.data.resize(size(), 0);
	for (size_t i = 0; i < size(); ++i) {
		for(size_t j = 0; j < data[i].size(); ++j) {
			res[i] += data[i][j] * vec[j];
		}
	}
	return res;
}

template <typename K>
Matrix<K> Matrix<K>::mul_mat(const Matrix<K>& vec) {
	auto [rows, cols] = shape();
	auto [rows2, cols2] = vec.shape();

	if (cols != rows2)
		throw std::invalid_argument("Matrix sizes error");
	Matrix<K> res;
	res.data.resize(rows,std::vector<K>(cols2, K{0}));

	for (size_t i = 0; i < rows; ++i) {
		for(size_t j = 0; j < cols2; ++j) {
			for (size_t k = 0; k < cols; ++k)
				res.data[i][j] += data[i][k] * vec.data[k][j];
		}
	}
	return res;
}

template <typename K>
bool Matrix<K>::square() const {
	auto [rows, cols] = shape();
	if (rows == cols)
		return true;
	return false;
}


template <typename K>
K Matrix<K>::trace() const {
	if (!square())
		throw std::invalid_argument("Matrix must be a square");
	K res = 0;
	for (size_t i = 0; i < size(); ++i)
		res += data[i][i];
	return res;
}

template <typename K>
Matrix<K> Matrix<K>::transpose() const {
	auto [rows, cols] = shape();

	Matrix<K> transpose;
	transpose.data.resize(cols, std::vector<K>(rows));
	for(size_t i = 0; i < rows; ++i)
		for(size_t j = 0; j < cols; ++j)
			transpose.data[j][i] = data[i][j];
	return transpose;
}

template <typename K>
Matrix<K> Matrix<K>::row_echelon() const {
	Matrix<K> res = *this;
	auto [rows, cols] = res.shape();
	size_t lead = 0;
	for (size_t r = 0; r < rows; ++r) {
		if (lead >= cols)
			break ;
		size_t i = r;
		while (i < rows && isZero(res[i][lead]))
			++i;
		if (i == rows) {
			++lead;
			--r;
			continue;
		}
		if (i != r)
			std::swap(res.data[i], res.data[r]);
		
		K pivot = res[r][lead];
		for (size_t j = lead; j < cols; ++j)
			res[r][j] /= pivot;
		for (size_t j = r + 1; j < rows; ++j) {
			K facteur = res[j][lead];
			for(size_t k = lead; k < cols; ++k)
				res[j][k] -= facteur * res[r][k];
			res[j][lead] = static_cast<K>(0);
		}
	}
	return res;
}

template <typename K>
Matrix<K> Matrix<K>::reduced_row_echelon() const {
	Matrix<K> res = row_echelon();
	auto [rows, cols] = res.shape();
	//RREF
	for (size_t i = 0; i < rows; ++i) {
		K pivot = static_cast<K>(0);
		size_t pivot_col = 0;
		for (size_t j = 0; j < cols; ++j) { //je cherche le pivot
			if (!isZero(res[i][j])) {
				pivot = res[i][j];
				pivot_col = j;
				break ;
			}
		}
		if (!isZero(pivot)) {  //je normalise tout ce qui est a droite du pivot
			for (size_t j = pivot_col; j < cols; ++j)
				res[i][j] /= pivot;
		}
	}
	for (int i = static_cast<int>(rows) - 1; i >= 0; --i) //je vais voir au dessus des pivot pour re-calculer la matrice
	{	
		size_t pivot_col = 0;
		while (pivot_col < cols && isZero(res[i][pivot_col]))
			++pivot_col;
		if (pivot_col == cols) //pas de pivots
			continue ;
		for (int j = i - 1; j >= 0; --j) {
			K facteur = res[j][pivot_col];
			if (!isZero(facteur)) {
				for(size_t k = pivot_col; k < cols; ++k)
					res[j][k] -= facteur * res[i][k];
			}
		}
	}
	return res;
}

template <typename K>
Matrix<K> Matrix<K>::row_echelon_det(size_t *swap) const {
	Matrix<K> res = *this;
	auto [rows, cols] = res.shape();
	size_t lead = 0;
	for (size_t r = 0; r < rows; ++r) {
		if (lead >= cols)
			break ;
		size_t i = r;
		while (i < rows && isZero(res[i][lead]))
			++i;
		if (i == rows) {
			++lead;
			--r;
			continue;
		}
		if (i != r)
		{
			std::swap(res.data[i], res.data[r]);
			(*swap)++;
		}
		K pivot = res[r][lead];
		for (size_t j = r + 1; j < rows; ++j) {
			K facteur = res[j][lead] / pivot;
			for(size_t k = lead; k < cols; ++k)
				res[j][k] -= facteur * res[r][k];
			res[j][lead] = static_cast<K>(0);
		}
	}
	return res;
}

template <typename K>
Matrix<K> identity(size_t n) {
	Matrix<K> res;
	res.data.resize(n);
	for (size_t i = 0; i < n; ++i)
		res.data[i].resize(n, 0);
	for (size_t i = 0; i < n; ++i)
		res[i][i] = 1;
	return res;
};


template <typename K>
K Matrix<K>::determinant() const {
	K det = 1;
	size_t swap = 0;
	auto [rows, cols] = shape();
	if (rows != cols || rows > 4)
		throw std::invalid_argument("This calculation only apply on a square matrix of 4x4 or less");
	Matrix<K> res = row_echelon_det(&swap);
	for (size_t i = 0; i < rows; i++) {
		if (isZero(res[i][i])) 
		{
			det = 0;
			break ;
		}
		det *= res[i][i];
	}
	if (swap % 2 != 0)
		det = -det;
	return det;
};

template <typename K> 
Matrix<K> Matrix<K>::submatrix(size_t r, size_t c) const {
	auto [rows, cols] = shape();
	if (rows <= 1 || cols <= 1)
    	throw std::invalid_argument("Cannot create submatrix for 1xN or Nx1");
	Matrix<K> res;
	res.data.resize(rows - 1);
	for (size_t i = 0; i < rows - 1; ++i)
		res.data[i].resize(cols - 1);
	size_t row_res = 0;
	for (size_t i = 0; i < rows; ++i) {
		if (i == r)
			continue ;
		size_t col_res = 0;
		for(size_t j = 0; j < cols; ++j) {
			if (j == c)
				continue ;
			res[row_res][col_res++] = (*this)[i][j];
		}
		row_res++;
	}
	return res;
};


template <typename K>
Matrix<K> Matrix<K>::inverse() const {
	if (!square())
		throw std::invalid_argument("Matrix must be a square");
	auto [rows, cols] = shape();
	Matrix<K> ident = identity<K>(rows);
	Matrix<K> res = *this;
	for (size_t i = 0; i < rows; ++i) {
		K pivot = res[i][i]; //recherche de pivot
		if (isZero(pivot)) // si le pivot est nul j'en cherche un non nul et je swap la ligne
		{
			bool found = false;
			for (size_t j = i + 1; j < rows; ++j){
				if (!isZero(res[j][i])) {
					std::swap(res[i], res[j]);
					std::swap(ident[i], ident[j]);
					pivot = res[i][i];
					found = true;
					break ;
				}
			}
			if (!found)
				throw std::invalid_argument("Singular Matrix");
		}
		for (size_t j = 0; j < rows; ++j) { // je normalise toute la ligne pour obtenir des 1.0 en diagonale
			res[i][j] /= pivot;
			ident[i][j] /= pivot;
		}
		for (size_t j = 0; j < rows; ++j) { //je vide toute la matrice sauf les pivot
			if (j == i)
				continue ;
			K factor = res[j][i];
			for (size_t k = 0; k < rows; ++k) {
				res[j][k] -= factor * res[i][k];
				ident[j][k] -= factor * ident[i][k];
			}
		}
	}
	return ident;
};

template <typename K>
size_t Matrix<K>::rank() const {
	auto [rows, cols] = shape();
	Matrix<K> RREF;
	size_t rank = 0;
	RREF = reduced_row_echelon();
	for (size_t i = 0; i < rows; ++i) {
		for (size_t j = 0; j < cols; ++j) { //Je parcourt la RREF jusqu'a trouver qqch different de 0 qui sera automatiquement mon pivot. Si je trouve qqch je change de ligne
			if (!isZero(RREF[i][j])) {
				rank++;
				break ;
			}
		}
	}
	return rank;
};

template <typename K>
Matrix<K> projection(float fov, float ratio, float near, float far) {
	K f = K(1) / std::tan(K(fov) / K(2));
	Matrix<K> proj({
		{f / ratio, K(0), K(0), K(0)},
		{K(0), f, K(0), K(0)},
		{K(0), K(0), -(far + near) / (far - near), K(1)},
		{K(0), K(0), -(2 * far * near) / (far - near), K(0)}
	});
	return proj;
};

template <typename K>
Matrix<K> Matrix<K>::operator*(K f) const {
	Matrix<K> res = *this;
	res.scl(f);
	return res;
};

template <typename K>
Matrix<K> Matrix<K>::operator*(const Matrix<K>& m) const {
	return mul_mat(m);
};

template <typename K>
Matrix<K> Matrix<K>::operator*(const Vector<K>& v) const {
	return mul_vec(v);
};

