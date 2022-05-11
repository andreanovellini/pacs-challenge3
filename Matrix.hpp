/*
 * Matrix.hpp
 *
 *  Created on: Apr 21, 2022
 *      Author: forma
 */

#ifndef MATRIX_MATRIX_HPP_
#define MATRIX_MATRIX_HPP_
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdlib>
#include <map>
#include <tuple>
#include <type_traits>
#include <vector>
#include <exception>
#include <numeric>
#include <algorithm>
#include <stdexcept>
#include <cassert>
#include <iostream>
#include "mmio.h"

namespace apsc
{
enum StorageOrder
{
  ROWMAJOR,
  COLUMNMAJOR
};

enum NormType
{
	ONE,
	INF,
	FRO
};

template <typename Scalar, StorageOrder Order=ROWMAJOR>
class Matrix
{
public:
	// Indexes are stored in an array of size 2 (unsigned)
  using COO = std::array<std::size_t, 2u>;

	// Default constructor
  Matrix() = default;

	// User-defined constructor
  Matrix(std::size_t nr, std::size_t nc) : nRows{nr}, nCols{nc} {};

	// Method to resize the matrix
  void resize(std::size_t nr,std::size_t nc);

	// Empties the matrix
  void reset();

	// Method to refactor the matrix in the compressed representation
  void makeCompressed();

	// Method to refactor the matrix in the uncompressed representation
  void decompress();

	// Number of non-zero elements
  std::size_t nnz() const;

	// Check if matrix is compressed or not
  bool is_compressed() const
	{
		return compressed_;
	}

	// Access operator (const version), read only!
  Scalar operator()(std::size_t i, std::size_t j) const;

	// Access operator (non-const version), read and write
  Scalar& operator()(std::size_t i,std::size_t j);
	
	// Overloading of the * operator: Matrix-vector product
	std::vector<Scalar> operator*(const std::vector<Scalar>& v) const;
	
	// Computes the norm of the matrix
	template<NormType normtype>	
	Scalar norm() const;
	
	// Reads matrix in matrix market format
	void readMatrix(const char* filename);	

	// Returns the inner size of the matrix (number of rows if ROWMAJOR, number of columns if COLUMNMAJOR)
  std::size_t nInner() const 
	{
		return Order == ROWMAJOR ? nRows : nCols;
	}

	// Returns the outer size of the matrix (number of columns if ROWMAJOR, number of rows if COLUMNMAJOR)
  std::size_t nOuter() const 
	{
		return Order == ROWMAJOR ? nCols : nRows;
	}
	
	// Returns vector of inner indexes
  auto getInnerIndexes() const 
	{
		return innerIndexes;
	}

	// Returns vector of outer indexes
  auto getOuterIndexes() const 
	{
		return outerIndexes;
	}

	// Returns the vector of values
  auto getValues() const 
	{
		return values;
	}
	
	// Prints the matrix	
	void print() const
	{
		for (auto i = 0u; i < nRows; ++i)
      {
      for (auto j = 0u; j < nCols; ++j)
          std::cout<<this->operator()(i,j)<<" ";
      std::cout<<std::endl;
      }
	}
	
	// Prints the vector of inner indexes
	void printInnerIndexes() const
	{
		if (!compressed_)
			return;

		for (const auto& i : innerIndexes)
			std::cout << i << " ";
	}

	// Prints the vector of outer indexes
  void printOuterIndexes() const 
	{
		if (!compressed_)
			return;

    for (const auto& i : outerIndexes) 
			std::cout << i << " ";
	}

	// Prints the vector of non-zero values
	void printValues() const
	{
    for (const auto& i : values) 
			std::cout << i << " ";
	}

	// Number of rows
  std::size_t nRows{0};
	// Number of columns
  std::size_t nCols{0};

private:
	// Returns the inner index based on the storage order of the matrix (ROWMAJOR or COLUMNMAJOR)
  static const std::size_t& getInner(const COO& c){
    if constexpr (Order == ROWMAJOR)
      return c[0];
    else
      return c[1];
  };

	// Returns the outer index based on the storage order of the matrix (ROWMAJOR or COLUMNMAJOR)
  static const std::size_t& getOuter(const COO& c){
    if constexpr (Order == ROWMAJOR)
      return c[1];
    else
      return c[0];
  };
	
	// Comparison method for the key-type of the map
  static constexpr auto compare = [](const COO& coo1, const COO& coo2) -> bool
  {
    return (std::tie(getInner(coo1), getOuter(coo1)) < std::tie(getInner(coo2),getOuter(coo2)));
  };
	
	// Data structure for the matrix (uncompressed). It' s a map: (i, j) -> A(i, j)
  using NCData = std::map<COO, Scalar, decltype(compare)>;
	
	// Boolean to indicate if the matrix is in the compressed or uncompressed form
  bool compressed_{false};
	// Vector of inner indexes
  std::vector<std::size_t> innerIndexes;
	// Vector of outer indexes
  std::vector<std::size_t> outerIndexes;
	// Vector of non zero values
  std::vector<Scalar> values;
  // Uncompressed data structure
	NCData ncData{compare};
};

template<typename Scalar, StorageOrder Order>
inline void Matrix<Scalar, Order>::resize(std::size_t nr, std::size_t nc)
{
	if (nr < nRows || nc < nCols)
	{
		throw std::runtime_error("I can only expand a Matrix with resize, sorry");
	}
	// I just have to change
	nRows = nr;
  nCols = nc;
  if (compressed_)
	{
		decompress();
		makeCompressed();
	}
}

template<typename Scalar, StorageOrder Order>
inline void Matrix<Scalar, Order>::reset()
{
	nRows = 0;
  nCols = 0;
  ncData.clear();
  innerIndexes.clear();
  outerIndexes.clear();
  values.clear();
  compressed_ = false;
}

template<typename Scalar, StorageOrder Order>
inline void Matrix<Scalar, Order>::makeCompressed()
{
	if (compressed_) 
		return;

	if constexpr (Order == ROWMAJOR)
	{
		innerIndexes.resize(nRows + 1u);
	} 
	else 
	{
		innerIndexes.resize(nCols + 1u);
	}

	std::fill(innerIndexes.begin(), innerIndexes.end(), 0u);
	outerIndexes.resize(nnz());
	values.resize(nnz());
	std::size_t counter{0u};
	for (const auto& [c, val] : ncData)
	{
		auto inner = getInner(c);
		auto outer = getOuter(c);
		++innerIndexes[inner + 1u];
		outerIndexes[counter] = outer;
		values[counter] = val;
		++counter;
	}

	// Fix inner
	for (auto i = 1u; i < innerIndexes.size(); ++i)
	{
		innerIndexes[i] += innerIndexes[i - 1u];
	}
	ncData.clear();
	compressed_ = true;
}

template<typename Scalar, StorageOrder Order>
inline void Matrix<Scalar, Order>::decompress()
{
	if (!compressed_) 
		return;
	for (auto in = 0u; in < nInner(); ++in)
	{
		for(auto k = innerIndexes[in]; k < innerIndexes[in + 1]; ++k)
		{
			auto out = outerIndexes[k];
			if constexpr(Order == ROWMAJOR)
				ncData.insert({{in,out}, values[k]});
			else
				ncData.insert({{out,in}, values[k]});
		}
	}
	innerIndexes.clear();
	outerIndexes.clear();
	values.clear();
	compressed_=false;
}

template<typename Scalar, StorageOrder Order>
inline Scalar  Matrix<Scalar, Order>::operator() (std::size_t i, std::size_t j) const
{
	assert(i<nRows && j<nCols);

	if(compressed_)
	{
		std::size_t inner = getInner({i, j});
		std::size_t outer = getOuter({i, j});
		std::size_t offset = innerIndexes[inner];
		std::size_t last = innerIndexes[inner + 1u];
		auto res = std::find(outerIndexes.begin() + offset, outerIndexes.begin() + last, outer);
		if(res == outerIndexes.begin() + last)
		{
			return Scalar{0};
		}
		else
		{
			return values[std::distance(outerIndexes.begin(), res)];
		}
	}
	else
	{
		auto res = ncData.find({i, j});
		if (res == ncData.end())
			return Scalar{0};
		else
		{
			return res->second;
		}
	}
}

template<typename Scalar, StorageOrder Order>
inline Scalar& Matrix<Scalar, Order>::operator()(std::size_t i, std::size_t j)
{
	if (compressed_)
	{
		assert(i < nRows && j < nCols);
		std::size_t inner = getInner({i, j});
		std::size_t outer = getOuter({i, j});
		std::size_t offset = innerIndexes[inner];
		std::size_t last = innerIndexes[inner + 1u];
		auto res = std::find(outerIndexes.begin() + offset, outerIndexes.begin() + last, outer);
		if (res == outerIndexes.begin() + last)
		{
			throw std::runtime_error("Cannot add element in compressed state");
		}
		else
		{
			return values[std::distance(outerIndexes.begin(), res)];
		}
	}
	else
	{
		auto res = ncData.find({i, j});
		if (res == ncData.end())
		{
			nRows = std::max(nRows, i + 1u);
			nCols = std::max(nCols,j + 1u);
			auto pos = ncData.insert({{i, j}, Scalar{0}});
			return pos.first->second;
		}
		else
		{
			return res->second;
		}
	}
}

template<typename Scalar, StorageOrder Order>
inline std::vector<Scalar> Matrix<Scalar, Order>::operator*(const std::vector<Scalar> &v) const
{
	if (v.size() == nCols)
	{
		std::vector<Scalar> ret(nRows, 0.0);
		
		if (compressed_)
		{
			for (auto in = 0u; in < nInner(); ++in)
			{
				for(auto k = innerIndexes[in]; k < innerIndexes[in + 1]; ++k)
				{
					auto out = outerIndexes[k];
					if constexpr(Order == ROWMAJOR)
						ret[in] += values[k] * v[out];
					else
						ret[out] += values[k] * v[in];
				}
			}
		}
		else
		{
			for (const auto& [c, val] : ncData)
				ret[c[0]] += val * v[c[1]];
		}	
		return ret;
	}
	else
		throw std::runtime_error("Incompatible dimensions");
}
template<typename Scalar, StorageOrder Order>
template<NormType normtype>
inline Scalar Matrix<Scalar, Order>::norm() const
{
	double norm = 0.0;
	if constexpr (normtype == ONE)
	{
		std::vector<double> ColSums(nCols, 0.0);
		if (compressed_)
		{
			for (auto in = 0u; in < nInner(); ++in)
			{
				for(auto k = innerIndexes[in]; k < innerIndexes[in + 1]; ++k)
				{
					auto out = outerIndexes[k];
					if constexpr(Order == ROWMAJOR)
						ColSums[out] += std::abs(values[k]); 
					else
						ColSums[in] += std::abs(values[k]);
				}
			}
		}
		else
		{
			for (const auto& [c, val] : ncData)
				ColSums[c[1]] += std::abs(val);
		}
		norm = *std::max_element(ColSums.cbegin(), ColSums.cend());
	}
	else if constexpr (normtype == INF)
	{
		std::vector<double> RowSums(nRows, 0.0);
		if (compressed_)
		{
			for (auto in = 0u; in < nInner(); ++in)
			{
				for(auto k = innerIndexes[in]; k < innerIndexes[in + 1]; ++k)
				{
					auto out = outerIndexes[k];
					if constexpr(Order == ROWMAJOR)
						RowSums[in] += std::abs(values[k]); 
					else
						RowSums[out] += std::abs(values[k]);
				}
			}
		}
		else
		{
			for (const auto& [c, val] : ncData)
				RowSums[c[0]] += std::abs(val);
		}
		norm = *std::max_element(RowSums.cbegin(), RowSums.cend());
	}
	else if constexpr (normtype == FRO)
	{
		if (compressed_)
			for (const auto& el : values)
				norm += el * el;
		else 
			for (const auto& [c, val] : ncData)
				norm += val * val;
		norm = std::sqrt(norm);
	}
	else
		throw std::runtime_error("Invalid norm type");

	return norm;
}

template<typename Scalar, StorageOrder Order>
inline void Matrix<Scalar, Order>::readMatrix(const char* filename)
{
	int ret_code;
	MM_typecode matcode;
	FILE* f;
	int M, N, nz;
	int i, *I, *J;
	double* val;

	if ((f = fopen(filename, "r")) == NULL)
		throw std::runtime_error("Could not open the file");

	if (mm_read_banner(f, &matcode) != 0)
		throw std::runtime_error("Could not read banner");

  /*  This is how one can screen matrix types if their application */
  /*  only supports a subset of the Matrix Market data types.      */
	if (mm_is_complex(matcode) && mm_is_matrix(matcode) && mm_is_sparse(matcode) )
			throw std::runtime_error("Sorry, this application does not support this Matrix Market type");

  /* find out size of sparse matrix .... */
  if ((ret_code = mm_read_mtx_crd_size(f, &M, &N, &nz)) !=0)
		throw std::runtime_error("Error reading the matrix");

	// Initialize the size of the matrix
	nRows = M;
	nCols = N;
			
	/* reseve memory for matrices */
  I = (int *) malloc(nz * sizeof(int));
  J = (int *) malloc(nz * sizeof(int));
  val = (double *) malloc(nz * sizeof(double));
	
  /* NOTE: when reading in doubles, ANSI C requires the use of the "l"  */
  /*   specifier as in "%lg", "%lf", "%le", otherwise errors will occur */
  /*  (ANSI C X3.159-1989, Sec. 4.9.6.2, p. 136 lines 13-15)            */
  for (i=0; i<nz; i++)
  {
      fscanf(f, "%d %d %lg\n", &I[i], &J[i], &val[i]);
      I[i]--;  /* adjust from 1-based to 0-based */
      J[i]--;

			ncData.insert({{I[i], J[i]}, val[i]});
  }

  if (f != stdin) fclose(f);
	
	return;
}

template<typename Scalar, StorageOrder Order>
inline std::size_t Matrix<Scalar, Order>::nnz() const
{
	if (compressed_)
		return outerIndexes.size();
	else
		return ncData.size();
}
}// end namespace apsc
#endif /* MATRIX_MATRIX_HPP_ */
