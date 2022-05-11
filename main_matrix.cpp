/*
 * main_matrix.cpp
 *
 *  Created on: Apr 22, 2022
 *      Author: forma
 */
#include <chrono>
#include <iostream>
#include <ostream>
#include "Matrix.hpp"

using timer = std::chrono::high_resolution_clock;

int main(int argc, char** argv)
{
  using namespace apsc;
	{
		std::cout << "***** MATRIX MARKET EXAMPLE *****" << std::endl;
		if (argc != 2)
		{
			std::cout << "Only the name of the file is needed." << std::endl;
			return 0;
		}
		else 
		{
			auto start = timer::now();
			auto end = timer::now();
			{
				// ROWMAJOR UNCOMPRESSED
				std::cout << "\n... Row-wise representation, uncompressed ...\n" << std::endl;
				Matrix<double, ROWMAJOR> A;
				A.readMatrix(argv[1]);
				// Print matrix dimensions
				std::cout << "nrows = " << A.nRows << ", ncols = " << A.nCols << std::endl;
				std::cout << "Number of non-zero elements: " << A.nnz() << std::endl;
				// Calculating norms
				std::cout << "One norm: " << A.norm<ONE>() << std::endl;
				std::cout << "Infinity norm: " << A.norm<INF>() << std::endl;
				std::cout << "Frobenius norm: " << A.norm<FRO>() << std::endl;
				// Matrix-vector product
				std::vector<double> v(A.nCols, 2.5);
				std::vector<double> result; 
				std::cout << "Computing A*v ..." << std::endl;
				start = timer::now();
				result = A * v;	
				end = timer::now();
				auto elapsed_time = std::chrono::duration_cast<std::chrono::microseconds>(end-start).count();
				std::cout << "Elapsed time to compute the product: " << elapsed_time << " microseconds" << std::endl;

				// ROWMAJOR COMPRESSED
				std::cout << "\n... Row-wise representation, compressed ...\n" << std::endl;
				// Compressing
				std::cout << "Compressing the matrix ..." << std::endl;
				A.makeCompressed();
				// Print matrix dimensions
				std::cout << "nrows = " << A.nRows << ", ncols = " << A.nCols << std::endl;
				std::cout << "Number of non-zero elements: " << A.nnz() << std::endl;
				// Calculating norms
				std::cout << "One norm: " << A.norm<ONE>() << std::endl;
				std::cout << "Infinity norm: " << A.norm<INF>() << std::endl;
				std::cout << "Frobenius norm: " << A.norm<FRO>() << std::endl;
				// Matrix-vector product
				std::cout << "Computing A*v ..." << std::endl;
				start = timer::now();
				result = A * v;	
				end = timer::now();
			  elapsed_time = std::chrono::duration_cast<std::chrono::microseconds>(end-start).count();
				std::cout << "Elapsed time to compute the product: " << elapsed_time << " microseconds" << std::endl;	
			}
			{
				// COLUMNMAJOR UNCOMPRESSED
				std::cout << "\n... Column-wise representation, uncompressed ...\n" << std::endl;
				Matrix<double, COLUMNMAJOR> A;
				A.readMatrix(argv[1]);
				// Print matrix dimensions
				std::cout << "nrows = " << A.nRows << ", ncols = " << A.nCols << std::endl;
				std::cout << "Number of non-zero elements: " << A.nnz() << std::endl;
				// Calculating norms
				std::cout << "One norm: " << A.norm<ONE>() << std::endl;
				std::cout << "Infinity norm: " << A.norm<INF>() << std::endl;
				std::cout << "Frobenius norm: " << A.norm<FRO>() << std::endl;
				// Matrix-vector product
				std::vector<double> v(A.nCols, 2.5);
				std::vector<double> result; 
				std::cout << "Computing A*v ..." << std::endl;
				auto start = timer::now();
				result = A * v;	
				auto end = timer::now();
				auto elapsed_time = std::chrono::duration_cast<std::chrono::microseconds>(end-start).count();
				std::cout << "Elapsed time to compute the product: " << elapsed_time << " microseconds" << std::endl;

				// COLUMNMAJOR COMPRESSED
				std::cout << "\n... Column-wise representation, compressed ...\n" << std::endl;
				// Compressing
				std::cout << "Compressing the matrix ..." << std::endl;
				A.makeCompressed();
				// Print matrix dimensions
				std::cout << "nrows = " << A.nRows << ", ncols = " << A.nCols << std::endl;
				std::cout << "Number of non-zero elements: " << A.nnz() << std::endl;
				// Calculating norms
				std::cout << "One norm: " << A.norm<ONE>() << std::endl;
				std::cout << "Infinity norm: " << A.norm<INF>() << std::endl;
				std::cout << "Frobenius norm: " << A.norm<FRO>() << std::endl;
				// Matrix-vector product
				std::cout << "Computing A*v ..." << std::endl;
				start = timer::now();
				result = A * v;	
				end = timer::now();
			  elapsed_time = std::chrono::duration_cast<std::chrono::microseconds>(end-start).count();
				std::cout << "Elapsed time to compute the product: " << elapsed_time << " microseconds" << std::endl;	

			}
		}
	}
	std::cout << "\n***** ROWMAJOR EXAMPLE *****" << std::endl;	
  {	
		std::cout << "\n... Simple row-wise example, starting from an uncompressed matrix...\n" << std::endl;
		// Creating the matrix
		std::cout << "Creating a 4 by 5 matrix ..." << std::endl;
    auto n = 4;
    auto m = 5;
    Matrix<double,ROWMAJOR> A(n, m);
		// Inserting two elements
    A(0, 0) = 1;
    A(3, 2) = -1;
		// Checking number of rows and number of columns
    std::cout << "nrows = " << A.nRows << ", ncols = " << A.nCols << std::endl;
		// Printing	
		std::cout << "Printing the matrix ..." << std::endl;
		A.print();
		std::cout << std::endl;	
		// Matrix-vector product (ROWMAJOR, uncompressed)
		std::cout << "Computing the matrix-vector product, A*v ..." << std::endl;
		std::cout << "Matrix A:" << std::endl;
		A.print();
		std::cout << std::endl;
		std::cout << "Vector v:" << std::endl;
		std::vector<double> v{1, 0, 4, 2, 3};
		for (const auto& el : v)
			std::cout << el << std::endl;
		std::cout << std::endl;
		std::cout << "Result:" << std::endl;
		auto result = A * v;
		for (const auto& el : result)
			std::cout << el << std::endl;
		std::cout << std::endl;
		// Compressing 
		std::cout << "Compressing the matrix ..." << std::endl;
    A.makeCompressed();
		// Printing compressed matrix relevant informations
		std::cout << "nrows = " << A.nRows << ", ncols = " << A.nCols << std::endl;
    std::cout << "Number of non-zero elements = " << A.nnz() << std::endl;
		// Printing inner indexes
		std::cout << "Inner indexes: ";
		A.printInnerIndexes();
		std::cout << std::endl;
		// Printing outer indexes
    std::cout << "Outer indexes: ";
		A.printOuterIndexes();
		std::cout << std::endl;
  	// Printing the values
		std::cout << "Values: ";
		A.printValues();
		std::cout << std::endl;
		// Matrix-vector product (ROWMAJOR, compressed)
		std::cout << "Computing the (compressed)matrix-vector product, A*v ..." << std::endl;
		std::cout << "Matrix A:" << std::endl;
		A.print();
		std::cout << std::endl;
		std::cout << "Vector v: " << std::endl;
		// std::vector<double> v{1, 0, 4, 2, 3}; // Not needed anymore
		for (const auto& el : v)
			std::cout << el << std::endl;
		std::cout << std::endl;
		std::cout << "Result: " << std::endl;
		result = A * v;
		for (const auto& el : result)
			std::cout << el << std::endl;
		std::cout << std::endl;
		// Decompressing the matrix	
		std::cout << "Decompressing the matrix ..." << std::endl;
    A.decompress();
		// Adding a new value (changes the matrix dimensions)
		std::cout << "Adding a new value ..." << std::endl;
    A(6, 6) = 10;
		// Printing the matrix
		std::cout << "The updated matrix is: " << std::endl;
		A.print();
		std::cout << std::endl;
		// Recompressing the matrix
		std::cout << "Recompressing the matrix: " << std::endl;
    A.makeCompressed();
		// Printing compressed matrix relevant informations
    std::cout << "nrows = " << A.nRows << ", ncols = " << A.nCols << std::endl;
    std::cout << "Number of non-zero elements = " << A.nnz() << std::endl;
		// Printing inner indexes
    std::cout << "Inner indexes: ";
		A.printInnerIndexes();
		std::cout << std::endl;
		// Printing outer indexes
    std::cout << "Outer indexes: ";
		A.printOuterIndexes();
		std::cout << std::endl;
		// Printing the values
    std::cout << "Values: ";
		A.printValues();
		std::cout << std::endl;
	}
	std::cout << "\n***** COLUMNMAJOR EXAMPLE *****" << std::endl;
  {	
		std::cout << "\n... Simple column-wise example, starting from an uncompressed matrix...\n" << std::endl;
		// Creating the matrix
		std::cout << "Creating a 4 by 5 matrix ..." << std::endl;
    auto n = 4;
    auto m = 5;
    Matrix<double,COLUMNMAJOR> A(n, m);
		// Inserting two elements
    A(0, 0) = 1;
    A(3, 2) = -1;
		// Checking number of rows and number of columns
    std::cout << "nrows = " << A.nRows << ", ncols = " << A.nCols << std::endl;
		// Printing	
		std::cout << "Printing the matrix ..." << std::endl;
		A.print();
		std::cout << std::endl;	
		// Matrix-vector product (COLUMNMAJOR, uncompressed)
		std::cout << "Computing the matrix-vector product, A*v ..." << std::endl;
		std::cout << "Matrix A:" << std::endl;
		A.print();
		std::cout << std::endl;
		std::cout << "Vector v:" << std::endl;
		std::vector<double> v{1, 0, 4, 2, 3};
		for (const auto& el : v)
			std::cout << el << std::endl;
		std::cout << std::endl;
		std::cout << "Result:" << std::endl;
		auto result = A * v;
		for (const auto& el : result)
			std::cout << el << std::endl;
		std::cout << std::endl;
		// Compressing 
		std::cout << "Compressing the matrix ..." << std::endl;
    A.makeCompressed();
		// Printing compressed matrix relevant informations
		std::cout << "nrows = " << A.nRows << ", ncols = " << A.nCols << std::endl;
    std::cout << "Number of non-zero elements = " << A.nnz() << std::endl;
		// Printing inner indexes
		std::cout << "Inner indexes: ";
		A.printInnerIndexes();
		std::cout << std::endl;
		// Printing outer indexes
    std::cout << "Outer indexes: ";
		A.printOuterIndexes();
		std::cout << std::endl;
  	// Printing the values
		std::cout << "Values: ";
		A.printValues();
		std::cout << std::endl;
		// Matrix-vector product (COLUMNMAJOR, compressed)
		std::cout << "Computing the (compressed)matrix-vector product, A*v ..." << std::endl;
		std::cout << "Matrix A:" << std::endl;
		A.print();
		std::cout << std::endl;
		std::cout << "Vector v: " << std::endl;
		// std::vector<double> v{1, 0, 4, 2, 3}; // Not needed anymore
		for (const auto& el : v)
			std::cout << el << std::endl;
		std::cout << std::endl;
		std::cout << "Result: " << std::endl;
		result = A * v;
		for (const auto& el : result)
			std::cout << el << std::endl;
		std::cout << std::endl;
		// Decompressing the matrix	
		std::cout << "Decompressing the matrix ..." << std::endl;
    A.decompress();
		// Adding a new value (changes the matrix dimensions)
		std::cout << "Adding a new value ..." << std::endl;
    A(6, 6) = 10;
		// Printing the matrix
		std::cout << "The updated matrix is: " << std::endl;
		A.print();
		std::cout << std::endl;
		// Recompressing the matrix
		std::cout << "Recompressing the matrix: " << std::endl;
    A.makeCompressed();
		// Printing compressed matrix relevant informations
    std::cout << "nrows = " << A.nRows << ", ncols = " << A.nCols << std::endl;
    std::cout << "Number of non-zero elements = " << A.nnz() << std::endl;
		// Printing inner indexes
    std::cout << "Inner indexes: ";
		A.printInnerIndexes();
		std::cout << std::endl;
		// Printing outer indexes
    std::cout << "Outer indexes: ";
		A.printOuterIndexes();
		std::cout << std::endl;
		// Printing the values
    std::cout << "Values: ";
		A.printValues();
		std::cout << std::endl;
	}
	return 0;
}


