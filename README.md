# A sparse matrix #
This is a simple  `c++` code to deal with sparse matrices.

The file  `Matrix.hpp` contains the implementation  of a template class  `Matrix<T, StorageOrder>`,  prototyping a sparse matrix that is able to switch from uncompressed to compressed state and viceversa. It also admits row-major and
column-major storage, according to the second template parameter (which is an enumerator).
In the case of the uncompressed state, the data is stored in COOMap format, the compressed
state employs CSR or CSC format, depending on the storage type.

The `main_matrix.cpp` does nothing much useful, apart from showing the main funcitonalities that have been implemented, such as:
- matrix-vector product of compressed and uncompressed matrices, both in a row-wise and column-wise representation;
- computing the `One`, `Infinity` and `Frobenius` norms of the matrix;
- reading the matrix from a `matrix market format`.

Moreover, it tests the computation time of the matrix-vector product, in all the possible ways mentioned above.
As a test case, it is used the matrix found in https://math.nist.gov/MatrixMarket/data/Harwell-Boeing/lns/lnsp_131.html.

## Build
To buid the program, it is enough to run the `make` command.

## Run
To run the main program, run the command `./main_matrix lnsp_131.mtx`
