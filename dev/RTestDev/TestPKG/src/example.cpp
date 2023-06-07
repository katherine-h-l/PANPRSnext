// [[Rcpp::depends(BH)]]

#include <Rcpp.h>
#include <iostream>
#include <sstream>
#include <boost/numeric/ublas/matrix_sparse.hpp>
#include <boost/multi_array.hpp>
#include <R_ext/Lapack.h>
#include <R_ext/Applic.h>
#include <R_ext/PrtUtil.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <R.h>
#include <Rmath.h>
// #include <omp.h>

namespace ublas = boost::numeric::ublas;

// [[Rcpp::export]]
void test_ublas()
{
    ublas::compressed_matrix<int> m1(10, 10); // 10x10 compressed matrix
    boost::multi_array<int, 2> m2(boost::extents[10][10]); // 10x10 multi_array

    std::istringstream in("1 1 4\n"
                          "1 3 5\n"
                          "2 1 6\n");

    // read from stream
    int val;
    for(size_t r, c; in >> r >> c >> val;) {
        m1(r, c) = val;
        m2[r][c] = val;
    }

    // print out
    for(size_t i = 0; i < m1.size1(); ++i)
    {
        for(size_t j = 0; j < m1.size2(); ++j)
             std::cout << m1(i,j) << ' ';
        std::cout << '\n';
    }

    // info on the storage
    std::cout << "Non-zeroes: " << m1.nnz() << '\n'
              << "Allocated storage for " << m1.nnz_capacity() << '\n';

    std::cout << "Memory usage for compressed matrix: " << sizeof(int) * m1.nnz_capacity() << " bytes\n";
    std::cout << "Memory usage for multi_array: " << sizeof(int) * m2.num_elements() << " bytes\n";
}
