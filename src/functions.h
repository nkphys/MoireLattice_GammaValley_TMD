#include "Matrix.h"
#include "tensor_type.h"


#ifndef define_functions
#define define_functions

Matrix<complex<double>> product(Matrix<complex<double>> A, Matrix<complex<double>> B){
    assert(A.n_col()==B.n_row());

    Matrix<complex<double>> C;
    C.resize(A.n_row(),B.n_col());

    for(int i=0; i<A.n_row(); i++) {
        for(int j=0; j<B.n_col(); j++) {

            C(i,j) = zero_complex;
            for(int k=0; k<B.n_row(); k++) {
            C(i,j) += A(i,k)*B(k,j);
            }

        }
    }

    return C;
} // ----------

#endif
