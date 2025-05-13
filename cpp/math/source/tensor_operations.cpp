#include "tensor_operations.hpp"

namespace math {


/* Cross product of two tensors of the second rank */
void cross(const vector_t<double,2>& A
          ,const vector_t<double,2>& B
               , vector_t<double,3>& res)
{
    // loop over 3rd dimension
    for (size_t k = 0; k < 3; ++k) {
        for (size_t i = 0; i < 3; ++i) {
            res[k][0][i] = A[k][1]*B[2][i] - A[k][2]*B[1][i];
            res[k][1][i] = A[k][2]*B[0][i] - A[k][0]*B[2][i];
            res[k][2][i] = A[k][0]*B[1][i] - A[k][1]*B[0][i];
        }
    }
}

/* Cross product of transpose tensor and ordinary tensor of the second rank */
void crossT0(const vector_t<double,2>& AT
            ,const vector_t<double,2>& B
                , vector_t<double,3>& res)
{
    // loop over 3rd dimension
    for (size_t k = 0; k < 3; ++k) {
        for (size_t i = 0; i < 3; ++i) {
            res[k][0][i] = AT[1][k]*B[2][i] - AT[2][k]*B[1][i];
            res[k][1][i] = AT[2][k]*B[0][i] - AT[0][k]*B[2][i];
            res[k][2][i] = AT[0][k]*B[1][i] - AT[1][k]*B[0][i];
        }
    }
}

/* Cross product of transpose tensor and ordinary tensor of the second rank */
void crossT0(const vector_t<double,2>& AT
    ,const vector<double>& b
        , vector_t<double,2>& res)
{
    // loop over 3rd dimension
    for (size_t k = 0; k < 3; ++k) {
        res[k][0] = AT[1][k]*b[2] - AT[2][k]*b[1];
        res[k][1] = AT[2][k]*b[0] - AT[0][k]*b[2];
        res[k][2] = AT[0][k]*b[1] - AT[1][k]*b[0];
    }
}

/* Cross product of transpose tensor and ordinary tensor of the second rank */
void cross0T(const vector_t<double,2>& A
            ,const vector_t<double,2>& BT
                , vector_t<double,3>& res)
{
    // loop over 3rd dimension
    for (size_t k = 0; k < 3; ++k) {
        for (size_t i = 0; i < 3; ++i) {
            res[k][0][i] = A[k][1]*BT[i][2] - A[k][2]*BT[i][1];
            res[k][1][i] = A[k][2]*BT[i][0] - A[k][0]*BT[i][2];
            res[k][2][i] = A[k][0]*BT[i][1] - A[k][1]*BT[i][0];
        }
    }
}

void cross(const vector_t<double,3>& A
         , const vector<double>& b
                ,vector_t<double,3>& res)
{
    for (size_t k = 0; k < 3; ++k) {
        for (size_t i = 0; i < 3; ++i) {
            res[k][i][0] = b[2]*A[k][i][1] - b[1]*A[k][i][2];
            res[k][i][1] = b[0]*A[k][i][2] - b[2]*A[k][i][0];
            res[k][i][2] = b[1]*A[k][i][0] - b[0]*A[k][i][1];
        }
    }
}


void rotation_tensor_transpose_diff(      vector_t<double,3>& dLTdv
                                  , const vector_t<double,3>& dLdv)
{
    for (size_t k = 0; k < 3; ++k) {
        for (size_t i = 0; i < 3; i++) {
            for (size_t j = 0; j < 3; ++j) {
                dLTdv[k][j][i] = dLdv[k][i][j];
            }
        }
    }
}


void zhilin_tensor_transpose_diff(        vector_t<double,3>& dBTdv
                                  , const vector_t<double,3>& dBdv)
{
    for (size_t k = 0; k < 3; ++k) {
        for (size_t i = 0; i < 3; i++) {
            for (size_t j = 0; j < 3; ++j) {
                dBTdv[k][j][i] = dBdv[k][i][j];
                // dBTdv[j][k][i] = dBdv[k][i][j];
            }
        }
    }
}

} // namespace math