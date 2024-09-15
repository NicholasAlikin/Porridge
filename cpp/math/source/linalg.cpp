#include "linalg.hpp"

namespace math {

vector<double> solve(const vector_t<double,2>& A, const vector<double>& b) {
    if (b.size() == 1)
        return { b[0]/A[0][0] };

    vector_t<double,2> L;
    vector<double> D(b.size());
    LDLT(A,D,L);
    vector<double> x(b.size());
    gauss_backward_tridown(x,transpose(L),b);
    gauss_backward_triup(x,L,x/D);
    
    return x;
}

vector<double> psolve(const vector_t<double,2>& A, const vector<double>& b) {
    /* A*x = b
       dim(A) = [n,m]
       dim(x) = [m,1]
       dim(b) = [n,1]
       x
    */
   throw 1;
   auto AT = transpose(A);
   return solve(dot(AT,A),dot(AT,b));
}


void LDLT(const vector_t<double,2>& K, vector_t<double>& D, vector_t<double,2>& L) {
    size_t n = K.size();
    L = eye<double>(n);
    vector<double> g = zeros<double>(n);

    D[0] = K[0][0];
    size_t mj = 0;
    for (size_t j = 1; j < n; ++j) {
        g[mj] = K[mj][j];
        for (size_t i = mj+1; i < j; ++i) {
            g[i] = K[i][j];
            for (size_t r = mj; r < i; ++r) {
                g[i] -= L[r][i]*g[r];
            }
        }
        for (size_t i = mj; i < j; ++i) {
            L[i][j] = g[i]/D[i];
        }
        D[j] = K[j][j];
        for (size_t r = mj; r < j; ++r) {
            D[j] -= L[r][j]*g[r];
        }
        
    }
}

/*Matrises like array of columns*/
void LDLT2(const vector_t<double,2>& K, vector_t<double>& D, vector_t<double,2>& L) {
    size_t n = K.size();
    L = eye<double>(n);
    vector<double> g = zeros<double>(n);

    D[0] = K[0][0];
    size_t mj = 0;
    auto k_col = K.begin()+1;//, k_col_end = K.end();
    for (size_t j = 1; j < n; ++j,++k_col) {
        g[mj] = (*k_col)[mj];
        for (size_t i = mj+1; i < j; ++i) {
            g[i] = K[i][j];
            for (size_t r = mj; r < i; ++r) {
                g[i] -= L[r][i]*g[r];
            }
        }
        for (size_t i = mj; i < j; ++i) {
            L[i][j] = g[i]/D[i];
        }
        D[j] = (*k_col)[j];
        for (size_t r = mj; r < j; ++r) {
            D[j] -= L[r][j]*g[r];
        }
        
    }
}

void gauss_backward_triup(vector<double>& x, const vector_t<double,2>& A, const vector<double>& b) {
    size_t n = x.size();
    double sum;
    for (int i = n-1; i > -1; --i) {
        sum = 0;
        for (size_t j = i+1; j < n; ++j) {
            sum += A[i][j]*x[j];
        }
        x[i] = (b[i] - sum) / A[i][i];
    }
}

void gauss_backward_tridown(vector<double>& x, const vector_t<double,2>& A, const vector<double>& b) {
    size_t n = x.size();
    double sum;
    for (size_t i = 0; i < n; ++i) {
        sum = 0;
        for (size_t j = 0; j < i; ++j) {
            sum += A[i][j]*x[j];
        }
        x[i] = (b[i] - sum) / A[i][i];
    }
}

vector_t<double,2> pinv(const vector_t<double,2>& A) {
    // TODO
    return A;
}

double det_(const vector_t<double,2>& mat) {
    vector<double> diag(mat.size());
    vector_t<double,2> L = zeros<double>(mat);
    LDLT(mat, diag,L);
    return product(diag);
}


} // namespace math