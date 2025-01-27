#include "linalg.hpp"

namespace math {

/*===========================
    LLT
===========================*/
void solve_llt(const vector_t<double,2>& A, const vector<double>& b, vector<double>& x, vector_t<double,2>& L, size_t n)
{
    LLT(A,L,n);
    
    // solve Ly=b, L - lower triangular matrix
    x[0] = b[0] / L[0][0];
    double tmp;
    for (size_t i = 1; i < n; ++i) {
        tmp = 0;
        for (size_t j = 0; j < i; ++j) {
            tmp += L[i][j]*x[j];
        }
        x[i] = (b[i] - tmp)/L[i][i];
    }
    // solve xL=y
    x[n-1] /= L[n-1][n-1];
    for (size_t j = n-1; j > 0; --j) {
        tmp = 0;
        for (size_t i = j; i < n; ++i) {
            tmp += L[i][j-1]*x[i];
        }
        x[j-1] = (x[j-1] - tmp)/L[j-1][j-1];
    }
}

void solve_llt(const vector_t<double,2>& A, const vector<double>& b, vector<double>& x, vector_t<double,2>& L)
{
    solve_llt(A,b,x,L,A.size());
}

void solve_llt(const vector_t<double,2>& A, const vector<double>& b, vector<double>& x, size_t n)
{
    vector_t<double,2> L = zeros<double>(A);
    solve_llt(A,b,x,L,n);
}

void solve_llt(const vector_t<double,2>& A, const vector<double>& b, vector<double>& x)
{   
    solve_llt(A,b,x,A.size());
}

void LLT(const vector_t<double,2>& A, vector_t<double,2>& L, size_t n) {
    /*Cholesky decomposition
    using iterators (~fast)*/
    double tmp;

    // iterators
    auto Arowi = A.begin();          // iterator on A row
    decltype(A[0].begin()) Arowi_col;    // iterator on A column
    decltype(L.begin()) Lrowi = L.begin(), Lrowk;
    decltype(L[0].begin()) Lrowi_col, Lrowk_col;


    for (size_t i = 0; i < n; ++i, ++Arowi, ++Lrowi) {
        Arowi_col = Arowi->begin();
        Lrowk = L.begin();
        for (size_t k = 0; k < i; ++k, ++Arowi_col, ++Lrowk) {
            tmp = *Arowi_col;
            Lrowi_col = Lrowi->begin();
            Lrowk_col = Lrowk->begin();
            for (size_t j = 0; j < k; ++j, ++Lrowi_col, ++Lrowk_col) {
                tmp -= (*Lrowi_col) * (*Lrowk_col);
            }
            *Lrowi_col = tmp/(*Lrowk_col);
        }

        tmp = *Arowi_col;
        Lrowi_col = Lrowi->begin();
        for (size_t j = 0; j < i; ++j,++Lrowi_col) {
            tmp -= (*Lrowi_col) * (*Lrowi_col);
        }
        *Lrowi_col = std::sqrt(tmp);
    }
}

void LLT(const vector_t<double,2>& A, vector_t<double,2>& L) {
    /*Cholesky decomposition
    using iterators (~fast)*/
    LLT(A,L,A.size());
}
/*===========================
    LLT END
===========================*/
/*===========================
    LDLT
===========================*/

// void solve_ldlt(const vector_t<double,2>& A, const vector<double>& b, vector<double>& x, vector_t<double,2>& L, vector<double>& D, size_t n) {
//     /* 
//     A * U = b
    
//     A = L * D * LT
//     L * V = b
//     LT * U = D-1 * V
//     */
//     LDLT(A,L,D,x,n);

//     auto itx = x.begin();
    
//     auto itb = b.begin();
    
//     auto itL_row = L.begin();
//     auto itL_row_end = itL_row + n;
//     decltype(L[0].begin()) itL_col, itL_col_end;

//     // solve Ly=b, L - lower triangular matrix
//     // x[0] = b[0];
//     *itx = *itb;
//     ++itb;
//     ++itL_row;
//     double tmp;
//     for (size_t i = 1; i < n; ++i) {
//         tmp = 0.0;
//         itL_col = itL_row->begin();
//         for (size_t j = 0; j < i; ++j) {
//             tmp += (*itL_col) * (*itx);
//             // tmp += L[i][j]*x[j];
//             ++itx;
//             ++itL_col;
//         }
//         // x[i] = b[i] - tmp;
//         *itx = *itb - tmp;
//         ++itb;
//         ++itL_row;
//         itx = x.begin();

//     }
    
//     // x /= D;
//     // auto itx = x.begin();
//     auto itx_end = itx + n;
//     auto itD = D.begin();
//     while (itx < itx_end) {
//         *itx /= *itD;
//         ++itx;
//         ++itD;
//     }
//     // solve xL=y
    
//     for (size_t i = n-1; i > 0; --i) { // x[i=0] will be calculated automatically
//         for (size_t j = 0; j < i; ++j) {
//             x[j] -= L[i][j]*x[i];
//         }
        
//     }
// }

void solve_ldlt(const vector_t<double,2>& A, const vector<double>& b, vector<double>& x, vector_t<double,2>& L, vector<double>& D, size_t n) {
    /* 
    A * U = b
    
    A = L * D * LT
    L * V = b
    LT * U = D-1 * V
    */
    LDLT(A,L,D,x,n);

    // solve Ly=b, L - lower triangular matrix
    x[0] = b[0];
    double tmp;
    for (size_t i = 1; i < n; ++i) {
        tmp = 0.0;
        for (size_t j = 0; j < i; ++j) {
            tmp += L[i][j]*x[j];
        }
        x[i] = b[i] - tmp;
    }
    
    // x /= D;
    auto itx = x.begin();
    auto itx_end = itx + n;
    auto itD = D.begin();
    while (itx < itx_end) {
        *itx /= *itD;
        ++itx;
        ++itD;
    }
    // solve xL=y
    for (size_t i = n-1; i > 0; --i) { // x[i=0] will be calculated automatically
        // tmp = 0.0;
        for (size_t j = 0; j < i; ++j) {
            x[j] -= L[i][j]*x[i];
        }
        // x[j-1] -= tmp;
    }
    // Slice slx(x.begin(),x.begin()+n);
    // std::cout << "#|Ax-b| = " << norm(dot(A,slx)-b) << ", |x| = " << norm(x) << std::endl;

}

void solve_ldlt(const vector_t<double,2>& A, const vector<double>& b, vector<double>& x) {
    size_t n = A.size();
    vector_t<double,2> L = zeros<double>(n,n);
    vector<double> D(n);
    solve_ldlt(A,b,x,L,D,n);
}

void solve_ldlt(const vector<double>& A, const vector<size_t>& diag, const vector<double>& b, vector<double>& x, vector<double>& LT, vector<double>& D, size_t n) {
    /* 
    A * U = b
    
    A = L * D * LT
    L * V = b
    LT * U = D-1 * V
    */
    LDLT(A,diag,LT,D,x,n);
    
    // size_t n = diag.size()-1;
    
    size_t mj;
    
    // upper triangular
    x[0] = b[0];
    double tmp;
    for (size_t j = 1; j < n; ++j) {
        tmp = 0;
        mj = j+1 - (diag[j+1]-diag[j]);
        for (size_t i = mj; i < j; ++i) {
            tmp += LT[diag[j]+j-i]*x[i];
        }
        x[j] = b[j] - tmp;
    }
    
    auto itx = x.begin();
    auto itx_end = itx + n;
    auto itD = D.begin();
    while (itx < itx_end) {
        *itx /= *itD;
        ++itx;
        ++itD;
    }
    // x /= D;

    // lower triangular
    // x.last() = b.last();
    // loop over columns
    for (size_t j = n-1; j > 0; --j) {
        mj = j+1 - (diag[j+1]-diag[j]);
        // std::cout << "j = " << j << ", mj = " << mj << '\n';
        // loop over rows of j-th column
        for (size_t i = mj; i < j; ++i) {
            // std::cout << "\ti = " << i << ", j = " << j <<", &L = " << diag[j]+j-i << ", Lij = " << L[diag[j]+j-i] << std::endl;
            x[i] -= LT[diag[j]+j-i]*x[j];
        }
        // x[j-1] = x[j-1];
    }
    

    // return x;
}

void solve_ldlt(const vector_t<double,1>& A, const vector<size_t>& diag, const vector<double>& b, vector<double>& x) {
    size_t n = diag.size()-1;
    vector<double> LT(A.size()), D(n);
    solve_ldlt(A,diag,b,x,LT,D,n);
}


void LDLT(const vector_t<double,2>& A, vector_t<double,2>& L, vector<double>& D, vector<double>& g, size_t n) {
    /*
    LT - upper trianguar unit matrix
    */
    // size_t n = A.size();
    
    decltype(D.begin()) itD_i, itD_j;
    
    auto itA_row = A.begin()+1;
    decltype(A[0].begin()) itA_row_col;

    decltype(g.begin()) itg_j, itg_r;
    double temp;

    decltype(L.begin()) itL_row;
    decltype(L[0].begin()) itL_row_col;

    
    D[0] = A[0][0];
    for (size_t i = 1; i < n; ++i) {
        itA_row_col = itA_row->begin();
        itg_j = g.begin();
        *itg_j = *itA_row_col;

        ++itA_row_col;
        ++itg_j;
        itL_row = L.begin()+1;
        for (size_t j = 1; j < i; ++j) {
            temp = *itA_row_col;
            itg_r = g.begin();
            itL_row_col = itL_row->begin();
            for (size_t r = 0; r < j; ++r) {
                temp -= (*itL_row_col) * (*itg_r);
                ++itg_r;
                ++itL_row_col;
            }
            *itg_j = temp;

            ++itA_row_col;
            ++itg_j;
            ++itL_row;
        }
        
        itL_row_col = itL_row->begin();
        itg_j = g.begin();
        itD_j = D.begin();
        for (size_t j = 0; j < i; ++j) {
            *itL_row_col = (*itg_j) / (*itD_j);
            ++itL_row_col;
            ++itg_j;
            ++itD_j;
        }

        temp = *itA_row_col;
        itL_row_col = itL_row->begin();
        itg_j = g.begin();
        for (size_t r = 0; r < i; ++r) {
            temp -= (*itL_row_col) * (*itg_j);
            ++itL_row_col;
            ++itg_j;
        }
        *itD_j = temp;

        ++itA_row;   
    }
}
/* // Legacy
void LDLT(const vector_t<double,2>& A, vector_t<double,2>& L, vector<double>& D, vector<double>& g, size_t n) {
    
    D[0] = A[0][0];
    for (size_t j = 1; j < n; ++j) {
        
        g[0] = A[0][j];
        for (size_t i = 0+1; i < j; ++i) {
            g[i] = A[i][j];
            for (size_t r = 0; r < i; ++r) {
                g[i] -= L[i][r]*g[r];
            }
        }
        
        for (size_t i = 0; i < j; ++i) {
            L[j][i] = g[i]/D[i];
        }
        D[j] = A[j][j];
        for (size_t r = 0; r < j; ++r) {
            D[j] -= L[j][r]*g[r];
        }
        
    }
}
*/

void LDLT(const vector_t<double,2>& A, vector_t<double,2>& L, vector<double>& D) {
    size_t n = A.size();
    vector<double> g(n);
    LDLT(A,L,D,g,n);
}

void LDLT(const vector<double>& A, const vector<size_t>& diag, vector<double>& LT, vector<double>& D, vector<double>& g, size_t n) {
    /* 
    A * U = b
    
    A = L * D * LT
    L * V = b
    LT * U = D-1 * V
    */
    // size_t n = diag.size()-1;
    
    D[0] = A[0];
    size_t mj;
    size_t mi;

    // loop over columns
    for (size_t j = 1; j < n; ++j) {
        mj = j+1 - (diag[j+1]-diag[j]);
        
        g[mj] = A[diag[j]+j-mj]; // A[mj][j];
        for (size_t i = mj+1; i < j; ++i) {
            g[i] = A[diag[j]+j-i];
            mi = i+1 - (diag[i+1]-diag[i]);
            for (size_t r = std::max(mj,mi); r < i; ++r) {
                g[i] -= LT[diag[i]+i-r]*g[r];
            }
        }
        for (size_t i = mj; i < j; ++i) {
            LT[diag[j]+j-i] = g[i]/D[i];
        }
        D[j] = A[diag[j]];
        for (size_t r = mj; r < j; ++r) {
            D[j] -= LT[diag[j]+j-r]*g[r];
        }
        
    }
}

void LDLT(const vector_t<double,1>& A, const vector<size_t>& diag, vector<double>& LT, vector<double>& D) {
    size_t n = A.size();
    vector<double> g(n);
    LDLT(A,diag,LT,D,g,n);
}


/*===========================
    LDLT END
===========================*/

/*===========================
    LU
===========================*/


void solve_lu(vector_t<double,2>& A, const vector<double>& b, vector<double>& x, size_t n) {
    /*Solve Ax=b for x using LU decomposition,
    where A is a sqare matrix of the general form
    (not necessary symmetric).
    */
    // LU decomposition
    LU(A,n);


    // solve Ly=b for y with lower triangular matrix L
    // vector<double> y(n);
    x[0] = b[0];
    double tmp_sum;
    size_t i,j;
    for (i = 1; i < n; ++i) {
        tmp_sum = 0;
        for (j = 0; j < i; ++j) {
            tmp_sum += x[j]*A[i][j];
        }
        x[i] = b[i] - tmp_sum;
    }

    // solve Ux=y for x with upper triangular matrix U
    x[n-1] /= A[n-1][n-1];
    for (i = n-1; i > 0; --i) {
        tmp_sum = 0;
        for (j = i; j < n; ++j) {
            tmp_sum += x[j]*A[i-1][j];
        }
        x[i-1] = (x[i-1] - tmp_sum) / A[i-1][i-1];
    }
}

void solve_lu(vector_t<double,2>& A, const vector<double>& b, vector<double>& x) {
    solve_lu(A,b,x,A.size());
}

void LU(vector_t<double,2>& A,  /*vector<size_t>& pos,*/ size_t n) {
    /*LU decomposition
    L is the lower triangular matrix with 1 on the diagonal
    U is the upper triangular matrix
    L and U are placed into A matrix as well
    
    pos - containts original matrix row positions in decomposed matrix,
        neaded for solve LU*x=b
    */
    
    // size_t k = 0;
    // for (auto posk: pos) {
    //     posk = k;
    //     ++k;
    // }

    for (size_t p = 1; p <= n-1; ++p) {
        /* Find pivot element */




        for (size_t i = p+1; i <= n; ++i) {
            // calculate L elements of k-th column
            A[i-1][p-1] /= A[p-1][p-1];
            
            // loop over i-th row upper triangular part of A - U 
            // lower triangular part mast be equal to 0
            for (size_t j = p+1; j <= n; ++j) {
                A[i-1][j-1] -= A[i-1][p-1] * A[p-1][j-1];
            }
        }
    }
}

void LU(vector_t<double,2>& A) {
    LU(A,A.size());
}


/*===========================
    LU END
===========================*/

void psolve(const vector_t<double,2>& A, const vector<double>& b, vector<double>& x, vector_t<double,2>& AAT
            , vector_t<double,2>& L, vector<double>& D) {
    /* A*x = b
       dim(A) = [n,m]
       dim(x) = [m,1]
       dim(b) = [n,1]
       x
    */
    dotT(A,A,AAT);
    solve_ldlt(AAT,b,D,L,x,AAT.size());
    fill(x.begin(),x.end(),0.0);
    dot(D, A, x);

    // solve_ldlt(dotT(A,A),b,D);
    // x = dot(D,A);
    
}

void psolve(const vector_t<double,2>& A, const vector<double>& b, vector<double>& x) {
    /* Solve A*x = b for x: |x| -> min
       dim(A) = [n,m]
       dim(x) = [m,1]
       dim(b) = [n,1] */
    solve_llt(dotT(A,A),b,x); // dotT(A,A) = [A.size() x A.size()]
    Slice slx(x.begin(),x.end()-1);
    x = dot(slx, A); 
}


void psolve_weight(const vector_t<double,2>& A, const vector<double>& b, const vector_t<double,2>& invW, vector<double>& x) {
    /* Solve: A*x = b for x: x*W*x -> min
       dim(A) = [n,m]
       dim(x) = [m,1]
       dim(b) = [n,1]
       invW = W^-1    
       
       x = invW * A^T * (A*invD*A^T)^-1 * b

       1) calculate A*invD*A^T
       2) solve (A*invD*A^T)*y = b  for y
       3) calculate y*A    
       4) calculate invW * (y*A)
    */
    vector<double> y(b.size());
    solve_ldlt(dotT(dot(A,invW),A),b,y);
    dot(invW,dot(y,A),x);
    
    
}

// double det_(const vector_t<double,2>& mat) {
//     vector_t<double,2> A = mat;
//     LU(A);
//     double res = 0;
//     size_t n = A.size();
//     for (size_t i = 0; i < n; ++i) {
//         res += A[i][i];
//     }
//     return res;
// }

void invLowTri(const vector_t<double,2>& A, vector_t<double,2>& invA) {
    size_t n = A.size();
    
    for (size_t i = 0; i < n; ++i) {
        invA[i][i] = 1.0/A[i][i];
        for (size_t j = 0; j < i; ++j) {
            double s = 0.0;
            for (size_t k = j; k < i; ++k) {
                s += A[i][k]*invA[k][j];
            }
            invA[i][j] = -s*invA[i][i];
        }
    }
}

void invLowTriUnit(const vector_t<double,2>& A, vector_t<double,2>& invA) {
    size_t n = A.size();

    for (size_t i = 0; i < n; ++i) {
        for (size_t j = 0; j < i; ++j) {
            double s = 0.0;
            for (size_t k = j; k < i; ++k) {
                s += A[i][k]*invA[k][j];
            }
            invA[i][j] = -s;
        }
    }
}

void invUpTri(const vector_t<double,2>& A, vector_t<double,2>& invA) {
    size_t n = A.size();
    
    for (int i = n-1; i > -1; --i) {
        invA[i][i] = 1.0/A[i][i];
        for (int j = n-1; j > i; --j) {
            double s = 0.0;
            for (int k = j; k > i; --k) {
                s += A[i][k]*invA[k][j];
            }
            invA[i][j] = -s*invA[i][i];
        }
    }
}

// void inv(const vector_t<double,2>& A, vector_t<double,2>& invA) {
//     vector_t<double,2> A_ = A;
//     LU(A_);
//     vector_t<double,2> invA = dot(invUpTri(A_),invLowTriUnit(A_));
//     return invA;
// }

void invSym(const vector_t<double,2>& A, vector_t<double,2>& invA
            , vector_t<double,2>& L, vector<double>& D
            , vector<double>& temp, vector_t<double,2>& invL, size_t n)
{
    /*Inverse symmetric matrix*/
    LDLT(A,L,D,temp,n);
    invLowTriUnit(L,invL);
    
    // calculate  D*LT
    /*
    auto L_row = L.begin(),
         L_row_end = L.end();
    decltype(L[0].begin()) L_col, L_col_end;
    auto Di = D.begin();
    double Di_val;
    size_t row = 0;
    while (L_row != L_row_end) {
        L_col = L_row->begin() + row;
        L_col_end = L_row->end();
        Di_val = *Di;
        while (L_col != L_col_end) {
            *L_col *= Di_val;
            ++L_col;
        }
        
        ++L_row;
        ++Di;
        ++row;
    }*/
    // calculate D*LT
    vector_t<double,2> L_ = L;
    for (size_t i = 0; i < n; ++i) {
        // loop over lover triangular matrix
        L[i][i] = D[i];
        for (size_t j = 0; j < i; ++j) {
            L[j][i] = L[i][j]*D[j];
            L[i][j] = 0.0;
        }
    }
    invUpTri(L,invA);
    
    invA = dotUL(invA,invL);
}

void invSym(const vector_t<double,2>& A, vector_t<double,2>& invA) {
    size_t n = A.size();
    vector_t<double,2> L = eye<double>(n);
    vector<double> D(n);
    vector<double> temp(n);
    vector_t<double,2> invL(L);
    invSym(A,invA,L,D,temp,invL,n);
}

// vector_t<double,2> pinv(const vector_t<double,2>& A) {
//     /*A+ = AT * (A,AT)^-1 */
//     auto AT = transpose(A);
//     vector_t<double,2> pinvA = dot(std::move(AT), inv(dot(A,AT)));
//     return pinvA;
// }


/* Cross product of two tensors of the second rank */
void cross(const vector_t<double,2>& A
          ,const vector_t<double,2>& B
               , vector_t<double,3>& res)
{
    // loop over 3rd dimension
    for (size_t k = 0; k < 3; ++k) {
        for (size_t i = 0; i < 3; ++i) {
            // first row
            res[k][i][0] = A[k][1]*B[2][i] - A[k][2]*B[1][i];
            // second row
            res[k][i][1] = A[k][2]*B[0][i] - A[k][0]*B[2][i];
            // third row
            res[k][i][2] = A[k][0]*B[1][i] - A[k][1]*B[0][i];
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
            // first row
            res[k][i][0] = AT[1][k]*B[2][i] - AT[2][k]*B[1][i];
            // second row
            res[k][i][1] = AT[2][k]*B[0][i] - AT[0][k]*B[2][i];
            // third row
            res[k][i][2] = AT[0][k]*B[1][i] - AT[1][k]*B[0][i];
        }
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
            // first row
            res[k][i][0] = A[k][1]*BT[i][2] - A[k][2]*BT[i][1];
            // second row
            res[k][i][1] = A[k][2]*BT[i][0] - A[k][0]*BT[i][2];
            // third row
            res[k][i][2] = A[k][0]*BT[i][1] - A[k][1]*BT[i][0];
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



void rotation_tensor_diff(const vector<double>& v
                        , vector_t<double,3>& dLdv
                        , const vector_t<double,2>& L
                        , const vector_t<double,2>& B)
{
    crossT0(B,L,dLdv);
}

void rotation_tensor_transpose_diff(const vector<double>& v
                                  , vector_t<double,3>& dLTdv
                                  , const vector_t<double,2>& L
                                  , const vector_t<double,2>& B)
{
    // -cross0T(B,L,dLTdv);
    // loop over 3rd dimension
    for (size_t k = 0; k < 3; ++k) {
        for (size_t i = 0; i < 3; ++i) {
            // first row
            dLTdv[k][i][0] = B[k][2]*L[i][1] - B[k][1]*L[i][2];
            // second row
            dLTdv[k][i][1] = B[k][0]*L[i][2] - B[k][2]*L[i][0];
            // third row
            dLTdv[k][i][2] = B[k][1]*L[i][0] - B[k][0]*L[i][1];
        }
    }
    
}


void zhilin_tensor_diff(const vector<double>& v
                      , vector_t<double,3>& dBdv
                      , const vector_t<double,2>& B
                      , const vector_t<double,3>& dLdv
                      , double eps)
{
    double norm_v = norm(v);
    if (norm_v < eps) {
        dBdv = {
            {{ 0.0,  0.0,  0.0},
             { 0.0,  0.0,  0.5},
             { 0.0, -0.5,  0.0}},

            {{ 0.0,  0.0, -0.5},
             { 0.0,  0.0,  0.0},
             { 0.5,  0.0,  0.0}},

            {{ 0.0,  0.5,  0.0},
             {-0.5,  0.0,  0.0},
             { 0.0,  0.0,  0.0}}
        };
        return;
    }

    vector<double> v_ = v/(-norm_v*norm_v); // e/|v|
    
    for (size_t k = 0; k < 3; ++k) {
        for (size_t i = 0; i < 3; ++i) {
            // v_ o v o v_ + E o v_ - BT o v_ - v_ o B - Z x v_
            dBdv[k][i][0] = (v_[k]*v[i]*v_[0]) + ((1.0 ? k==i : 0.0)*v_[0]) - (B[i][k]*v_[0]) - (v_[k]*B[i][0]) - (v_[2]*dLdv[k][i][1] - v_[1]*dLdv[k][i][2]);
            dBdv[k][i][1] = (v_[k]*v[i]*v_[1]) + ((1.0 ? k==i : 0.0)*v_[1]) - (B[i][k]*v_[1]) - (v_[k]*B[i][1]) - (v_[0]*dLdv[k][i][2] - v_[2]*dLdv[k][i][0]);
            dBdv[k][i][2] = (v_[k]*v[i]*v_[2]) + ((1.0 ? k==i : 0.0)*v_[2]) - (B[i][k]*v_[2]) - (v_[k]*B[i][2]) - (v_[1]*dLdv[k][i][0] - v_[0]*dLdv[k][i][1]);
        }
    }

}

void zhilin_tensor_transpose_diff(const vector<double>& v
                                , vector_t<double,3>& dBTdv
                                , const vector_t<double,2>& B
                                , const vector_t<double,3>& dLTdv
                                , double eps)
{
    double norm_v = norm(v);
    if (norm_v < eps) {
        dBTdv = {
            {{ 0.0,  0.0,  0.0},
             { 0.0,  0.0, -0.5},
             { 0.0,  0.5,  0.0}},

            {{ 0.0,  0.0,  0.5},
             { 0.0,  0.0,  0.0},
             {-0.5,  0.0,  0.0}},

            {{ 0.0, -0.5,  0.0},
             { 0.5,  0.0,  0.0},
             { 0.0,  0.0,  0.0}}
        };
        return;
    }

    vector<double> v_ = v/(-norm_v*norm_v); // e/|v|
    
    for (size_t k = 0; k < 3; ++k) {
        for (size_t i = 0; i < 3; ++i) {
            // v_ o v o v_ + E o v_ - BT o v_ - v_ o B - Z x v_                  
            dBTdv[k][i][0] = (v_[k]*v[i]*v_[0]) + ((1.0 ? k==i : 0.0)*v_[0]) - (B[k][i]*v_[0]) - (v_[k]*B[0][i]) + (v_[2]*dLTdv[k][i][1] - v_[1]*dLTdv[k][i][2]);
            dBTdv[k][i][1] = (v_[k]*v[i]*v_[1]) + ((1.0 ? k==i : 0.0)*v_[1]) - (B[k][i]*v_[1]) - (v_[k]*B[1][i]) + (v_[0]*dLTdv[k][i][2] - v_[2]*dLTdv[k][i][0]);
            dBTdv[k][i][2] = (v_[k]*v[i]*v_[2]) + ((1.0 ? k==i : 0.0)*v_[2]) - (B[k][i]*v_[2]) - (v_[k]*B[2][i]) + (v_[1]*dLTdv[k][i][0] - v_[0]*dLTdv[k][i][1]);
        }
    }

}

} // namespace math