#pragma once


/* Tensor operations.hpp 
    Constaints tensor operation.
*/


#include "linalg.hpp"


namespace math {

template <Vector V1, Vector V2>
auto dyad(const V1& v1, const V2& v2)
		-> vector_t<general_type_t<V1,V2>,2>
{
	vector_t<general_type_t<V1,V2>,2> res = {
		v1[0]*v2,
		v1[1]*v2,
		v1[2]*v2
	};
	return res;
}

template <Vector V>
auto dyad(const V& v1, const V& v2, const V& v3)
		-> vector_t<typename V::basic_value_type,3>
{
	vector_t<typename V::basic_value_type,3> res = zeros<typename V::basic_value_type>(v1.size(),v2.size(),v3.size());
    for (size_t k = 0; k < v1.size(); ++k) {
        for (size_t i = 0; i < v2.size(); ++i) {
            double v1k_v2i = v1[k]*v2[i];
            for (size_t j = 0; j < v3.size(); ++j) {
                res[k][i][j] = v1k_v2i*v3[j];
            }
        }
    }
    return res;
}

template <Matrix M, Vector V>
requires HaveGeneralType<M,V>
auto dyad(const M& mat, const V& vec)
		-> vector_t<general_type_t<M,V>,3>
{
    auto mat_sz = size(mat);
    vector_t<general_type_t<M,V>,3> res = zeros<general_type_t<M,V>>(mat_sz[0],mat_sz[1],vec.size());
    for (size_t k = 0; k < mat_sz[0]; ++k) {
        for (size_t i = 0; i < mat_sz[1]; ++i) {
            for (size_t j = 0; j < vec.size(); ++j) {
                res[k][i][j] = mat[k][i]*vec[j];
            }
        }
    }
    return res;
}

template <Matrix M, Vector V>
requires HaveGeneralType<M,V>
auto dyadT0(const M& matT, const V& vec)
		-> vector_t<general_type_t<M,V>,3>
{
    auto mat_sz = size(matT);
    vector_t<general_type_t<M,V>,3> res = zeros<general_type_t<M,V>>(mat_sz[1],mat_sz[0],vec.size());
    for (size_t k = 0; k < mat_sz[0]; ++k) {
        for (size_t i = 0; i < mat_sz[1]; ++i) {
            for (size_t j = 0; j < vec.size(); ++j) {
                res[k][i][j] = matT[i][k]*vec[j];
            }
        }
    }
    return res;
}

template <Matrix M, Vector V>
requires HaveGeneralType<M,V>
auto dyad(const V& vec, const M& mat)
		-> vector_t<general_type_t<M,V>,3>
{
    auto mat_sz = size(mat);
    vector_t<general_type_t<M,V>,3> res = zeros<general_type_t<M,V>>(vec.size(),mat_sz[0],mat_sz[1]);
    for (size_t k = 0; k < vec.size(); ++k) {
        for (size_t i = 0; i < mat_sz[0]; ++i) {
            for (size_t j = 0; j < mat_sz[1]; ++j) {
                res[k][i][j] = vec[k]*mat[i][j];
            }
        }
    }
    return res;
}

template <Matrix M, Vector V>
requires HaveGeneralType<M,V>
auto dyad0T(const V& vec, const M& matT)
		-> vector_t<general_type_t<M,V>,3>
{
    auto mat_sz = size(matT);
    vector_t<general_type_t<M,V>,3> res = zeros<general_type_t<M,V>>(vec.size(),mat_sz[1],mat_sz[0]);
    for (size_t k = 0; k < vec.size(); ++k) {
        for (size_t i = 0; i < mat_sz[0]; ++i) {
            for (size_t j = 0; j < mat_sz[1]; ++j) {
                res[k][i][j] = vec[k]*matT[j][i];
            }
        }
    }
    return res;
}


namespace detail {
template <NotVectorLike T>
double rotation_tensor_helper1(const T& x, double eps = 1e-9) {
	return (x < eps) ? 0.5 : (1. - std::cos(x))/(x*x);
}
template <NotVectorLike T>
double rotation_tensor_helper2(const T& x, double eps = 1e-9) {
	return (x < eps) ? 1.0 : std::sin(x)/x;
}
template <NotVectorLike T>
double rotation_tensor_helper3(const T& x, double eps = 1e-9) {
	return (x < eps) ? 1.0/6.0 : (x - std::sin(x))/(x*x*x);
}
} // namespace detail

template <Vector V>
void rotation_tensor(const V& v, vector_t<double,2>& L, double eps = 1e-12) {
	/* L = E*cos(|v|) + vv*f1(|v|) + spin(v)*f2(|v|) */
	
	double abs_v = norm(v);
	double cos_v = cos(abs_v);
    double f1 = detail::rotation_tensor_helper1(abs_v,eps);
	double f2 = detail::rotation_tensor_helper2(abs_v,eps);
    double v_f1[] = {v[0]*f1, v[1]*f1, v[2]*f1};
    double v_f2[] = {v[0]*f2, v[1]*f2, v[2]*f2};

	
	L[0][0] = cos_v + v[0]*v_f1[0];
    L[0][1] =         v[0]*v_f1[1] - v_f2[2];
    L[0][2] =       + v[0]*v_f1[2] + v_f2[1];

    L[1][0] =         v[1]*v_f1[0] + v_f2[2];
    L[1][1] = cos_v + v[1]*v_f1[1];
    L[1][2] =       + v[1]*v_f1[2] - v_f2[0];

    L[2][0] =         v[2]*v_f1[0] - v_f2[1];
    L[2][1] =       + v[2]*v_f1[1] + v_f2[0];
    L[2][2] = cos_v + v[2]*v_f1[2];
    
}

template <Vector V>
vector_t<double,2> rotation_tensor(const V& v, double eps = 1e-12) {
	/* L = E*cos(|v|) + vv*f1(|v|) + spin(v)*f2(|v|) */
	
	
	vector_t<double,2> L = zeros<double>(v.size(),v.size());
    rotation_tensor(v,L,eps);
    return L;
}


template <Vector V>
void zhilin_tensor(const V& v, vector_t<double,2>& B, double eps = 1e-12) {
	/* B = E*sin(|v|)/|v| + vv*(|v| - sin(|v|))/|v|^3 + skew(v)*(1 - cos(|v|))/|v|^2 */
	
	double abs_v = norm(v);
	double f1 = detail::rotation_tensor_helper1(abs_v,eps);
	double f2 = detail::rotation_tensor_helper2(abs_v,eps);
    double f3 = detail::rotation_tensor_helper3(abs_v,eps);
    double v_f1[] = {v[0]*f1, v[1]*f1, v[2]*f1};
    double v_f3[] = {v[0]*f3, v[1]*f3, v[2]*f3};

    B[0][0] = f2    + v[0]*v_f3[0];
    B[0][1] =         v[0]*v_f3[1] - v_f1[2];
    B[0][2] =       + v[0]*v_f3[2] + v_f1[1];

    B[1][0] =         v[1]*v_f3[0] + v_f1[2];
    B[1][1] = f2    + v[1]*v_f3[1];
    B[1][2] =       + v[1]*v_f3[2] - v_f1[0];

    B[2][0] =         v[2]*v_f3[0] - v_f1[1];
    B[2][1] =       + v[2]*v_f3[1] + v_f1[0];
    B[2][2] = f2    + v[2]*v_f3[2];
    
}

template <Vector V>
vector_t<double,2> zhilin_tensor(const V& v, double eps = 1e-12) {
	/* L = E*cos(|v|) + vv*f1(|v|) + spin(v)*f2(|v|) */
	
	
	vector_t<double,2> B = zeros<double>(v.size(),v.size());
    zhilin_tensor(v,B,eps);
    return B;
}

/* Vector invariant of the dot product of 2 tensors (square matrises, size=3):
 mat1 and transpose(mat2) */
template <Matrix M1, ArithmeticVectorsLike<M1> M2, Vector V>
void vector_invariant(const M1& mat1, const M2& mat2T, V& res) {
    res[0] = dot(mat1[1],mat2T[2]) - dot(mat1[2],mat2T[1]);
    res[1] = dot(mat1[2],mat2T[0]) - dot(mat1[0],mat2T[2]);
    res[2] = dot(mat1[0],mat2T[1]) - dot(mat1[1],mat2T[0]);
}

template <Matrix M1, ArithmeticVectorsLike<M1> M2>
auto vector_invariant(const M1& mat1, const M2& mat2T)
        -> vector<general_type_t<M1,M2>>
{   
    vector<general_type_t<M1,M2>> res = {
		math::dot(mat1[1],mat2T[2]) - math::dot(mat1[2],mat2T[1]),
        math::dot(mat1[2],mat2T[0]) - math::dot(mat1[0],mat2T[2]),
        math::dot(mat1[0],mat2T[1]) - math::dot(mat1[1],mat2T[0]),
	};
	return res;
}

template <Matrix M, Vector V>
void vector_invariant(const M& mat, V& res) {    
    res[0] = mat[1][2] - mat[2][1];
    res[1] = mat[2][0] - mat[0][2];
    res[2] = mat[0][1] - mat[1][0];
}

template <Matrix M>
auto vector_invariant(const M& mat)
        -> vector<typename M::basic_value_type>
{    
    vector<typename M::basic_value_type> res = {
		mat[1][2] - mat[2][1],
        mat[2][0] - mat[0][2],
        mat[0][1] - mat[1][0]
	};
	return res;
}


template <Vector V>
auto skew_tensor(const V& vec)
		-> vector_t<typename V::basic_value_type,2>
{
	vector_t<typename V::basic_value_type,2> res = {
		{   0   ,-vec[2], vec[1]},
		{ vec[2],   0   ,-vec[0]},
		{-vec[1], vec[0],   0   }
	};
	return res;
}

/* Cross product of two tensors of the second rank */
void cross(const vector_t<double,2>& A
          ,const vector_t<double,2>& B
               , vector_t<double,3>& res);
/* Cross product of transpose and ordinary tensors of the second rank */
void crossT0(const vector_t<double,2>& AT
            ,const vector_t<double,2>& B
                 , vector_t<double,3>& res);
/* Cross product of transpose tensors of the second rank and vector  */
void crossT0(const vector_t<double,2>& AT
            ,const vector_t<double,1>& b
            , vector_t<double,2>& res);
/* Cross product of ordinary and transpose tensors of the second rank */
void cross0T(const vector_t<double,2>& A
            ,const vector_t<double,2>& BT
                 , vector_t<double,3>& res);
void cross(const vector_t<double,3>& A
         , const vector<double>& b
                ,vector_t<double,3>& res);

template <Vector V, Matrix M, Matrix M1>
void cross(const V& v, const M& m, M1& res) {
    auto v1 = v[0]
        ,v2 = v[1]
        ,v3 = v[2];
    auto m1 = m[0].begin(), m1_end = m[0].end()
        ,m2 = m[1].begin()
        ,m3 = m[2].begin();
    auto res1 = res[0].begin()
        ,res2 = res[1].begin()
        ,res3 = res[2].begin();
    
    while (m1 < m1_end) {
        *res1 = v2* (*m3) - v3* (*m2);
        *res2 = v3* (*m1) - v1* (*m3);
        *res3 = v1* (*m2) - v2* (*m1);
        ++res1; ++res2; ++res3;
        ++m1;   ++m2;   ++m3;
    }
}



template <Vector V>
void rotation_tensor_diff(const V& v
                        , vector_t<double,3>& dLdv
                        , const vector_t<double,2>& L
                        , const vector_t<double,2>& B)
{
    crossT0(B,L,dLdv);
}


template <Vector V>
void rotation_tensor_transpose_diff(V& v
                                  , vector_t<double,3>& dLTdv
                                  , const vector_t<double,2>& L
                                  , const vector_t<double,2>& B)
{
    // -cross0T(B,L,dLTdv);
    // loop over 3rd dimension
    for (size_t k = 0; k < 3; ++k) {
        for (size_t i = 0; i < 3; ++i) {
            dLTdv[k][0][i] = B[k][2]*L[i][1] - B[k][1]*L[i][2];
            dLTdv[k][1][i] = B[k][0]*L[i][2] - B[k][2]*L[i][0];
            dLTdv[k][2][i] = B[k][1]*L[i][0] - B[k][0]*L[i][1];
        }
    }
}



void rotation_tensor_transpose_diff(      vector_t<double,3>& dLTdv
                                  , const vector_t<double,3>& dLdv);


template <Vector V>
void zhilin_tensor_diff(const V& v
                      , vector_t<double,3>& dBdv
                      , const vector_t<double,2>& B
                      , const vector_t<double,3>& dLdv
                      , double eps = 1e-12)
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

    vector<double> v_ = v/(norm_v*norm_v);
    
    for (size_t k = 0; k < 3; ++k) {
        for (size_t i = 0; i < 3; ++i) {
            // - ( sum(en o v_ o en, n=1...3) + E o v_ - BT o v_ - v_ o B - dLdv x v_)
            dBdv[k][i][0] = (1.0 ? k==0 : 0.0)*v_[i] + ((1.0 ? k==i : 0.0)*v_[0]) - (B[i][k]*v_[0]) - (v_[k]*B[i][0]) - (v_[2]*dLdv[k][i][1] - v_[1]*dLdv[k][i][2]);
            dBdv[k][i][1] = (1.0 ? k==1 : 0.0)*v_[i] + ((1.0 ? k==i : 0.0)*v_[1]) - (B[i][k]*v_[1]) - (v_[k]*B[i][1]) - (v_[0]*dLdv[k][i][2] - v_[2]*dLdv[k][i][0]);
            dBdv[k][i][2] = (1.0 ? k==2 : 0.0)*v_[i] + ((1.0 ? k==i : 0.0)*v_[2]) - (B[i][k]*v_[2]) - (v_[k]*B[i][2]) - (v_[1]*dLdv[k][i][0] - v_[0]*dLdv[k][i][1]);
        }
    }
}


template <Vector V>
void zhilin_tensor_transpose_diff(const V& v
                                , vector_t<double,3>& dBTdv
                                , const vector_t<double,2>& B
                                , const vector_t<double,3>& dLTdv
                                , double eps = 1e-12) 
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

    vector<double> v_ = v/(norm_v*norm_v);
    
    for (size_t k = 0; k < 3; ++k) {
        for (size_t i = 0; i < 3; ++i) {
            // - ( sum(en o v_ o en, n=1...3) + E o v_ - B o v_ - v_ o BT + dLTdv x v_)
            // dBTdv[k][0][i] = (1.0 ? k==0 : 0.0)*v_[i] + ((1.0 ? k==i : 0.0)*v_[0]) - (B[i][k]*v_[0]) - (v_[k]*B[i][0]) - (v_[2]*dLdv[k][i][1] - v_[1]*dLdv[k][i][2]);
            // dBTdv[k][1][i] = (1.0 ? k==1 : 0.0)*v_[i] + ((1.0 ? k==i : 0.0)*v_[1]) - (B[i][k]*v_[1]) - (v_[k]*B[i][1]) - (v_[0]*dLdv[k][i][2] - v_[2]*dLdv[k][i][0]);
            // dBTdv[k][2][i] = (1.0 ? k==2 : 0.0)*v_[i] + ((1.0 ? k==i : 0.0)*v_[2]) - (B[i][k]*v_[2]) - (v_[k]*B[i][2]) - (v_[1]*dLdv[k][i][0] - v_[0]*dLdv[k][i][1]);

            dBTdv[k][i][0] = (1.0 ? k==0 : 0.0)*v_[i] + ((1.0 ? k==i : 0.0)*v_[0]) - (B[k][i]*v_[0]) - (v_[k]*B[0][i]) + (v_[2]*dLTdv[k][i][1] - v_[1]*dLTdv[k][i][2]);
            dBTdv[k][i][1] = (1.0 ? k==1 : 0.0)*v_[i] + ((1.0 ? k==i : 0.0)*v_[1]) - (B[k][i]*v_[1]) - (v_[k]*B[1][i]) + (v_[0]*dLTdv[k][i][2] - v_[2]*dLTdv[k][i][0]);
            dBTdv[k][i][2] = (1.0 ? k==2 : 0.0)*v_[i] + ((1.0 ? k==i : 0.0)*v_[2]) - (B[k][i]*v_[2]) - (v_[k]*B[2][i]) + (v_[1]*dLTdv[k][i][0] - v_[0]*dLTdv[k][i][1]);

        }
    }
}


void zhilin_tensor_transpose_diff(        vector_t<double,3>& dBTdv
                                  , const vector_t<double,3>& dBdv);


} // namespace math