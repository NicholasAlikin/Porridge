#include "harmonic_balance.hpp"

namespace math {


HBM::HBM(const Matrix_t& mass
       , const Matrix_t& damp
       , const Matrix_t& stif
       , const DFT& dft
       , const std::function<Vector_t(const Vector_t&,double,const DFT&)>& funcnl
       , const std::function<Vector_t(double,const DFT&)>& funcex)
    : Corrector(mass, damp, stif, dft)
    , Mass(kron(dft.derivative2_basic,mass))
    , Damp(kron(dft.derivative_basic,damp))
    , Stif(kron(eye<double>(2*dft.H+1),stif))
    , funcnl(funcnl), funcex(funcex)
{

}

void HBM::correction_system_response(Vector_t& fun, Matrix_t& jac, const Vector_t& y) {

    Vector_t x(y.begin(),y.end()-1);
    double freq = *(y.end()-1);
    Matrix_t L = linear_system_dynamic_reaction(freq);
    Vector_t fnl = funcnl(x,freq,dft);
    Vector_t fex = funcex(freq,dft);
    fun = correction_system(x,L,fnl,fex);
    jac = corretion_system_jac(x,freq,L,fnl);
}

Matrix_t HBM::linear_system_dynamic_reaction(double freq) {
    Matrix_t L = (freq*freq)*Mass + freq*Damp + Stif;
    return L;
}

Vector_t HBM::correction_system(const Vector_t& x, const Matrix_t& L, const Vector_t& fnl, const Vector_t& fex) {
    return fex - dot(L,x) - fnl;
}

Matrix_t HBM::corretion_system_jac(Vector_t& x, double freq, const Matrix_t& L, const Vector_t& fnl)
{
    Matrix_t jac = L;

    Vector_t fnl_k;
    typename Vector_t::iterator it_fnl_k;
    typename Matrix_t::iterator it_jac, end_jac = jac.end();
    typename Vector_t::const_iterator it_fnl;
    for (size_t k = 0; k < x.size(); ++k) {
        x[k] += _dx;
        fnl_k = funcnl(x,freq,dft);
        it_fnl_k = fnl_k.begin();
        it_jac = jac.begin();
        it_fnl = fnl.begin();
        while (it_jac < end_jac) {
            (*it_jac)[k] += (*it_fnl_k - *it_fnl)/_dx; // jac[i][k] += (fnl_k[i] - fnl[i])/_dx;
            ++it_jac;   ++it_fnl_k;   ++it_fnl;
        }
        x[k] -= _dx;
    }
    
    return jac;
}

} // namespace math