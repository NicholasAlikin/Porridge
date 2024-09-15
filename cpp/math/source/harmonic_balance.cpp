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

void HBM::system_response(Vector_t& fun, Matrix_t& jac, const Vector_t& y) const {
    Vector_t x(y.begin(),y.end()-1);
    double freq = *(y.end()-1);
    Matrix_t L = linear_system_dynamic_reaction(freq);
    Vector_t fnl = funcnl(x,freq,dft);
    Vector_t fex = funcex(freq,dft);
    fun = system(x,L,fnl,fex);
    system_jac(jac, x,freq,L,fnl);
}

void HBM::system_response_ext(Vector_t& fun, Matrix_t& jac, const Vector_t& y) const {
    Vector_t x(y.begin(),y.end()-1);
    double freq = *(y.end()-1);
    Matrix_t L = linear_system_dynamic_reaction(freq);
    Vector_t fnl = funcnl(x,freq,dft);
    Vector_t fex = funcex(freq,dft);
    Slice fun_sl(fun.begin(),fun.begin()+x.size());
    fun_sl = system(x,L,fnl,fex);
    system_jac(jac, x,freq,L,fnl);
    Vector_t jac_freq;
    system_jac_freq(jac_freq, x,freq,fnl,fex);

    size_t last = x.size();
    auto it_jac = jac.begin();
    auto it_jac_freq = jac_freq.begin(), end_jac_freq = jac_freq.end();
    while (it_jac_freq < end_jac_freq) {
        (*it_jac)[last] = *it_jac_freq;
        ++it_jac;   ++it_jac_freq;
    }
}

void HBM::system_response(Vector_t& fun, Matrix_t& jac, const Vector_t& y, const Vector_t&) const {
    system_response(fun,jac, y);
}
void HBM::system_response_ext(Vector_t& fun, Matrix_t& jac, const Vector_t& y, const Vector_t&) const {
    system_response_ext(fun,jac, y);
}


Matrix_t HBM::linear_system_dynamic_reaction(double freq) const {
    Matrix_t L = (freq*freq)*Mass + freq*Damp + Stif;
    return L;
}

Vector_t HBM::system(const Vector_t& x, const Matrix_t& L, const Vector_t& fnl, const Vector_t& fex) const {
    return fex - dot(L,x) - fnl;
}

void HBM::system_jac(Matrix_t& jac, Vector_t& x, double freq, const Matrix_t& L, const Vector_t& fnl) const
{
    Vector_t fnl_k;
    typename Vector_t::iterator it_fnl_k;
    typename Matrix_t::iterator it_jac;
    typename Matrix_t::const_iterator it_L;
    typename Vector_t::const_iterator it_fnl, end_fnl = fnl.end();
    for (size_t k = 0; k < x.size(); ++k) {
        x[k] += _dx;
        fnl_k = funcnl(x,freq,dft);
        it_fnl_k = fnl_k.begin();
        it_jac = jac.begin();
        it_fnl = fnl.begin();
        it_L = L.begin();
        while (it_fnl < end_fnl) {
            (*it_jac)[k] = (*it_L)[k] + (*it_fnl_k - *it_fnl)/_dx; // jac[i][k] += (fnl_k[i] - fnl[i])/_dx;
            ++it_jac;   ++it_fnl_k;   ++it_fnl;   ++it_L;
        }
        x[k] -= _dx;
    }
}


void HBM::system_jac_freq(Vector_t& jac, const Vector_t& x, double freq, const Vector_t& fnl, const Vector_t& fex) const {
    jac = dot((2*freq)*Mass + Damp,x) +
            ((funcnl(x,freq+_dx,dft) - funcex(freq+_dx,dft))
            -(fnl                    - fex                 ))/_dx;
}

} // namespace math