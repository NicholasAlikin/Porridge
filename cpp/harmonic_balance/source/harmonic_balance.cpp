#include "harmonic_balance.hpp"

namespace npath {

HBM::HBM(const math::vector_t<double,2>& mass
       , const math::vector_t<double,2>& damp
       , const math::vector_t<double,2>& stif
       , const DFT& dft
       , const std::function<math::vector<double>(const math::vector<double>&,double,const DFT&)>& funcnl)
    : dft(dft), ndof(mass.size())
    , Mass(math::kron(dft.derivative2_basic,mass))
    , Damp(math::kron(dft.derivative_basic,damp))
    , Stif(math::kron(math::eye<double>(2*dft.H+1),stif))
    , funcnl(funcnl)
{}

HBM::HBM(const math::vector_t<double,2>& mass
       , const math::vector_t<double,2>& damp
       , const math::vector_t<double,2>& stif
       , const DFT& dft
       , const std::function<math::vector<double>(const math::vector<double>&,double,const DFT&)>& funcnl
       , const std::function<math::vector<double>(double,const DFT&)>& funcex)
    : dft(dft), ndof(mass.size())
    , Mass(math::kron(dft.derivative2_basic,mass))
    , Damp(math::kron(dft.derivative_basic,damp))
    , Stif(math::kron(math::eye<double>(2*dft.H+1),stif))
    , funcnl(funcnl), funcex(funcex)
{}

void HBM::system_response(math::vector<double>& fun, math::vector_t<double,2>& jac, math::vector<double>& y) {
    math::vector<double> x(y.begin(),y.end()-1);
    double freq = *(y.end()-1);
    math::vector_t<double,2> L = linear_system_dynamic_reaction(freq);
    math::vector<double> fnl = funcnl(x,freq,dft);
    math::vector<double> fex = funcex(freq,dft);
    fun = system(x,L,fnl,fex);
    system_jac(jac, x,freq,L,fnl);
}

void HBM::system_response_extended(math::vector<double>& fun, math::vector_t<double,2>& jac, math::vector<double>& y) {
    math::vector<double> x(y.begin(),y.end()-1);
    double freq = y.last();

    math::vector_t<double,2> L = linear_system_dynamic_reaction(freq);
    math::vector<double> fnl = funcnl(x,freq,dft);
    math::vector<double> fex = funcex(freq,dft);
    
    math::Slice fun_sl(fun.begin(),fun.begin()+x.size());
    fun_sl = system(x,L,fnl,fex);
    system_jac(jac, x,freq,L,fnl);
    math::vector<double> jac_freq;
    system_jac_freq(jac_freq, x,freq,fnl,fex);

    size_t last = x.size();
    auto it_jac = jac.begin();
    auto it_jac_freq = jac_freq.begin(), end_jac_freq = jac_freq.end();
    while (it_jac_freq < end_jac_freq) {
        (*it_jac)[last] = *it_jac_freq;
        ++it_jac;   ++it_jac_freq;
    }
}

void HBM::system_response(math::vector<double>& fun, math::vector_t<double,2>& jac, math::vector<double>& y, const math::vector<double>&) {
    system_response(fun,jac, y);
}
void HBM::system_response_extended(math::vector<double>& fun, math::vector_t<double,2>& jac, math::vector<double>& y, const math::vector<double>&) {
    system_response_extended(fun,jac, y);
}

void HBM::calculate_response_norm(const math::vector<double>& y
                                      , math::vector<double>& ynorm) {
    auto it = ynorm.begin();
    math::Slice<decltype(y.begin()),decltype(y.end())> yk;
    
    for (size_t k = 0; k < ndof; ++k, ++it) {
        yk.new_slice(y.begin()+k,y.end()-2,ndof);
        *it = norm(yk);
    }
}


math::vector_t<double,2> HBM::linear_system_dynamic_reaction(double freq) const {
    math::vector_t<double,2> L = (freq*freq)*Mass + freq*Damp + Stif;
    return L;
}

math::vector<double> HBM::system(const math::vector<double>& x, const math::vector_t<double,2>& L, const math::vector<double>& fnl, const math::vector<double>& fex) const {
    return fex - math::dot(L,x) - fnl;
}

void HBM::system_jac(math::vector_t<double,2>& jac, math::vector<double>& x, double freq, const math::vector_t<double,2>& L, const math::vector<double>& fnl) const
{
    // std::cout << "jac = " << size(jac)[0] << ", " << size(jac)[1] <<"; x = " << x.size() << "; L = " << size(L)[0] << ", " << size(L)[1] << "; fnl = " << fnl.size() << std::endl;
    math::vector<double> fnl_k;
    typename math::vector<double>::iterator it_fnl_k;
    typename math::vector_t<double,2>::iterator it_jac;
    typename math::vector_t<double,2>::const_iterator it_L;
    typename math::vector<double>::const_iterator it_fnl, end_fnl = fnl.end();
    for (size_t k = 0; k < x.size(); ++k) {
        x[k] += _dx;
        fnl_k = funcnl(x,freq,dft);
        it_fnl_k = fnl_k.begin();
        it_jac = jac.begin();
        it_fnl = fnl.begin();
        it_L = L.begin();
        while (it_fnl < end_fnl) {
            // std::cout << "k = " << k << std::endl;
            (*it_jac)[k] = (*it_L)[k] + (*it_fnl_k - *it_fnl)/_dx; // jac[i][k] += (fnl_k[i] - fnl[i])/_dx;
            ++it_jac;   ++it_fnl_k;   ++it_fnl;   ++it_L;
        }
        x[k] -= _dx;
    }
}


void HBM::system_jac_freq(math::vector<double>& jac, const math::vector<double>& x, double freq, const math::vector<double>& fnl, const math::vector<double>& fex) const {
    jac = dot((2*freq)*Mass + Damp,x) +
            ((funcnl(x,freq+_dx,dft) - funcex(freq+_dx,dft))
            -(fnl                    - fex                 ))/_dx;
}

size_t HBM::response_norm_size() {
    return ndof;
}

double HBM::fun_norm(const math::vector<double>& fun,
                    const math::vector<double>& y) {
    return math::norm(fun);
}

void HBM::continuation_initialization(math::vector_t<double,2>& jac
                                     ,math::vector<double>& ynorm) {
    
    // set ynorm size
    ynorm.resize(ndof);
}

} // namespace npath