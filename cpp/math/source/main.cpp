#include "dft.hpp"
#include "slice.hpp"
#include "harmonic_balance.hpp"
#include <functional>

using namespace math;



template <typename T>
struct Debug {
    Debug() = delete;
};


// struct S {
//     vector<int> v;
//     S(vector<int>&& v);
// };

// S::S(vector<int>&& v): v(std::move(v)) {}

Vector_t funcnl(const Vector_t& x, double w, const DFT& dft, double f0) {
    Vector_t x_time = dot(dft.backward,f0*x);
    Vector_t dxdt_time = dot(dft.backward,dot(dft.derivative,f0*w*x));
    Vector_t fnl = x_time*x_time*x_time;
    return dot(dft.forward, fnl);
}

Vector_t funcex(double w, const DFT& dft, double f0) {
    Vector_t f = zeros<double>((2*dft.H+1)*dft.ndof);
    f[1] = f0;
    return f;
}

int main() {

    /*Matrix_t K = {{ 5,-4, 1, 0}
                 ,{-4, 6,-4, 1}
                 ,{ 1,-4, 6,-4}
                 ,{ 0, 1,-4, 5}};
    Vector_t f = {0,1,0,0};
    // vector<double> D(K.size());
    // auto L = K;
    // LDLT(K,D,L);
    // std::cout << D << '\n' << L << '\n' << transpose(L);
    std::cout << solve(K,f);*/
    int n = 20;
    auto m1 = eye<int>(n);
    auto m2 = eye<double>(n);
    std::cout << det(m2) << std::endl;
    std::cout << det(m1) << std::endl;

    /*Matrix_t m = {{ 1   }};
    Matrix_t d = {{ 0.1 }};
    Matrix_t k = {{ 1   }};
    size_t H = 3, N = 64, ndof = 1;
    double w = 3, T = 2*pi/w;
    // double dt = T/N;
    // auto t = linspace<double>(0,T-dt,N);
    
    DFT dft(H,N,ndof);

    
    double fex0 = 1;
    double fnl0 = 1.5;

    HBM hbm(m,d,k,dft,
            [fnl0](const Vector_t& x, double w, const DFT& dft) { return funcnl(x,w,dft,fnl0); },
            [fex0](double w, const DFT& dft) { return funcex(w,dft,fex0); } );
    
    Vector_t y0 = zeros<double>((2*H+1)*ndof+1);
    y0[y0.size()-1] = w;
    // std::cout << hbm.linear_system_dynamic_reaction(w);
    auto y = hbm.correction(y0);
    Slice x(y.begin(),y.end()-1,1);
    std::cout << y << norm(x);*/

    


    
    
    

    /*size_t H = 1, N = 8, ndof = 2;
    double w = 2;
    double T = 2*pi/w;
    double dt = T/N;
    auto t = linspace<double>(0,T-dt,N);
    
    // auto f = 1 - 3*cos(w*t) + 228.5*sin(w*t);
    auto f = zeros<double>(t.size()*2);
    Slice fsl(f.begin(),f.end(),ndof);
    Slice fsl2(f.begin()+1,f.end(),ndof);
    fsl =   1 - 3*cos(w*t) + 228.5*sin(w*t);
    fsl2 = 22 + 18.35*cos(w*t) - 99.115*sin(w*t);

    // std::cout << f << std::endl;
    // DFT dft;
    Matrix_t T_forward = DFT::calculate_forward(H,N,ndof);
    Matrix_t T_backward = DFT::calculate_backward(H,N,ndof);
    auto ampl = dot(T_forward,f);
    std::cout << ampl << std::endl;
    auto f1 = dot(T_backward,ampl);
    std::cout << norm(f - f1);*/
    // std::cout << for_each(v, std::cos<double>);
    // std::cout << T_forward << "\n\n\n";
    // std::cout << T_backward;
    

}