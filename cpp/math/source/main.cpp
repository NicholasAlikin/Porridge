#include "linalg.hpp"
#include "dft.hpp"
#include "slice.hpp"
#include "corrector.hpp"

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

int main() {
    auto m = eye<double>(2)*5;
    auto d = eye<double>(2)*0.1;
    auto s = eye<double>(2)*2;
    auto v = repmat<double>(1,16);
    Slice sl(v.begin(),v.end(),2);
    std::cout << norm(sl);
    // DFT dft(1,2,3);
    // Corrector c(m,d,s,dft);
    // Corrector c2(m,d,s,8,1);
    // std::cout << c2.dft.ndof;
    
    // c.temp();
    // S s(std::move(v));
    // std::cout << v << '\n';
    // std::cout << s.v << '\n';
    // vector<int> v2 = {228,1337};
    // Slice sl(v.begin(),v.end(),2);
    // auto it = sl.begin();
    // sl = v2;
    // std::cout << v;
    // Slice sl(v.cbegin(),v.cend(),228);
    // auto it = sl.begin();
    // auto end = sl.end();

    // Debug<typename decltype(sl)::iterator::pointer>();

    // while (it < end) {
    //     std::cout << *it << ' ';
    //     ++it;
    // }
    
    

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
    Matrix_t T_forward = DFT::calculate_forward(N,H,ndof);
    Matrix_t T_backward = DFT::calculate_backward(N,H,ndof);
    auto ampl = dot(T_forward,f);
    std::cout << ampl << std::endl;
    auto f1 = dot(T_backward,ampl);
    std::cout << f - f1;*/
    // std::cout << for_each(v, std::cos<double>);
    // std::cout << T_forward << "\n\n\n";
    // std::cout << T_backward;
    

}