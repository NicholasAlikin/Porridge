#pragma once
#include "linalg.hpp"

namespace npath {

class DFT {

public:
    size_t H = 0;
    size_t N = 0;
    size_t ndof = 1;

    // math::vector_t<double,2> forward_basic;     // for ndof = 1
    math::vector_t<double,2> forward_basic;
    math::vector_t<double,2> backward_basic;
    math::vector_t<double,2> derivative_basic;  // for ndof = 1
    math::vector_t<double,2> derivative2_basic; // for ndof = 1

    // math::vector_t<double,2> forward;
    // math::vector_t<double,2> backward;
    // math::vector_t<double,2> derivative;
    // math::vector_t<double,2> derivative2;


    static math::vector_t<double,2> calculate_forward_basic(size_t H, size_t N);
    static math::vector_t<double,2> calculate_backward_basic(size_t H, size_t N);
    static math::vector_t<double,2> calculate_derivative_basic(size_t H, size_t N);
    static math::vector_t<double,2> calculate_derivative2_basic(size_t H, size_t N);
    static math::vector_t<double,2> calculate_derivative2_basic(const math::vector_t<double,2>& T_derivative_basic);


    static math::vector_t<double,2> calculate_forward(size_t H, size_t N, size_t ndof);
    static math::vector_t<double,2> calculate_backward(size_t H, size_t N, size_t ndof);
    static math::vector_t<double,2> calculate_derivative(size_t H, size_t N, size_t ndof);


    DFT() = default;
    DFT(size_t H, size_t N, size_t ndof = 1);

private:
    math::vector<double> transform(const math::vector_t<double,2>& transform_matrix, size_t res_size, const math::vector<double>& vec) const;
    void __derivative(math::vector<double>& res) const;
    void __derivative2(math::vector<double>& res) const;
public:
    math::vector<double> forward(const math::vector<double>& vec) const;
    math::vector<double> backward(const math::vector<double>& vec) const;
    math::vector<double> derivative(const math::vector<double>& vec) const;
    math::vector<double> derivative(math::vector<double>&& vec) const;
    math::vector<double> derivative2(const math::vector<double>& vec) const;
    math::vector<double> derivative2(math::vector<double>&& vec) const;
    

};

} // namespace npath