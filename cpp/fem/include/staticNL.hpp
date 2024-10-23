// #pragma once

// #include "corrector.hpp"
// #include <functional>
// #include "fem.hpp"

// namespace fem {


// class StaticNL: public math::Corrector {
// private:
// public:
    
    
//     StaticNL() = default;

//     math::Vector_t analysis();

//     void system_response(math::Vector_t& fun, math::Matrix_t& jac, const math::Vector_t& y) const override;
//     void system_response(math::Vector_t& fun, math::Matrix_t& jac, const math::Vector_t& y, const math::Vector_t& Dy) const override;
//     void system_response_ext(math::Vector_t& fun, math::Matrix_t& jac, const math::Vector_t& y) const override;
//     void system_response_ext(math::Vector_t& fun, math::Matrix_t& jac, const math::Vector_t& y, const math::Vector_t& Dy) const override;

//     // void correction_system_response(Vector_t& fun, Matrix_t& jac, const Vector_t& y) override;
//     // void continuation_system_response(Vector_t& fun, Matrix_t& jac, const Vector_t& y, const Vector_t& Dy) override;

    
    
//     math::Vector_t system(const math::Vector_t& x, const math::Matrix_t& L, const math::Vector_t& fnl, const math::Vector_t& fex) const;
//     void system_jac(math::Matrix_t& jac, math::Vector_t& x, double freq, const math::Matrix_t& L, const math::Vector_t& fnl) const;
//     void system_jac_load(math::Vector_t& jac, const math::Vector_t& x, double freq, const math::Vector_t& fnl, const math::Vector_t& fex) const;
// };

// } // namespace math