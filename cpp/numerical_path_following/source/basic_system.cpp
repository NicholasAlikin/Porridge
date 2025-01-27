#include "basic_system.hpp"

namespace npath
{

void BasicSystem::continuation_initialization(math::vector_t<double,2>& jac
                                             ,math::vector<double>& ynorm) {
      // ynorm.resize(response_norm_size());
}

void BasicSystem::process_total_increment(math::vector<double>&
                                        , math::vector<double>&) {
}

void BasicSystem::process_total_increment(math::vector<double>& Dy
                                        , math::vector<double> &dy
                                  , const math::vector<double>& predictor
                                  ,                    double) {
}

void BasicSystem::solution_not_found() {
}

void BasicSystem::step_back(double new_step) {
}



} // namespace npath
