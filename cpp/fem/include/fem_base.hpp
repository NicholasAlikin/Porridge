/*Base FEM module with all declarations*/

#pragma once

#include "linalg.hpp"
#include "tensor_operations.hpp"


namespace fem {

using block_t = math::vector< math::Slice<typename math::vector<double>::iterator,
						                  typename math::vector<double>::iterator> >;
using const_block_t = math::vector< math::Slice<typename math::vector<double>::const_iterator,
												typename math::vector<double>::const_iterator> >;

class Model;
class Assemble;
struct ModelTraits;
struct AnalysisTraits;

} // namespace fem