#include "element_traits.hpp"

namespace fem {

math::vector_t<double, 2> ElementTraits<ElementType::BEAM2D>::stiffness(double E, double I, double l) {
	math::vector_t<double,2> K = {{ 12.,  6., -12.,  6.},
								  {  6.,  4.,  -6.,  2.},
								  {-12., -6.,  12., -6.},
								  {  6.,  2.,  -6.,  4.}};
	K *= E*I/(l*l*l);
	return K;
}

math::vector_t<double, 2> ElementTraits<ElementType::BEAM2D>::mass(double rho, double A, double l) {
	math::vector_t<double,2> M = {{ 156., 22.,   54., -13.},
								  {  22.,  4.,   13.,  -3.},
								  {  54., 13.,  156., -22.},
								  { -13., -3.,  -22.,   4.}};
	M *= rho*A*l/420.;
	return M;
}

}// namespace fem