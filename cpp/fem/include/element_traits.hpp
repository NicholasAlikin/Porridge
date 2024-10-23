#pragma once

#include "linalg.hpp"

namespace fem {
	

struct ElementType {
	static const int BEAM2D = 1;
};

template <int ElemType>
struct ElementTraits {
	static math::vector_t<double,2> stiffness();
	static math::vector_t<double,2> mass;
	static math::vector_t<double,2> damp;
};

template <>
struct ElementTraits<ElementType::BEAM2D> {
	static const size_t nnodes = 2;
	static const size_t ndofs_node = 2;
	static const size_t ndofs = nnodes*ndofs_node;
	
	static math::vector_t<double,2> stiffness(double E, double I, double l);
	static math::vector_t<double,2> mass(double rho, double A, double l);

};

struct MaterialType {
	static const int ISOTROPIC = 1;
};

template <int MetlType>
struct MaterialTraits;

template <>
struct MaterialTraits<MaterialType::ISOTROPIC> {
	/* E, mu, rho */
	static const int num_of_consts = 3;
};

	
} // namespace fem