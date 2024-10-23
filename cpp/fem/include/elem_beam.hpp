#pragma once

#include "fem.hpp"

namespace fem {

struct ElemBEAM2D: Elem {
	using Elem::Elem;
	enum prop { I,A };
	enum matl { E,RHO };
	enum prms { L };


	
	static math::vector_t<double,2> stiffnessStatic(const math::vector<double>& property
											 , const math::vector<double>& material
											 , const math::vector<double>& parameters);
	math::vector_t<double,2> stiffness(const math::vector<double>& property
											 , const math::vector<double>& material) const override;
	
	static math::vector_t<double,2> massStatic(const math::vector<double>& property
	 								    , const math::vector<double>& material
									    , const math::vector<double>& parameters);
	math::vector_t<double,2> mass(const math::vector<double>& property
	 								    , const math::vector<double>& material) const override;
	
	void calc_parameters(const math::vector<Node>& elem_nodes) override;

	size_t nnodes() const override;
	size_t ndofs_node() const override;
};

using block_t = math::vector< math::Slice<typename math::vector<double>::iterator,
						                  typename math::vector<double>::iterator> >;

struct ElemBEAM: Elem {
	using Elem::Elem;
	enum prop { Iy,Iz,Ik,A,ky,kz,Jp };
	enum matl { E,mu,rho };
	enum prms { L };
	
	enum dof { uax,uay,uaz,tax,tay,taz,
			   ubx,uby,ubz,tbx,tby,tbz};
			   
	
	static math::vector_t<double,2> stiffnessStatic(const math::vector<double>& property
											 , const math::vector<double>& material
											 , const math::vector<double>& parameters);
	math::vector_t<double,2> stiffness(const math::vector<double>& property
											 , const math::vector<double>& material) const override;
	
	static math::vector_t<double,2> massStatic(const math::vector<double>& property
	 								    , const math::vector<double>& material
									    , const math::vector<double>& parameters);
	math::vector_t<double,2> mass(const math::vector<double>& property
	 								    , const math::vector<double>& material) const override;
	
	void calc_parameters(const math::vector<Node>& elem_nodes) override;

	size_t nnodes() const override;
	size_t ndofs_node() const override;
};

struct ElemBEAMLD: ElemNL {
	using ElemNL::ElemNL;
	using dof = ElemBEAM::dof;
	using prop = ElemBEAM::prop;
	using matl = ElemBEAM::matl;
	using prms = ElemBEAM::prms;
	
	void tangentStiffness_internalLoad(math::vector_t<double,2>& K
								   	    , math::vector_t<double,1>& internal_load
								   ,const math::vector_t<double,1>& property
	 							   ,const math::vector_t<double,1>& material
                                   ,const math::vector_t<double,2>& R0
								   ,const math::vector_t<double,2>& basis0
                                   ,const math::vector_t<double,3>& Rsum
								   ,const math::vector_t<double,1>& q ) const;

	void calc_parameters(const math::vector<Node>& elem_nodes) override;
};

} // namespace fem 