#pragma once

#include "fem.hpp"

namespace fem {
	
using block_t = math::vector< math::Slice<typename math::vector<double>::iterator,
						                  typename math::vector<double>::iterator> >;

struct ElemPLATE: Elem {
	using Elem::Elem;
	enum prop { h,K6rot };
	enum matl { E,mu,rho };
	enum prms { L1,L2,L3,L4 };
	
	enum dof { uax,uay,uaz,tax,tay,taz,
			   ubx,uby,ubz,tbx,tby,tbz,
               ucx,ucy,ucz,tcx,tcy,tcz,
			   udx,udy,udz,tdx,tdy,tdz};
			   
	

	
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

struct ElemPLATELD: ElemNL {
	using ElemNL::ElemNL;
	using dof = ElemPLATE::dof;
	using prop = ElemPLATE::prop;
	using matl = ElemPLATE::matl;
	enum prms { L1,L2,L3,L4,
	            ax,ay,az,
				bx,by,bz,
				cx,cy,cz,
				dx,dy,dz };
	
	void tangentStiffness_internalLoad(math::vector_t<double,2>& K
											 , math::vector_t<double,1>& internal_load
										,const math::vector_t<double,1>& property
	 								    ,const math::vector_t<double,1>& material
                                        ,const math::vector_t<double,2>& R0
										,const math::vector_t<double,2>& basis0
                                        ,const math::vector_t<double,3>& Rsum
										,const math::vector_t<double,1>& q) const;

	void calc_parameters(const math::vector<Node>& elem_nodes) override;
	
	size_t nnodes() const override;
	size_t ndofs_node() const override;
};

} // namespace fem 