#pragma once

#include "elements.hpp"

namespace fem {


struct ElemBEAM2D: Elem {
	static constexpr short NNODES = 2;
	static constexpr short NDOFS_NODE = 2;
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
	static constexpr short NNODES = 2;
	static constexpr short NDOFS_NODE = 6;
	using Elem::Elem;
	enum prop { Iy,Iz,Ik,A,ky,kz,Jp };
	enum matl { E,mu,rho };
	enum prms { // orient vector
				vx,vy,vz
				// length
				,L 
				 // R0
				,R01x, R01y, R01z
				,R02x, R02y, R02z 
				,R03x, R03y, R03z};
	
	enum dof { uax,uay,uaz,tax,tay,taz,
			   ubx,uby,ubz,tbx,tby,tbz};
	
	ElemBEAM(size_t ID, size_t propID, size_t matlID
			   ,const math::vector<size_t>& nodes
			   ,const math::vector<double>& orient_vector);
	
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
	
	double prms_L() const;
	block_t prms_R0();
	const_block_t prms_R0() const;
	math::Slice<math::vector<double>::const_iterator
	           ,math::vector<double>::const_iterator> prms_orientVec() const;
	math::Slice<math::vector<double>::iterator
	           ,math::vector<double>::iterator> prms_orientVec();
	
};

struct ElemBEAMLD: ElemNL {
	static constexpr short NNODES = 2;
	static constexpr short NDOFS_NODE = 6;
	using ElemNL::ElemNL;
	using dof = ElemBEAM::dof;
	using prop = ElemBEAM::prop;
	using matl = ElemBEAM::matl;
	
	enum prms { // orient vector
				vx,vy,vz
				// length
				,L 
				 // R0
				,R01x, R01y, R01z // == ex0
				,R02x, R02y, R02z 
				,R03x, R03y, R03z
				
				
				// releases shifts for both nodes
				,release_shift1, release_shift2};
	
	ElemBEAMLD(size_t ID, size_t propID, size_t matlID
			   ,const math::vector<size_t>& nodes
			   ,const math::vector<double>& orient_vector);
	
	void tangentStiffness_innerLoad(math::vector_t<double,2>& K
								  , math::vector_t<double,1>& inner_load
							, const math::vector_t<double,1>& property
	 						, const math::vector_t<double,1>& material
                            , const math::vector_t<double,3>& Rsum
							, const math::vector_t<double,1>& q ) const override;

	void tangentMass_inertiaLoad(math::vector_t<double,2>& M
							  , math::vector_t<double,1>& inert_load
						, const math::vector_t<double,1>& property
					    , const math::vector_t<double,1>& material
    	                , const math::vector_t<double,3>& Rsum
						, const math::vector_t<double,1>& q
						, const math::vector_t<double,1>& dqdt
						, const math::vector_t<double,1>& d2qdt2 ) const override;

	void tangentMassGyro_inertiaLoad(math::vector_t<double,2>& M
							 	   , math::vector_t<double,2>& G
							 	   , math::vector_t<double,1>& inert_load
							 , const math::vector_t<double,1>& property
	 						 , const math::vector_t<double,1>& material
                        	 , const math::vector_t<double,3>& Rsum
							 , const math::vector_t<double,1>& q
							 , const math::vector_t<double,1>& dqdt
							 , const math::vector_t<double,1>& d2qdt2 ) const override;
private:
	void tangentStiffness_innerLoad(math::vector_t<double,2>& K
								  , math::vector_t<double,1>& inner_load
							, const math::vector_t<double,1>& property
	 						, const math::vector_t<double,1>& material
							, const math::vector_t<double,1>& q ) const override {};

	void tangentMass_inertiaLoad(math::vector_t<double,2>& M
							  , math::vector_t<double,1>& inert_load
						, const math::vector_t<double,1>& property
					    , const math::vector_t<double,1>& material
						, const math::vector_t<double,1>& q
						, const math::vector_t<double,1>& dqdt
						, const math::vector_t<double,1>& d2qdt2 ) const override {};
public:			

	void calc_parameters(const math::vector<Node>& elem_nodes) override;
	void calc_parameters(const math::vector<Node>& nodes_info
						,const ElemReleases& releases) override;
	
	size_t nnodes() const override;
	size_t ndofs_node() const override;
	
	double prms_L() const;
	block_t prms_R0();
	const_block_t prms_R0() const;
	math::Slice<math::vector<double>::const_iterator
	           ,math::vector<double>::const_iterator> prms_basis0_ex0() const;
	math::Slice<math::vector<double>::iterator
	           ,math::vector<double>::iterator> prms_basis0_ex0();
	math::Slice<math::vector<double>::const_iterator
	           ,math::vector<double>::const_iterator> prms_orientVec() const;
	math::Slice<math::vector<double>::iterator
	           ,math::vector<double>::iterator> prms_orientVec();
	size_t prms_shift1() const;
	size_t prms_shift2() const;
};

/* Corotation beam element with swap Euler vector when it is qual to 2pi*/
struct ElemBEAMLD2: ElemBEAMLD {

	void tangentStiffness_innerLoad(math::vector_t<double,2>& K
								  , math::vector_t<double,1>& inner_load
							, const math::vector_t<double,1>& property
	 						, const math::vector_t<double,1>& material
							, const math::vector_t<double,1>& q ) const override;

	void tangentMass_inertiaLoad(math::vector_t<double,2>& M
							  , math::vector_t<double,1>& inert_load
						, const math::vector_t<double,1>& property
					    , const math::vector_t<double,1>& material
						, const math::vector_t<double,1>& q
						, const math::vector_t<double,1>& dqdt
						, const math::vector_t<double,1>& d2qdt2 ) const override {};
						
private:
	void tangentStiffness_innerLoad(math::vector_t<double,2>& K
								  , math::vector_t<double,1>& inner_load
							, const math::vector_t<double,1>& property
	 						, const math::vector_t<double,1>& material
                            , const math::vector_t<double,3>& Rsum
							, const math::vector_t<double,1>& q ) const override {};

	void tangentMass_inertiaLoad(math::vector_t<double,2>& M
							  , math::vector_t<double,1>& inert_load
						, const math::vector_t<double,1>& property
					    , const math::vector_t<double,1>& material
    	                , const math::vector_t<double,3>& Rsum
						, const math::vector_t<double,1>& q
						, const math::vector_t<double,1>& dqdt
						, const math::vector_t<double,1>& d2qdt2 ) const override {};
public:
	
};

} // namespace fem 