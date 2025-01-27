/*Elements declarations*/
#pragma once

#include "fem_base.hpp"
#include "nodes.hpp"

#include <cassert>

namespace fem {

struct BaseElement;
struct Elem;
struct ElemNL;
struct ElemReleases;

struct BaseElement {
	static size_t EmptyID;

	size_t nodes_dofs; // elem id
	size_t propID;
	size_t matlID;

	math::vector<size_t> nodes; // nodes IDs
	math::vector<double> parameters;
	
	BaseElement() = default;
	BaseElement(const BaseElement& other);
	BaseElement(BaseElement&& other);
	BaseElement(size_t nodes_dofs, size_t propID, size_t matlID
			   ,const math::vector<size_t>& nodes);
	
	
	virtual void calc_parameters(const math::vector<Node>& nodes_info) {};
	virtual void calc_parameters(const math::vector<Node>& nodes_info
								,const ElemReleases& releases) {};

	virtual size_t nnodes() const = 0;
	virtual size_t ndofs_node() const = 0;
	size_t ndofs() const;
	
	
	virtual ~BaseElement() = default;


	
	virtual math::vector_t<double,2> stiffness(const math::vector<double>& property
											 , const math::vector<double>& material) const = 0;
	virtual math::vector_t<double,2> mass(const math::vector<double>& property
	 								    , const math::vector<double>& material) const = 0;


	virtual void tangentStiffness_innerLoad(math::vector_t<double,2>& K
										  , math::vector_t<double,1>& inner_load
									, const math::vector_t<double,1>& property
	 								, const math::vector_t<double,1>& material
                            		, const math::vector_t<double,3>& Rsum
									, const math::vector_t<double,1>& q ) const = 0;

	virtual void tangentMass_inertiaLoad(math::vector_t<double,2>& M
									  , math::vector_t<double,1>& inert_load
								, const math::vector_t<double,1>& property
	 						    , const math::vector_t<double,1>& material
                                , const math::vector_t<double,3>& Rsum
								, const math::vector_t<double,1>& q
								, const math::vector_t<double,1>& dqdt
								, const math::vector_t<double,1>& d2qdt2 ) const = 0;

	virtual void tangentMassGyro_inertiaLoad(math::vector_t<double,2>& M
									  , math::vector_t<double,2>& G
									  , math::vector_t<double,1>& inert_load
								, const math::vector_t<double,1>& property
	 						    , const math::vector_t<double,1>& material
                                , const math::vector_t<double,3>& Rsum
								, const math::vector_t<double,1>& q
								, const math::vector_t<double,1>& dqdt
								, const math::vector_t<double,1>& d2qdt2 ) const = 0;

	virtual void tangentStiffness_innerLoad(math::vector_t<double,2>& K
										  , math::vector_t<double,1>& inner_load
									, const math::vector_t<double,1>& property
	 								, const math::vector_t<double,1>& material
									, const math::vector_t<double,1>& q ) const = 0;

	virtual void tangentMass_inertiaLoad(math::vector_t<double,2>& M
									  , math::vector_t<double,1>& inert_load
								, const math::vector_t<double,1>& property
	 						    , const math::vector_t<double,1>& material
								, const math::vector_t<double,1>& q
								, const math::vector_t<double,1>& dqdt
								, const math::vector_t<double,1>& d2qdt2 ) const = 0;

	template <std::derived_from<BaseElement> El, std::random_access_iterator It>
	static auto prms_R0(It&& prmsIt)
				-> std::conditional_t<std::is_const_v<std::remove_reference_t<typename It::reference>>, const_block_t
																			,block_t>
	{
		std::conditional_t<std::is_const_v<std::remove_reference_t<typename It::reference>>,
				const_block_t,block_t> R0 = {
			{prmsIt + El::prms::R01x, prmsIt + El::prms::R02x}
		   ,{prmsIt + El::prms::R02x, prmsIt + El::prms::R03x}
		   ,{prmsIt + El::prms::R03x, prmsIt + El::prms::R03z+1}
		};
		return R0;
	}


	template <typename It>
    requires requires(It it) {
        {*it} -> std::same_as<size_t&>;
    }
    static size_t min_elem_dof(It beg, It end) {
        /*Min element greater than 0.*/
        size_t value = *beg;
        ++beg;
        while (beg < end) {
            if ( ((*beg) > 0 ) &&
                ( ((*beg) < value) || (value == 0) )
            )
                value = *beg;
            ++beg;
        }
        return value;
    }
};

// Linear elements
struct Elem: BaseElement {
	using BaseElement::BaseElement;
};
// Non-linear elements
struct ElemNL: BaseElement {
	using BaseElement::BaseElement;

private:
	static math::vector_t<double,2> stiffnessStatic(const math::vector<double>& property
											      , const math::vector<double>& material
											      , const math::vector<double>& parameters) { return {}; };
	math::vector_t<double,2> stiffness(const math::vector<double>& property
									 , const math::vector<double>& material) const { return {}; };
	
	static math::vector_t<double,2> massStatic(const math::vector<double>& property
	 								    	 , const math::vector<double>& material
									    	 , const math::vector<double>& parameters) { return {}; };
	math::vector_t<double,2> mass(const math::vector<double>& property
	 							, const math::vector<double>& material) const { return {}; };
	
	// virtual math::vector<double> internal_load(const math::vector<double>& property
	// 										 , const math::vector<double>& material
	// 										 , const math::vector<double>& parameters) const = 0;
    
};


struct ElemReleases {
	struct NodeRelease {
		size_t node_localid;  // element node local id
		math::vector<size_t> dofs; // element node local dofs, which are neaded to be released
		size_t shift = 0; // number of released other elements nodes in this node
	};
	size_t elem_id;   // elements global id
	math::vector<ElemReleases::NodeRelease> nodes;
	bool operator<(const ElemReleases& other) const;
};
} // namespace fem