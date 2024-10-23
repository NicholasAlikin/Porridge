#pragma once

#include "linalg.hpp"
#include <iostream>
#include <array>
#include "element_traits.hpp"

namespace fem {



struct Point;
struct Line;

struct BaseNode;
struct Node;
struct NodeConstraint;
struct NodeLoad;

struct GlobalDofs;

struct BaseElement;

template <int ElemType>
struct Element;

struct Elements;


struct Point {
    static size_t EmptyID;

};

struct BaseNode {
	static const uint8_t DIM = 3;
    static const uint8_t DOFS = 6;
	
	size_t id;
	
	bool operator<(const BaseNode& other) const;
};

struct Node: public BaseNode {
    math::vector<double> xyz;
};

struct NodeConstraint: public BaseNode {
    math::vector<uint8_t> dofs;
};

struct NodeLoad: public BaseNode {
	math::vector<double> dofs;
};

// template <typename NodeLike>
// class NodeContainer {
// /*
	// Container for store Nodes like objects: Node, NodeConstraint, NodeLoad, etc.
	// For container using std::unordered_map: key - Node id, value - Node object.
	
	// O(1) to get node by id.
	// Iterator to go throw the all Nodes in container.
	
	// Nodes are stored sorted by id.
// */
// private:
	// std::multimap<size_t,NodeLike> container;
// public:
	
// };

struct GlobalDofs {
	size_t ndofs;
	math::vector_t<size_t,2> ID;
	math::vector_t<size_t,2> ElemID;
	math::vector<size_t> colhs;
	math::vector<size_t> diags;
	size_t band_width;
	math::vector<double> stiffness;
	
	math::vector_t<double,3> plate_basis;

	GlobalDofs() = default;
	GlobalDofs(size_t ndofs, const math::vector_t<size_t,2>& ID);
	GlobalDofs(size_t ndofs, math::vector_t<size_t,2>&& ID);
	GlobalDofs(const GlobalDofs& other);
	GlobalDofs(GlobalDofs&& other);
};

/*
Elements are defined using virtual functions.

*/

struct BaseElement {
	static size_t EmptyID;

	size_t ID; // elem id
	size_t propID;
	size_t matlID;

	math::vector<size_t> nodes; // nodes IDs
	math::vector<double> parameters;
	
	BaseElement() = default;
	BaseElement(const BaseElement& other);
	BaseElement(BaseElement&& other);
	BaseElement(size_t ID=0, size_t propID=0, size_t matlID=0
			   ,const math::vector<size_t>& nodes = {}
			   ,const math::vector<double>& parameters = {});
	
	
	virtual void calc_parameters(const math::vector<Node>& elem_nodes) = 0;

	virtual size_t nnodes() const = 0;
	virtual size_t ndofs_node() const = 0;
	size_t ndofs() const;
	
	void set_nodes(const math::vector<size_t>& nodes);
	void set_nodes(math::vector<size_t>&& nodes);

	virtual ~BaseElement() = default;
};

// Linear elements
struct Elem: BaseElement {
	using BaseElement::BaseElement;

	static math::vector_t<double,2> stiffnessStatic(const math::vector<double>& property
											 , const math::vector<double>& material
											 , const math::vector<double>& parameters);
	virtual math::vector_t<double,2> stiffness(const math::vector<double>& property
											 , const math::vector<double>& material) const = 0;
	
	static math::vector_t<double,2> massStatic(const math::vector<double>& property
	 								    , const math::vector<double>& material
									    , const math::vector<double>& parameters);
	virtual math::vector_t<double,2> mass(const math::vector<double>& property
	 								    , const math::vector<double>& material) const = 0;
};
// Non-linear elements
struct ElemNL: BaseElement {
	using BaseElement::BaseElement;
	
	virtual void tangentStiffness_internalLoad(math::vector_t<double,2>& K
											 , math::vector_t<double,1>& internal_load
										,const math::vector_t<double,1>& property
	 								    ,const math::vector_t<double,1>& material
                                        ,const math::vector_t<double,2>& R0
										,const math::vector_t<double,2>& basis0
                                        ,const math::vector_t<double,3>& Rsum
										,const math::vector_t<double,1>& q ) const = 0;
	// virtual math::vector<double> internal_load(const math::vector<double>& property
	// 										 , const math::vector<double>& material
	// 										 , const math::vector<double>& parameters) const = 0;
};

template <typename It>
requires requires(It it) {
	{*it} -> std::same_as<size_t&>;
}
size_t min_elem_dof(It beg, It end) {
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


GlobalDofs parse_nodal_data(const math::vector<Node>& nodes_info
                          , const math::vector<NodeConstraint>& constraints_info
						  , math::vector<BaseElement*>& elements );
										
void store_load_vector(math::vector<double> load
									 , const GlobalDofs& ID
									 , const math::vector<NodeLoad>& loads_info);
									 
void assemble_precomputing(GlobalDofs& GDofs, const math::vector<BaseElement*>& elements);
math::vector_t<double,2> assemble(GlobalDofs& GDofs, const math::vector<Elem*>& elements
				, const math::vector_t<double,2>& properties
				, const math::vector_t<double,2>& materials
				, math::vector_t<double>& band
				, math::vector_t<double,2> (Elem::* element_matrix)(const math::vector<double>&,const math::vector<double>&) const);

math::vector_t<double,2> assembleNL(GlobalDofs& GDofs, const math::vector<ElemNL*>& elements
				, const math::vector_t<double,2>& properties
				, const math::vector_t<double,2>& materials
				, math::vector_t<double>& band
				, math::vector_t<double,2> (ElemNL::* element_matrix)(       math::vector_t<double,2>&
                                                                      ,      math::vector_t<double,1>&
                                                                      ,const math::vector_t<double,1>&
                                                                      ,const math::vector_t<double,1>&
                                                                      ,const math::vector_t<double,1>&
                                                                      ,const math::vector_t<double,2>&
                                                                      ,const math::vector_t<double,2>&
                                                                      ,const math::vector_t<double,3>&
                                                                      ,const math::vector_t<double,1>&) const);



void solveLD(GlobalDofs& GDofs, const math::vector<BaseElement*>& elements
				, const math::vector_t<double,2>& properties
				, const math::vector_t<double,2>& materials
				, math::vector_t<double>& band);

} // namespace fem