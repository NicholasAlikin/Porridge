#pragma once

#include "model.hpp"
#include "elements.hpp"
#include "nodes.hpp"

namespace fem {


struct Point;
struct Line;


struct Point {
    static size_t EmptyID;

};


struct StateQ {
	math::vector<double> displacement;
};

struct StateQV: StateQ {
	math::vector<double> velocity;
};

struct StateQVA: StateQV {
	math::vector<double> acceleration;
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





/*
Elements are defined using virtual functions.

*/





										
// void store_load_vector(math::vector<double>& load
// 									 , const Assemble& nodes_dofs
// 									 , const math::vector<NodeLoad>& loads_info);
									 
// void assemble_precomputing(Assemble& GDofs, const math::vector<BaseElement*>& elements
// 											, const math::vector<ElemReleases>& releases = {});
// math::vector_t<double,2> assemble(Assemble& GDofs, const math::vector<BaseElement*>& elements
// 				, const math::vector_t<double,2>& properties
// 				, const math::vector_t<double,2>& materials
// 				, math::vector_t<double>& band
// 				, math::vector_t<double,2> (BaseElement::* element_matrix)(const math::vector<double>&,const math::vector<double>&) const);

// void assembleNL(const Assemble& GDofs, const math::vector<BaseElement*>& elements
// 				, const math::vector_t<double,2>& properties
// 				, const math::vector_t<double,2>& materials
// 				, math::vector<double>& band
// 				, math::vector<double>& load
// 				, void (BaseElement::* element_matrix_load)(    math::vector_t<double,2>& /* K */
// 													,      math::vector_t<double,1>& /* internal_load */ 
// 													,const math::vector_t<double,1>& /* property */ 
// 													,const math::vector_t<double,1>& /* material */ 
// 													,const math::vector_t<double,3>& /* Rsum */ 
// 													,const math::vector_t<double,1>& /* q */  ) const
// 				, const math::vector<double>& q
// 				, const math::vector_t<double,3>& Rsum
// 				, math::vector_t<double,2>& matrix_local
// 				, math::vector<double>& load_local
// 				, math::vector<double>& q_elem);

// void assembleNL(const Assemble& GDofs, const math::vector<BaseElement*>& elements
// 				, const math::vector_t<double,2>& properties
// 				, const math::vector_t<double,2>& materials
// 				, math::vector_t<double,2>& matrix
// 				, math::vector<double>& load
// 				, void (BaseElement::* element_matrix_load)(    math::vector_t<double,2>& /* K */
// 													,      math::vector_t<double,1>& /* internal_load */ 
// 													,const math::vector_t<double,1>& /* property */ 
// 													,const math::vector_t<double,1>& /* material */ 
// 													,const math::vector_t<double,3>& /* Rsum */ 
// 													,const math::vector_t<double,1>& /* q */  ) const
// 				, const math::vector<double>& q
// 				, const math::vector_t<double,3>& Rsum);


} // namespace fem