/*Nodes structures*/
#pragma once

#include "linalg.hpp"
#include <functional>

namespace fem {

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

struct NodeLoadVar: public NodeLoad {
	std::function<void(math::vector<double>&, double)> calc;
	void operator()(double var);
};


/* All releases in node */
struct NodeAllReleases {
	static const uint8_t IsReleased = 1;
	static const uint8_t IsNotReleased = 0;
	
	// node with releases (from all elements)
	size_t node_globalid;  // element node global id
	// all released dofs
	math::vector<size_t> dofs; // element node local dofs, which are neaded to be released
	
	// Release(size_t elem_id, uint8_t node_id, const math::vector<uint8_t>& dofs);
	// Release(size_t elem_id, uint8_t node_id, math::vector<uint8_t>&& dofs);
	// Release(const Release& other);
	// Release(Release&& other);
	bool operator<(const NodeAllReleases& other) const;
	// sort by elements id or element node global id
	// struct {
	// 	bool operator() (const Release& a, const Release& b) {
	// 		return a.elem_id < b.elem_id;
	// 	}
	// } ReleaseLessElemets;
	
	// struct {
	// 	bool operator() (const Release& a, const Release& b) {
	// 		return a.node_globalid < b.node_globalid;
	// 	}
	// } ReleaseLessNodes;
};

} // namespace fem