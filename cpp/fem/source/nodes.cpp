#include "nodes.hpp"

namespace fem {

bool BaseNode::operator<(const BaseNode& other) const {
	return id < other.id;
}

void NodeLoadVar::operator()(double var) {
    calc(dofs,var);
}

bool NodeAllReleases::operator<(const NodeAllReleases &other) const {
    return node_globalid < other.node_globalid;
}


} // namespace fem