#include "elements.hpp"

namespace fem {

BaseElement::BaseElement(const BaseElement& other)
		: nodes_dofs(other.nodes_dofs),propID(other.propID),matlID(other.matlID)
		 ,nodes(other.nodes),parameters(other.parameters) {}
BaseElement::BaseElement(BaseElement&& other)
		: nodes_dofs(other.nodes_dofs),propID(other.propID),matlID(other.matlID)
		 ,nodes(std::move(other.nodes)),parameters(std::move(other.parameters)) {}
BaseElement::BaseElement(size_t nodes_dofs, size_t propID, size_t matlID
					   , const math::vector<size_t> &nodes)
		: nodes_dofs(nodes_dofs),propID(propID),matlID(matlID),nodes(nodes) {}

size_t BaseElement::ndofs() const {
	return nnodes()*ndofs_node();
}

bool ElemReleases::operator<(const ElemReleases &other) const {
    return elem_id < other.elem_id;
}


} // namespace fem