#include "analyzes.hpp"

namespace fem {

math::vector_t<double,1> staticLD(const math::vector<fem::Node>& nodes_info
				, const math::vector<fem::NodeConstraint>& constraints_info
				, const math::vector<fem::NodeLoad>& loads_info
				, math::vector<BaseElement*>& elements
				, const math::vector_t<double,2>& properties
				, const math::vector_t<double,2>& materials
				, size_t load_steps
                , double epsq = 1e-5)
{
	/*Static analysis with large displacements, linear elasticity.
	
	1. Inital iteration like for linear static analysis
	2. Load step
	2.1. Newton-Raphson iterations while displacement and/or load convergence.

	*/

	// precomputing
	fem::GlobalDofs GDofs = fem::parse_nodal_data(nodes_info, constraints_info, elements);
	
    assemble_precomputing(GDofs,elements);
    // initialazing stiffness matrix, load vector and displacement vector
	math::vector_t<double> stif(math::sum(GDofs.colhs))
                        ,  loadExt(GDofs.ndofs)
                        ,  q(GDofs.ndofs)
                        ,  loadInt(GDofs.ndofs)
                        ,  dq;
    
    // Calculate exact load vector
    store_load_vector(loadExt,GDofs,loads_info);
    double loadExt_abs = norm(loadExt);
    double loadExt_step_abs = loadExt_abs/load_steps;
    // load /= load_abs;

    // Inital value of q equal to 0, so zeros iteration like static linear analysis.
    
    // assemble
    loadExt *= loadExt_step_abs/loadExt_abs;
    q = math::solve2(stif,loadExt,GDofs.diags);

    
    for (size_t iter = 0; iter < load_steps; ++iter) {
        // assemble stiffness matrix and internal load vector
        assembleNL();
        // update load external load vector
        dq = math::solve2(stiff,loadExt*() - loadInt,GDofs.diags);
        q += dq;
        while (math::norm(dq) < epsq) {
            assembleNL();
            dq = math::solve2(stiff,loadExt*() - loadInt,GDofs.diags);
            q += dq;
        }
    }

}



} // namespace fem