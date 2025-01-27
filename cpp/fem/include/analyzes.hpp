#pragma once

#include "fem.hpp"
#include "linalg.hpp"

namespace fem {

struct AnalysisTraits {
	static void setup_analysis(const Model& model, Assemble& assemble);

	static math::vector_t<double,1> staticLD(Model& model
											,Assemble& assemble
											,size_t load_steps
											,double epsq = 1e-5
											,double epsload = 1e-5);
	static math::vector_t<double,1> staticLD2(Model& model
											,Assemble& assemble
											,size_t load_steps
											,double epsq = 1e-5
											,double epsload = 1e-5);

	static void update_Rsum( const Model& model
							,const Assemble& assemble
							,math::vector_t<double,3>& Rsum
							,math::vector<double>& q
							,math::vector<double>& temp_theta
							,math::vector_t<double,2>& temp_rotTensor
							,math::vector_t<double,2>& temp_Rsumi);
	static math::vector_t<double,3> setup_Rsum(const Model& model);
private:
	static void update_Rsumi_releases(math::vector<double>& q
					, typename math::vector_t<double,2>::iterator Rsumi
					, typename math::vector<size_t>::const_iterator gdof
					, typename math::vector<size_t>::const_iterator subdof
					, 		   math::vector<double>& theta
					, math::vector_t<double,2>& temp_rotTensor
					, math::vector_t<double,2>& temp_Rsumi);
	static void update_Rsumi_releases(math::vector<double>& q
					, typename math::vector_t<double,2>::iterator Rsumi
					, typename math::vector<size_t>::const_iterator gdof
					, 		   math::vector<double>& theta
					, math::vector_t<double,2>& temp_rotTensor
					, math::vector_t<double,2>& temp_Rsumi);
	static void update_Rsumi(math::vector<double>& q
					, typename math::vector_t<double,3>::iterator Rsumi
					, typename math::vector<size_t>::const_iterator gdof
					, 		   math::vector<double>& theta
					, math::vector_t<double,2>& temp_rotTensor
					, math::vector_t<double,2>& temp_Rsumi);

};


math::vector_t<double,1> staticLD2(Model& model
								, Assemble& assemble
				                , size_t load_steps
								, double epsq = 1e-5
								, double epsload = 1e-5);

void update_loadExt_vector(Model& model
						, const Assemble& assemble
						, math::vector<double>& loadExt
						, double& cur_loadExt_norm
						, double new_loadExt_norm);


} // namespace fem