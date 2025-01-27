#include "elem_beam.hpp"

namespace fem {


math::vector_t<double,2> ElemBEAM2D::stiffnessStatic(const math::vector<double>& property
	 								 		 , const math::vector<double>& material
									 		 , const math::vector<double>& parameters) {
	math::vector_t<double,2> K = {{ 12.,  6., -12.,  6.},
								  {  6.,  4.,  -6.,  2.},
								  {-12., -6.,  12., -6.},
								  {  6.,  2.,  -6.,  4.}};
	K *= material[ElemBEAM2D::matl::E]
	    *property[ElemBEAM2D::prop::I]
	    /std::pow(
		parameters[ElemBEAM2D::prms::L]
		      ,3);
	return K;
}
math::vector_t<double,2> ElemBEAM2D::stiffness(const math::vector<double>& property
											 , const math::vector<double>& material) const
{
	return stiffnessStatic(property,material,parameters);
}

math::vector_t<double,2> ElemBEAM2D::massStatic(const math::vector<double>& property
	 								 	, const math::vector<double>& material
									 	, const math::vector<double>& parameters)
{
	math::vector_t<double,2> M = {{ 156., 22.,   54., -13.},
								  {  22.,  4.,   13.,  -3.},
								  {  54., 13.,  156., -22.},
								  { -13., -3.,  -22.,   4.}};
	M *= material[ElemBEAM2D::matl::RHO]
	    *property[ElemBEAM2D::prop::A]
		*parameters[ElemBEAM2D::prms::L]/420.;
	return M;
}

math::vector_t<double,2> ElemBEAM2D::mass(const math::vector<double>& property
	 								 	, const math::vector<double>& material) const
{
	return massStatic(property,material,parameters);	
}

void ElemBEAM2D::calc_parameters(const math::vector<Node>& nodes_info) {
	parameters = {
		// length
		math::norm(nodes_info[nodes[0]].xyz - nodes_info[nodes[1]].xyz)
	};
}

size_t ElemBEAM2D::nnodes() const {
	return ElemBEAM2D::NNODES;
}
size_t ElemBEAM2D::ndofs_node() const {
	return ElemBEAM2D::NDOFS_NODE;
}

ElemBEAM::ElemBEAM(size_t ID, size_t propID, size_t matlID
			   ,const math::vector<size_t>& nodes
			   ,const math::vector<double>& orient_vector)
		: ElemBEAM(ID,propID,matlID,nodes) {
	parameters = orient_vector / math::norm(orient_vector);
}

math::vector_t<double,2> ElemBEAM::stiffnessStatic(const math::vector<double>& property
	 								 		 , const math::vector<double>& material
									 		 , const math::vector<double>& parameters)
{
	double E = material[ElemBEAM::matl::E];
    double mu = material[ElemBEAM::matl::mu];
	double G = E/(2*(1+mu));
    double Iy = property[ElemBEAM::prop::Iy];
    double Iz = property[ElemBEAM::prop::Iz];
    double Ik = property[ElemBEAM::prop::Ik];
    double A = property[ElemBEAM::prop::A];
	
    double L = parameters[ElemBEAM::prms::L];
    double L2 = L*L;
    double L3 = L2*L;

#if 0
	double ky = property[ElemBEAM::prop::ky];
	double kz = property[ElemBEAM::prop::kz];
    math::vector_t<double,2> K = {
		{ (A*E)/L,                                           0,                                           0,         0,                                                        0,                                                        0, -(A*E)/L,                                           0,                                           0,         0,                                                        0,                                                        0},
		{       0,  (12*A*E*G*Iz*ky)/(L*(12*E*Iz + A*G*L2*ky)),                                           0,         0,                                                        0,                    (6*A*E*G*Iz*ky)/(12*E*Iz + A*G*L2*ky),        0, -(12*A*E*G*Iz*ky)/(L*(12*E*Iz + A*G*L2*ky)),                                           0,         0,                                                        0,                    (6*A*E*G*Iz*ky)/(12*E*Iz + A*G*L2*ky)},
		{       0,                                           0,  (12*A*E*G*Iy*kz)/(L*(12*E*Iy + A*G*L2*kz)),         0,                   -(6*A*E*G*Iy*kz)/(12*E*Iy + A*G*L2*kz),                                                        0,        0,                                           0, -(12*A*E*G*Iy*kz)/(L*(12*E*Iy + A*G*L2*kz)),         0,                   -(6*A*E*G*Iy*kz)/(12*E*Iy + A*G*L2*kz),                                                        0},
		{       0,                                           0,                                           0,  (G*Ik)/L,                                                        0,                                                        0,        0,                                           0,                                           0, -(G*Ik)/L,                                                        0,                                                        0},
		{       0,                                           0,      -(6*A*E*G*Iy*kz)/(12*E*Iy + A*G*L2*kz),         0,  (4*E*Iy*(3*E*Iy + A*G*L2*kz))/(L*(12*E*Iy + A*G*L2*kz)),                                                        0,        0,                                           0,       (6*A*E*G*Iy*kz)/(12*E*Iy + A*G*L2*kz),         0, -(2*E*Iy*(6*E*Iy - A*G*L2*kz))/(L*(12*E*Iy + A*G*L2*kz)),                                                        0},
		{       0,       (6*A*E*G*Iz*ky)/(12*E*Iz + A*G*L2*ky),                                           0,         0,                                                        0,  (4*E*Iz*(3*E*Iz + A*G*L2*ky))/(L*(12*E*Iz + A*G*L2*ky)),        0,      -(6*A*E*G*Iz*ky)/(12*E*Iz + A*G*L2*ky),                                           0,         0,                                                        0, -(2*E*Iz*(6*E*Iz - A*G*L2*ky))/(L*(12*E*Iz + A*G*L2*ky))},
		{-(A*E)/L,                                           0,                                           0,         0,                                                        0,                                                        0,  (A*E)/L,                                           0,                                           0,         0,                                                        0,                                                        0},
		{       0, -(12*A*E*G*Iz*ky)/(L*(12*E*Iz + A*G*L2*ky)),                                           0,         0,                                                        0,                   -(6*A*E*G*Iz*ky)/(12*E*Iz + A*G*L2*ky),        0,  (12*A*E*G*Iz*ky)/(L*(12*E*Iz + A*G*L2*ky)),                                           0,         0,                                                        0,                   -(6*A*E*G*Iz*ky)/(12*E*Iz + A*G*L2*ky)},
		{       0,                                           0, -(12*A*E*G*Iy*kz)/(L*(12*E*Iy + A*G*L2*kz)),         0,                    (6*A*E*G*Iy*kz)/(12*E*Iy + A*G*L2*kz),                                                        0,        0,                                           0,  (12*A*E*G*Iy*kz)/(L*(12*E*Iy + A*G*L2*kz)),         0,                    (6*A*E*G*Iy*kz)/(12*E*Iy + A*G*L2*kz),                                                        0},
		{       0,                                           0,                                           0, -(G*Ik)/L,                                                        0,                                                        0,        0,                                           0,                                           0,  (G*Ik)/L,                                                        0,                                                        0},
		{       0,                                           0,      -(6*A*E*G*Iy*kz)/(12*E*Iy + A*G*L2*kz),         0, -(2*E*Iy*(6*E*Iy - A*G*L2*kz))/(L*(12*E*Iy + A*G*L2*kz)),                                                        0,        0,                                           0,       (6*A*E*G*Iy*kz)/(12*E*Iy + A*G*L2*kz),         0,  (4*E*Iy*(3*E*Iy + A*G*L2*kz))/(L*(12*E*Iy + A*G*L2*kz)),                                                        0},
		{       0,       (6*A*E*G*Iz*ky)/(12*E*Iz + A*G*L2*ky),                                           0,         0,                                                        0, -(2*E*Iz*(6*E*Iz - A*G*L2*ky))/(L*(12*E*Iz + A*G*L2*ky)),        0,      -(6*A*E*G*Iz*ky)/(12*E*Iz + A*G*L2*ky),                                           0,         0,                                                        0,  (4*E*Iz*(3*E*Iz + A*G*L2*ky))/(L*(12*E*Iz + A*G*L2*ky))}
	};
#else
	math::vector_t<double,2> K = {
		{ (A*E)/L,              0,              0,       0,             0,             0, -(A*E)/L,              0,             0,       0,            0,            0},
		{       0,   (12*E*Iz)/L3,              0,       0,             0,   (6*E*Iz)/L2,         0, -(12*E*Iz)/L3,             0,       0,            0,  (6*E*Iz)/L2},
		{       0,              0,   (12*E*Iy)/L3,       0,  -(6*E*Iy)/L2,             0,         0,             0, -(12*E*Iy)/L3,       0, -(6*E*Iy)/L2,            0},
		{       0,              0,              0,  G*Ik/L,             0,             0,         0,             0,             0, -G*Ik/L,            0,            0},
		{       0,              0,   -(6*E*Iy)/L2,       0,    (4*E*Iy)/L,             0,         0,             0,   (6*E*Iy)/L2,       0,   (2*E*Iy)/L,            0},
		{       0,    (6*E*Iz)/L2,              0,       0,             0,    (4*E*Iz)/L,         0,  -(6*E*Iz)/L2,             0,       0,            0,   (2*E*Iz)/L},
		{-(A*E)/L,              0,              0,       0,             0,             0,   (A*E)/L,             0,             0,       0,            0,            0},
		{       0,  -(12*E*Iz)/L3,              0,       0,             0,  -(6*E*Iz)/L2,         0,  (12*E*Iz)/L3,             0,       0,            0, -(6*E*Iz)/L2},
		{       0,              0,  -(12*E*Iy)/L3,       0,   (6*E*Iy)/L2,             0,         0,             0,  (12*E*Iy)/L3,       0,  (6*E*Iy)/L2,            0},
		{       0,              0,              0, -G*Ik/L,             0,             0,         0,             0,             0,  G*Ik/L,            0,            0},
		{       0,              0,   -(6*E*Iy)/L2,       0,    (2*E*Iy)/L,             0,         0,             0,   (6*E*Iy)/L2,       0,   (4*E*Iy)/L,            0},
		{       0,    (6*E*Iz)/L2,              0,       0,             0,    (2*E*Iz)/L,         0,  -(6*E*Iz)/L2,             0,       0,            0,   (4*E*Iz)/L}
	};
#endif
	return K;
}

math::vector_t<double,2> ElemBEAM::stiffness(const math::vector<double>& property
											 , const math::vector<double>& material) const
{
	return stiffnessStatic(property,material,parameters);
}

math::vector_t<double,2> ElemBEAM::massStatic(const math::vector<double>& property
	 								 		 , const math::vector<double>& material
									 		 , const math::vector<double>& parameters)
{
	double rho = material[ElemBEAM::matl::rho];
	double A = property[ElemBEAM::prop::A];
	double Jp = property[ElemBEAM::prop::Jp];
    double L = parameters[ElemBEAM::prms::L];
    double L2 = L*L;
    double L3 = L2*L;

    
    math::vector_t<double,2> M = {
		{(A*L*rho)/3,                  0,                  0,            0,                  0,                 0, (A*L*rho)/6,                  0,                  0,            0,                 0,                  0},
		{          0,    (13*A*L*rho)/35,                  0,            0,                  0, (11*A*L2*rho)/210,           0,     (9*A*L*rho)/70,                  0,            0,                 0, -(13*A*L2*rho)/420},
		{          0,                  0,    (13*A*L*rho)/35,            0, -(11*A*L2*rho)/210,                 0,           0,                  0,     (9*A*L*rho)/70,            0, (13*A*L2*rho)/420,                  0},
		{          0,                  0,                  0, (Jp*L*rho)/3,                  0,                 0,           0,                  0,                  0, (Jp*L*rho)/6,                 0,                  0},
		{          0,                  0, -(11*A*L2*rho)/210,            0,     (A*L3*rho)/105,                 0,           0,                  0, -(13*A*L2*rho)/420,            0,   -(A*L3*rho)/140,                  0},
		{          0,  (11*A*L2*rho)/210,                  0,            0,                  0,    (A*L3*rho)/105,           0,  (13*A*L2*rho)/420,                  0,            0,                 0,    -(A*L3*rho)/140},
		{(A*L*rho)/6,                  0,                  0,            0,                  0,                 0, (A*L*rho)/3,                  0,                  0,            0,                 0,                  0},
		{          0,     (9*A*L*rho)/70,                  0,            0,                  0, (13*A*L2*rho)/420,           0,    (13*A*L*rho)/35,                  0,            0,                 0, -(11*A*L2*rho)/210},
		{          0,                  0,     (9*A*L*rho)/70,            0, -(13*A*L2*rho)/420,                 0,           0,                  0,    (13*A*L*rho)/35,            0, (11*A*L2*rho)/210,                  0},
		{          0,                  0,                  0, (Jp*L*rho)/6,                  0,                 0,           0,                  0,                  0, (Jp*L*rho)/3,                 0,                  0},
		{          0,                  0,  (13*A*L2*rho)/420,            0,    -(A*L3*rho)/140,                 0,           0,                  0,  (11*A*L2*rho)/210,            0,    (A*L3*rho)/105,                  0},
		{          0, -(13*A*L2*rho)/420,                  0,            0,                  0,   -(A*L3*rho)/140,           0, -(11*A*L2*rho)/210,                  0,            0,                 0,     (A*L3*rho)/105}
	};
	return M;
}

math::vector_t<double,2> ElemBEAM::mass(const math::vector<double>& property
											 , const math::vector<double>& material) const
{
	return massStatic(property,material,parameters);
}

void ElemBEAM::calc_parameters(const math::vector<Node>& nodes_info) {
	// R0 123 312 231
	auto t10 = nodes_info[nodes[1]].xyz - nodes_info[nodes[0]].xyz;
	double L = math::norm(t10); // length
	t10 /= L;
	auto t30 = prms_orientVec();
	auto t20 = math::cross(t30,t10);
	parameters = {
		// orient vector
		t30[0],t30[1],t30[2]

		// length
		,L

		// R0
		,t10[0],t20[0],t30[0]
		,t10[1],t20[1],t30[1]
		,t10[2],t20[2],t30[2]
		// ,t10[0],t10[1],t10[2]
		// ,t20[0],t20[1],t20[2]
		// ,t30[0],t30[1],t30[2]
	};
}

size_t ElemBEAM::nnodes() const {
	return 2;
}
size_t ElemBEAM::ndofs_node() const {
	return 6;
}

double ElemBEAM::prms_L() const {
	return parameters[ElemBEAM::prms::L];
}


block_t ElemBEAM::prms_R0() {
	return BaseElement::prms_R0<ElemBEAM>(parameters.begin());
}
const_block_t ElemBEAM::prms_R0() const {
	return BaseElement::prms_R0<ElemBEAM>(parameters.begin());
}

math::Slice<math::vector<double>::const_iterator
	       ,math::vector<double>::const_iterator> ElemBEAM::prms_orientVec() const {
	auto it = parameters.begin();
	math::Slice v(it + ElemBEAM::prms::vx, it + ElemBEAM::prms::vz+1);
	return v;
}

math::Slice<math::vector<double>::iterator
	       ,math::vector<double>::iterator> ElemBEAM::prms_orientVec() {
	auto it = parameters.begin();
	math::Slice v(it + ElemBEAM::prms::vx, it + ElemBEAM::prms::vz+1);
	return v;
}


ElemBEAMLD::ElemBEAMLD(size_t ID, size_t propID, size_t matlID
			   ,const math::vector<size_t>& nodes
			   ,const math::vector<double>& orient_vector)
		: ElemBEAMLD(ID,propID,matlID,nodes) {
	parameters = orient_vector / math::norm(orient_vector);
}

void ElemBEAMLD::tangentStiffness_innerLoad(math::vector_t<double,2>& K
											 , math::vector_t<double,1>& inner_load
										,const math::vector_t<double,1>& property
	 								    ,const math::vector_t<double,1>& material
                                        ,const math::vector_t<double,3>& Rsum
										,const math::vector_t<double,1>& q ) const
{
	/*
    Rsum[0] == R0a, Rsum[1] == R0b - total rotation tensors
    */
	auto itq = q.begin();
	math::Slice ua(itq+dof::uax,itq+dof::tax);
	math::Slice ub(itq+dof::ubx,itq+dof::tbx);
	math::Slice ta(itq+dof::tax,itq+dof::ubx);
	math::Slice tb(itq+dof::tbx,itq+dof::tbz+1);
	
	
#ifdef BEAM_RELEASES
	size_t release_shift1 = prms_shift1()*3
		 , release_shift2 = prms_shift2()*3;
	math::Slice Rsuma(Rsum[nodes[0]].begin() + release_shift1, Rsum[nodes[0]].begin() + release_shift1+3);
	math::Slice Rsumb(Rsum[nodes[1]].begin() + release_shift2, Rsum[nodes[1]].begin() + release_shift2+3);

	math::vector<double> theta1  = math::vector_invariant(Rsumb,Rsuma); // vector invariant of (R0b * RoaT) tensor
	double theta1_abs = 0.5 * std::asin(0.5*math::norm(theta1));
	
	theta1 *= -0.25 / math::detail::rotation_tensor_helper2(2.*theta1_abs);
	

    // middle element point total rotation tensor
	math::vector_t<double,2> Rm = math::dot(math::rotation_tensor(theta1),Rsuma);

#else
	math::vector<double> theta1  = math::vector_invariant(Rsum[nodes[1]],Rsum[nodes[0]]); // vector invariant of (R0b * RoaT) tensor

    double theta1_abs = 0.5 * std::asin(0.5*math::norm(theta1));
	
	theta1 *= -0.25 / math::detail::rotation_tensor_helper2(2.*theta1_abs);
	

    // middle element point total rotation tensor
	math::vector_t<double,2> Rm = math::dot(math::rotation_tensor(theta1),Rsum[nodes[0]]);
#endif
	
	// block of rotation tensor from 0 to actual position
	const_block_t R0 = prms_R0();
    math::vector_t<double,2> Rall_block = math::dot(Rm,R0); 
	// math::vector_t<double,2> RallT_block = math::transpose(Rall_block); 

	
	// Small displacement vector (phisical items same as in global dispalcement vector q)
	math::vector<double> z(q.size());
	auto itz = z.begin();
	math::Slice xia(  itz+dof::uax,itz+dof::uaz+1);
	math::Slice xib(  itz+dof::ubx,itz+dof::ubz+1);
	math::Slice betaa(itz+dof::tax,itz+dof::taz+1);
	math::Slice betab(itz+dof::tbx,itz+dof::tbz+1);
	
	double L = parameters[ElemBEAM::prms::L];
	auto ex0 = prms_basis0_ex0();
	xia = (ua-ub)/2. + (math::dot(Rm,ex0) - ex0)*(L/2.);
	xib -= xia; // xib was zero vector
	
	// betaa = ta - theta1;
	math::sub(ta.begin(),theta1.begin(),betaa.begin(),betaa.end());
	// betab = tb + theta1;
	math::sum(tb.begin(),theta1.begin(),betab.begin(),betab.end());
	
	// stiffness matrix
	auto K0 = ElemBEAM::stiffnessStatic(property,material,parameters); // inital (for linear element)
	auto K2 = math::zeros<double>(K0); // in 2nd state 
	
	// calculate K2 = Rall * K0 * RallT
	static const size_t blocksize = BaseNode::DIM;
	
	block_t K0block(blocksize);
	block_t K2block(blocksize); // upper triangular blocks
	block_t Kblock_sub(blocksize); // lower triangular blocks
								// // which equal to transposed upper triangular blocks
	math::vector_t<double,2> block_sub = math::zeros<double>(blocksize,blocksize);
	
	
	// loop over n blocks of K matrix
	size_t n_blocks = ndofs()/blocksize;
	for (size_t i = 0; i < n_blocks; ++i) {
		for (size_t j = i; j < n_blocks; ++j) {
			
			math::matrix_block_set(K0,i,j,blocksize,blocksize,K0block);
			math::matrix_block_set(K2,i,j,blocksize,blocksize,K2block);
			
			// block calculation
			math::dotT(K0block, Rall_block, block_sub);
			math::dot(Rall_block, block_sub, K2block);
			// K2block = math::dot(Rall_block, math::dot(K0block, RallT_block) ); // Kblock_ij = Rall_block * K0block_ij * RallT_block
			
			if (i == j) continue;
			math::matrix_block_set(K2,j,i,blocksize,blocksize,Kblock_sub);
			math::transpose(K2block,Kblock_sub);
		}
	}
	
	// calculate HTK2 = HT * K2 and K1 = HTK2*H
	/*	HTK2 = (here K == K2) 
	 [  1/2(K1_1 - K3_1), 1/2(K1_2 - K3_2), 1/2(K1_3 - K3_3), 1/2(K1_4 - K3_4)
		K2_1, 				K2_2, 			K2_3, 				K2_4
		1/2(K3_1 - K1_1), 1/2(K3_2 - K1_2), 1/2(K3_3 - K1_3), 1/2(K3_4 - K1_4)
		K4_1, 				K4_2, 			K4_3, 				K4_4			]
	*/

#if 0
	K = K2;
	block_t HTK2block1 = std::move(K0block), HTK2block2 = std::move(K2block);
	block_t tmpblock1 = std::move(Kblock_sub);
	block_t tmpblock2(blocksize);
	for (size_t k = 0; k < n_blocks; ++k) {
		math::matrix_block_set(K2,0,k,blocksize,blocksize,tmpblock1);
		math::matrix_block_set(K2,2,k,blocksize,blocksize,tmpblock2);
		math::matrix_block_set(K,0,k,blocksize,blocksize,HTK2block1);
		HTK2block1 = 0.5*(tmpblock1-tmpblock2);
		
		math::matrix_block_set(K,2,k,blocksize,blocksize,HTK2block2);
		HTK2block2 = -1*HTK2block1; // TODO
	}
#else
	
    static const math::vector_t<double,2> H = {
        { 0.5,  0,   0,   0,   0,   0, -0.5,  0,   0,   0,   0,   0 },
        {  0,  0.5,  0,   0,   0,   0,   0, -0.5,  0,   0,   0,   0 },
        {  0,   0,  0.5,  0,   0,   0,   0,   0, -0.5,  0,   0,   0 },
        {  0,   0,   0,  1.0,  0,   0,   0,   0,   0,   0,   0,   0 },
        {  0,   0,   0,   0,  1.0,  0,   0,   0,   0,   0,   0,   0 },
        {  0,   0,   0,   0,   0,  1.0,  0,   0,   0,   0,   0,   0 },
        {-0.5,  0,   0,   0,   0,   0,  0.5,  0,   0,   0,   0,   0 },
        {  0, -0.5,  0,   0,   0,   0,   0,  0.5,  0,   0,   0,   0 },
        {  0,   0, -0.5,  0,   0,   0,   0,   0,  0.5,  0,   0,   0 },
        {  0,   0,   0,   0,   0,   0,   0,   0,   0,  1.0,  0,   0 },
        {  0,   0,   0,   0,   0,   0,   0,   0,   0,   0,  1.0,  0 },
        {  0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,  1.0}
    }; // H == HT
	
	/* Fill K with 0 */
	{
		auto itK = K.begin();
		auto itK_end = itK + ndofs();
		while (itK < itK_end) {
			auto itK_col = itK->begin();
			auto itK_col_end = itK_col + ndofs();
			while (itK_col < itK_col_end) {
				*itK_col = 0.0;
				++itK_col;
			}
			++itK;
		}
	}
	math::dot(H,K2,K);

#endif
	// // 2d variant
	// block_t HTK2block1 = std::move(K0block), HTK2block2 = std::move(K2block);
	// block_t tmpblock1 = std::move(Kblock_sub);
	// block_t tmpblock2(blocksize);
	// // loop over block columns
	// for (size_t k = 0; k < n_blocks; ++k) {
	// 	// 1st row
	// 	math::matrix_block_set(K2,0,k,blocksize,blocksize,tmpblock1);
	// 	math::matrix_block_set(K2,2,k,blocksize,blocksize,tmpblock2);
	// 	math::matrix_block_set(K ,0,k,blocksize,blocksize,HTK2block1);
	// 	HTK2block1 = 0.5*(tmpblock1-tmpblock2);

	// 	// 2nd row
	// 	math::matrix_block_set(K2,1,k,blocksize,blocksize,tmpblock1);
	// 	math::matrix_block_set(K ,1,k,blocksize,blocksize,HTK2block1);
	// 	HTK2block1 = tmpblock1;

	// 	// 3rd row
	// 	math::matrix_block_set(K ,2,k,blocksize,blocksize,tmpblock1);
	// 	tmpblock1 = -1*HTK2block1;

	// 	// 4th row
	// 	math::matrix_block_set(K2,3,k,blocksize,blocksize,tmpblock1);
	// 	math::matrix_block_set(K ,3,k,blocksize,blocksize,HTK2block1);
	// 	HTK2block1 = tmpblock1;
	// }
	// std::cout << "\nK = \n" << K << '\n';
	

	// std::cout << "HTK2 - math::dot(H,K2) = \n" << K - math::dot(H,K2) << std::endl;
	// vector of internal loads
	math::dot(K,z,inner_load);
	
	// stiffness matrix K1
	// auto K = math::zeros<double>(K0);
	
#if 0
	block_t Kblock1 = std::move(HTK2block1), Kblock2 = std::move(HTK2block2);
	
	math::matrix_block_set(K,0,1,blocksize,blocksize,Kblock1);
	math::matrix_block_set(K,1,0,blocksize,blocksize,Kblock2);
	Kblock2 = math::transpose(Kblock1);
	
	math::matrix_block_set(K,2,1,blocksize,blocksize,Kblock1);
	math::matrix_block_set(K,1,2,blocksize,blocksize,Kblock2);
	Kblock2 = math::transpose(Kblock1);
	
	math::matrix_block_set(K,0,3,blocksize,blocksize,Kblock1);
	math::matrix_block_set(K,3,0,blocksize,blocksize,Kblock2);
	Kblock2 = math::transpose(Kblock1);
	
	math::matrix_block_set(K,2,3,blocksize,blocksize,Kblock1);
	math::matrix_block_set(K,3,2,blocksize,blocksize,Kblock2);
	Kblock2 = math::transpose(Kblock1);
	
	math::matrix_block_set(K,0,0,blocksize,blocksize,Kblock1);
	math::matrix_block_set(K,2,2,blocksize,blocksize,Kblock2);
	Kblock1 = (Kblock1 + Kblock2)/2.;
	Kblock2 = Kblock1;
	
	math::matrix_block_set(K,0,2,blocksize,blocksize,Kblock1);
	math::matrix_block_set(K,2,0,blocksize,blocksize,Kblock2);
	Kblock1 = (Kblock1 + Kblock2)/2.;
	Kblock2 = Kblock1;
#else
	K = math::dot(K,H);
#endif

	// std::cout << "HTK2H - math::dot(H,K2,H) = \n" << K - math::dot(math::dot(H,K2),H) << std::endl;
	
	
	
	// // 2nd part of the stiffness matrix K2block
	block_t Kblock1 = std::move(K0block), Kblock2 = std::move(K2block);
	auto itloadInt = inner_load.begin();
	math::Slice Fa(itloadInt+dof::uax, itloadInt+dof::tax);
	math::Slice Fb(itloadInt+dof::ubx, itloadInt+dof::tbx);
	auto S1 = math::skew_tensor(Fa*0.25);
	auto S2 = math::skew_tensor(Fb*0.25);
	
	for (size_t k = 0; k < 2; ++k) {
		math::matrix_block_set(K,0,2*k+1,blocksize,blocksize,Kblock1);
		math::matrix_block_set(K,2*k+1,0,blocksize,blocksize,Kblock2);
		Kblock1 -= S1;
		Kblock2 += S1;
		
		math::matrix_block_set(K,2,2*k+1,blocksize,blocksize,Kblock1);
		math::matrix_block_set(K,2*k+1,2,blocksize,blocksize,Kblock2);
		Kblock1 -= S2;
		Kblock2 += S2;
	}
	
	// math::vector_t<double,2> S = math::zeros<double>(12,12);
	// for (size_t i = 0; i < 3; ++i) {
	// 	for (size_t j = 0; j < 3; ++j) {
	// 		S[i][j+3]   = -S1[i][j];
	// 		S[i][j+9]   = -S1[i][j];
	// 		S[i+6][j+3] = -S2[i][j];
	// 		S[i+6][j+9] = -S2[i][j];
	// 	}
	// }
	// K += (S + math::transpose(S));

	// std::cout 	
	// // << "det(K) = " << math::det_<math::vector_t<double,2>>(K)
	// 			// << ", det(K0) = " << math::det_<math::vector_t<double,2>>(K0)
	// 			<< "K: det(Rm) = " << math::det_<math::vector_t<double,2>>(Rm)
	// 			<< ", det(R0) = " << math::det_<const_block_t>(R0)
	// 			<< '\n';
}

/* Calculate tangent mass matrix `M` and inertia load vector `inert_load`*/
void ElemBEAMLD::tangentMass_inertiaLoad(math::vector_t<double,2>& M
									  , math::vector_t<double,1>& inert_load
								, const math::vector_t<double,1>& property
	 						    , const math::vector_t<double,1>& material
                                , const math::vector_t<double,3>& Rsum
								, const math::vector_t<double,1>& q				/* displacements */
								, const math::vector_t<double,1>& dqdt			/* velocities */
								, const math::vector_t<double,1>& d2qdt2 		/* accelerations */) const
{
	/* Nodal displacements */
	auto itq = q.begin();
	math::Slice ua(itq+dof::uax, itq+dof::uaz+1);
	math::Slice ub(itq+dof::ubx, itq+dof::ubz+1);
	math::Slice ta(itq+dof::tax, itq+dof::taz+1);
	math::Slice tb(itq+dof::tbx, itq+dof::tbz+1);

	itq = dqdt.begin();
	math::Slice duadt(itq+dof::uax, itq+dof::uaz+1);
	math::Slice dubdt(itq+dof::ubx, itq+dof::ubz+1);
	math::Slice dtadt(itq+dof::tax, itq+dof::taz+1);
	math::Slice dtbdt(itq+dof::tbx, itq+dof::tbz+1);

	
	/* Rotation and Zhilin tensors and their derivative of Euler vecto r*/
	math::vector_t<double,2> La = math::rotation_tensor(ta);
	math::vector_t<double,2> Lb = math::rotation_tensor(tb);
	math::vector_t<double,2> Ba = math::zhilin_tensor(ta);
	math::vector_t<double,2> Bb = math::zhilin_tensor(tb);
	math::vector_t<double,3> dLdta = math::zeros<double>(BaseNode::DIM,BaseNode::DIM,BaseNode::DIM);
	math::vector_t<double,3> dLdtb = math::zeros<double>(BaseNode::DIM,BaseNode::DIM,BaseNode::DIM);
	math::vector_t<double,3> dBdta = math::zeros<double>(BaseNode::DIM,BaseNode::DIM,BaseNode::DIM);
	math::vector_t<double,3> dBdtb = math::zeros<double>(BaseNode::DIM,BaseNode::DIM,BaseNode::DIM);
	math::vector_t<double,3> dLTdta = math::zeros<double>(BaseNode::DIM,BaseNode::DIM,BaseNode::DIM);
	math::vector_t<double,3> dLTdtb = math::zeros<double>(BaseNode::DIM,BaseNode::DIM,BaseNode::DIM);
	math::vector_t<double,3> dBTdta = math::zeros<double>(BaseNode::DIM,BaseNode::DIM,BaseNode::DIM);
	math::vector_t<double,3> dBTdtb = math::zeros<double>(BaseNode::DIM,BaseNode::DIM,BaseNode::DIM);
	
	math::rotation_tensor_diff(ta,dLdta,La,Ba);
	math::rotation_tensor_diff(tb,dLdtb,Lb,Bb);
	math::zhilin_tensor_diff(ta,dBdta,Ba,dLdta);
	math::zhilin_tensor_diff(tb,dBdtb,Bb,dLdtb);

	// components of dNdt
	math::vector_t<double,2> dtadt_dLdta = math::dot(dtadt,dLdta);
	math::vector_t<double,2> dtadt_dBdta = math::dot(dtadt,dBdta);
	math::vector_t<double,2> dtbdt_dLdtb = math::dot(dtbdt,dLdtb);
	math::vector_t<double,2> dtbdt_dBdtb = math::dot(dtbdt,dBdtb);
	
	
	/* Total rotetion tensors */
	static const size_t blocksize = BaseNode::DIM;
	/* Calculate rotation tensor R1 = [R0a*R0...] */
	math::vector_t<double,2> R1a = math::zeros<double>(blocksize,blocksize);
	math::vector_t<double,2> R1b = math::zeros<double>(blocksize,blocksize);
	const_block_t R0 = prms_R0();

#ifdef BEAM_RELEASES
	size_t release_shift1 = prms_shift1()*3
		 , release_shift2 = prms_shift2()*3;
	math::Slice Rsuma(Rsum[nodes[0]].begin() + release_shift1, Rsum[nodes[0]].begin() + release_shift1+3);
	math::Slice Rsumb(Rsum[nodes[1]].begin() + release_shift2, Rsum[nodes[1]].begin() + release_shift2+3);
	
	math::dot(Rsuma,R0,R1a);
	math::dot(Rsumb,R0,R1b);
#else
	math::dot(Rsum[nodes[0]],R0,R1a);
	math::dot(Rsum[nodes[1]],R0,R1b);
#endif
	


	/* Calculate simultaneously:
		o Mstat = R1   * M0    * R1T 
	   	o M     = N    * Mstat * NT 
	   	o dMdt  = dNdt * Mstat * NT  + N * Mstat * dNTdt
	*/
	// Inital mass matrix - for linear element
	auto M0 = ElemBEAM::massStatic(property,material,parameters);
	
	// constant for current step mass matrix
	auto Mstat = math::zeros<double>(M0);
	// Derivative mass matrix by time
	auto dMdt = math::zeros<double>(M0);
	
	/* Fill M with 0.0 */
	{
		auto itM = M.begin();
		auto itM_end = itM + ndofs();
		while (itM < itM_end) {
			auto itM_col = itM->begin();
			auto itM_col_end = itM_col + ndofs();
			while (itM_col < itM_col_end) {
				*itM_col = 0.0;
				++itM_col;
			}
			++itM;
		}
	}
	// M like block definition
	block_t block_M0(blocksize)		// M0 block
	       ,block_Mstat(blocksize)	// Mstat block
	       ,block_M(blocksize)		// M block
	       ,block_dMdt(blocksize)	// dMdt block
		   ,block_sub2(blocksize);	// sub M like block
	math::vector_t<double,2> block_sub = math::zeros<double>(blocksize,blocksize); 	// sub M like block

	// pointers on transformation matrices
	math::vector_t<double,2> *pR1_left  // pointers on R1a 
							,*pR1_right	// or R1b
							,*pN_left
							,*pN_right
							,*pdNdt_left
							,*pdNdt_right;
	
	// loop over n blocks of M like matrices
	size_t n_blocks = ndofs()/blocksize;  // 4
	for (size_t i = 0; i < n_blocks; ++i) {
		switch (i)
		{
			case 0:
				pR1_left   = &R1a;
				pN_left    = &La;
				pdNdt_left = &dtadt_dLdta;
				break;
			case 1:
				pR1_left   = &R1a;
				pN_left    = &Ba;
				pdNdt_left = &dtadt_dBdta;
				break;
			case 2:
				pR1_left   = &R1b;
				pN_left    = &Lb;
				pdNdt_left = &dtbdt_dLdtb;
				break;
			case 3:
				pR1_left   = &R1b;
				pN_left    = &Bb;
				pdNdt_left = &dtbdt_dBdtb;
				break;
			default:
				assert(false);
		}
		
		for (size_t j = i; j < n_blocks; ++j) {
			
			// block declaretion
			math::matrix_block_set(M0   ,i,j,blocksize,blocksize,block_M0);
			math::matrix_block_set(Mstat,i,j,blocksize,blocksize,block_Mstat);
			math::matrix_block_set(M    ,i,j,blocksize,blocksize,block_M);
			math::matrix_block_set(dMdt ,i,j,blocksize,blocksize,block_dMdt);

			// block calculation
			switch (j)
			{
				case 0:
					pR1_right   = &R1a;
					pN_right    = &La;
					pdNdt_right = &dtadt_dLdta;
					break;
				case 1:
					pR1_right   = &R1a;
					pN_right    = &Ba;
					pdNdt_right = &dtadt_dBdta;
					break;
				case 2:
					pR1_right   = &R1b;
					pN_right    = &Lb;
					pdNdt_right = &dtbdt_dLdtb;
					break;
				case 3:
					pR1_right   = &R1b;
					pN_right    = &Bb;
					pdNdt_right = &dtbdt_dBdtb;
					break;
				default:
					assert(false);
		}
			// Mstat = R1   * M0    * R1T
			math::dotT(block_M0, 	*pR1_right, block_sub);
			math::dot(*pR1_left, 	block_sub, block_Mstat);
			// M     = N    * Mstat * NT
			math::dotT(block_Mstat, *pN_right, 	block_sub); 	
			math::dot(*pN_left,  	block_sub, block_M);
			// dMdt  = dNdt * Mstat * NT +
			math::dotT(block_Mstat, *pN_right, 	block_sub);
			math::dot(*pdNdt_left,  block_sub, block_dMdt);
			//          + N * Mstat * dNTdt
			math::dotT(block_Mstat, *pdNdt_right, 	block_sub);
			math::dot(*pN_left,  	block_sub, 	block_dMdt);   	
			
			
			// if not diagonal ij element - get ji element by transpose ij
			if (i == j) continue;
			math::matrix_block_set(Mstat,j,i,blocksize,blocksize,block_sub2);
			math::transpose(block_Mstat,block_sub2);
			
			math::matrix_block_set(M,j,i,blocksize,blocksize,block_sub2);
			math::transpose(block_M,block_sub2);

			math::matrix_block_set(dMdt,j,i,blocksize,blocksize,block_sub2);
			math::transpose(block_dMdt,block_sub2);
		}
	}

	/* Calculate pseudo velocities vector */
	math::vector<double> pseudo_v(ndofs());
	math::vector<double> pseudo_v0(ndofs());
	auto itpseudo_v = pseudo_v0.begin();
	math::Slice pseudo_va(itpseudo_v+dof::uax, itpseudo_v+dof::uaz+1);
	math::Slice pseudo_vb(itpseudo_v+dof::ubx, itpseudo_v+dof::ubz+1);
	math::Slice pseudo_wa(itpseudo_v+dof::tax, itpseudo_v+dof::taz+1);
	math::Slice pseudo_wb(itpseudo_v+dof::tbx, itpseudo_v+dof::tbz+1);

	math::dot(duadt,La,pseudo_va);
	math::dot(dtadt,Ba,pseudo_wa);
	math::dot(dubdt,Lb,pseudo_vb);
	math::dot(dtbdt,Bb,pseudo_wb);

	math::dot(Mstat,pseudo_v0,pseudo_v);

	itpseudo_v = pseudo_v.begin();
	pseudo_va.new_slice(itpseudo_v+dof::uax, itpseudo_v+dof::uaz+1);
	pseudo_vb.new_slice(itpseudo_v+dof::ubx, itpseudo_v+dof::ubz+1);
	pseudo_wa.new_slice(itpseudo_v+dof::tax, itpseudo_v+dof::taz+1);
	pseudo_wb.new_slice(itpseudo_v+dof::tbx, itpseudo_v+dof::tbz+1);


	/* Calculate componsents of Finert */	
	math::fill(inert_load.begin(),inert_load.begin()+ndofs(),0.0);
	math::dot(M,d2qdt2,inert_load);
#if 0
	math::dot(dMdt,dqdt,inert_load); // M*d2qdt2 + dMdt*dqdt
	
	/* Calculate componsents of dTdq */
	// clculate dLTdv and dBTdv and place reults to dLdv and dBdv
	math::rotation_tensor_transpose_diff(ta,dLTdta,La,Ba);
	math::rotation_tensor_transpose_diff(tb,dLTdtb,Lb,Bb);
	math::zhilin_tensor_transpose_diff(ta,dBTdta,Ba,dLdta);
	math::zhilin_tensor_transpose_diff(tb,dBTdtb,Bb,dLdtb);
	
	// !!!!!!!
	// change preudo_v ->   -preudo_v
	pseudo_v *= -1.0;
	// !!!!!!!
	
	auto itfinert = inert_load.begin();
	math::Slice inert_load_ta(itfinert+dof::tax,itfinert+dof::taz+1);
	math::Slice inert_load_tb(itfinert+dof::tbx,itfinert+dof::tbz+1);

	// dTdta
	math::fill(block_sub.begin(),block_sub.end(),0.0);
	math::dot(dLTdta,duadt,block_sub);
	math::dot(block_sub,pseudo_va,inert_load_ta);
	
	math::fill(block_sub.begin(),block_sub.end(),0.0);
	math::dot(dBTdta,dtadt,block_sub);
	math::dot(block_sub,pseudo_wa,inert_load_ta);

	// dTdtb
	math::fill(block_sub.begin(),block_sub.end(),0.0);
	math::dot(dLTdtb,dubdt,block_sub);
	math::dot(block_sub,pseudo_vb,inert_load_tb);
	
	math::fill(block_sub.begin(),block_sub.end(),0.0);
	math::dot(dBdtb,dtbdt,block_sub);
	math::dot(block_sub,pseudo_wb,inert_load_tb);
#endif
	// std::cout 	
	// // << ", det(M) = " << math::det_<math::vector_t<double,2>>(M)
	// 			// << ", det(M0) = " << math::det_<math::vector_t<double,2>>(M0)
	// 			<< "M: det(R1a) = " << math::det_<math::vector_t<double,2>>(R1a)
	// 			<< ", det(R1b) = " << math::det_<math::vector_t<double,2>>(R1b)
	// 			<< ", det(R0) = " << math::det_<const_block_t>(R0)
	// 			<< '\n';
	
	// std::cout << ", |M| = " << math::norm(M) << ", |dMdt| = " << math::norm(dMdt) << '\n';
	
	/*
	abc * d = ab(c*d)
	(abc * d)^T = d * (abc)^T = d * cda = (d*c)da =(d*c)(ad)^T 

	(abc*d)*e = ((abc*d)*e)^T = e*(abc*d)^T 
	*/ 


}

/* Calculate tangent mass matrix `M`, mass- and gyro- effects `G`  and inertia load vector `inert_load`*/
void ElemBEAMLD::tangentMassGyro_inertiaLoad(math::vector_t<double,2>& M
							 	   , math::vector_t<double,2>& G
							 	   , math::vector_t<double,1>& inert_load
							 , const math::vector_t<double,1>& property
	 						 , const math::vector_t<double,1>& material
                        	 , const math::vector_t<double,3>& Rsum
							 , const math::vector_t<double,1>& q
							 , const math::vector_t<double,1>& dqdt
							 , const math::vector_t<double,1>& d2qdt2 ) const
{
	/* Nodal displacements */
	auto itq = q.begin();
	math::Slice ua(itq+dof::uax, itq+dof::uaz+1);
	math::Slice ub(itq+dof::ubx, itq+dof::ubz+1);
	math::Slice ta(itq+dof::tax, itq+dof::taz+1);
	math::Slice tb(itq+dof::tbx, itq+dof::tbz+1);

	itq = dqdt.begin();
	math::Slice duadt(itq+dof::uax, itq+dof::uaz+1);
	math::Slice dubdt(itq+dof::ubx, itq+dof::ubz+1);
	math::Slice dtadt(itq+dof::tax, itq+dof::taz+1);
	math::Slice dtbdt(itq+dof::tbx, itq+dof::tbz+1);

	
	/* Rotation and Zhilin tensors and their derivative of Euler vecto r*/
	math::vector_t<double,2> La = math::rotation_tensor(ta);
	math::vector_t<double,2> Lb = math::rotation_tensor(tb);
	math::vector_t<double,2> Ba = math::zhilin_tensor(ta);
	math::vector_t<double,2> Bb = math::zhilin_tensor(tb);
	math::vector_t<double,3> dLdta = math::zeros<double>(BaseNode::DIM,BaseNode::DIM,BaseNode::DIM);
	math::vector_t<double,3> dLdtb = math::zeros<double>(BaseNode::DIM,BaseNode::DIM,BaseNode::DIM);
	math::vector_t<double,3> dBdta = math::zeros<double>(BaseNode::DIM,BaseNode::DIM,BaseNode::DIM);
	math::vector_t<double,3> dBdtb = math::zeros<double>(BaseNode::DIM,BaseNode::DIM,BaseNode::DIM);
	math::vector_t<double,3> dLTdta = math::zeros<double>(BaseNode::DIM,BaseNode::DIM,BaseNode::DIM);
	math::vector_t<double,3> dLTdtb = math::zeros<double>(BaseNode::DIM,BaseNode::DIM,BaseNode::DIM);
	math::vector_t<double,3> dBTdta = math::zeros<double>(BaseNode::DIM,BaseNode::DIM,BaseNode::DIM);
	math::vector_t<double,3> dBTdtb = math::zeros<double>(BaseNode::DIM,BaseNode::DIM,BaseNode::DIM);
	
	math::rotation_tensor_diff(ta,dLdta,La,Ba);
	math::rotation_tensor_diff(tb,dLdtb,Lb,Bb);
	math::zhilin_tensor_diff(ta,dBdta,Ba,dLdta);
	math::zhilin_tensor_diff(tb,dBdtb,Bb,dLdtb);

	// components of dNdt
	math::vector_t<double,2> dtadt_dLdta = math::dot(dtadt,dLdta);
	math::vector_t<double,2> dtadt_dBdta = math::dot(dtadt,dBdta);
	math::vector_t<double,2> dtbdt_dLdtb = math::dot(dtbdt,dLdtb);
	math::vector_t<double,2> dtbdt_dBdtb = math::dot(dtbdt,dBdtb);
	
	
	/* Total rotetion tensors */
	static const size_t blocksize = BaseNode::DIM;
	/* Calculate rotation tensor R1 = [R0a*R0...] */
	math::vector_t<double,2> R1a = math::zeros<double>(blocksize,blocksize);
	math::vector_t<double,2> R1b = math::zeros<double>(blocksize,blocksize);
	const_block_t R0 = prms_R0();

#ifdef BEAM_RELEASES
	size_t release_shift1 = prms_shift1()*3
		 , release_shift2 = prms_shift2()*3;
	math::Slice Rsuma(Rsum[nodes[0]].begin() + release_shift1, Rsum[nodes[0]].begin() + release_shift1+3);
	math::Slice Rsumb(Rsum[nodes[1]].begin() + release_shift2, Rsum[nodes[1]].begin() + release_shift2+3);
	
	math::dot(Rsuma,R0,R1a);
	math::dot(Rsumb,R0,R1b);
#else
	math::dot(Rsum[nodes[0]],R0,R1a);
	math::dot(Rsum[nodes[1]],R0,R1b);
#endif
	


	/* Calculate simultaneously:
		o Mstat = R1   * M0    * R1T 
	   	o M     = N    * Mstat * NT 
	   	o dMdt  = dNdt * Mstat * NT  + N * Mstat * dNTdt
	*/
	// Inital mass matrix - for linear element
	auto M0 = ElemBEAM::massStatic(property,material,parameters);
	
	// constant for current step mass matrix
	auto Mstat = math::zeros<double>(M0);
	auto MstatNT = math::zeros<double>(M0); // Mstat * NT
	auto Wsub = math::zeros<double>(M0);
	// Derivative mass matrix by time
	auto dMdt = math::zeros<double>(M0);
	
	/* Fill M and G with 0.0 */
	{
		auto itM = M.begin();
		auto itM_end = itM + ndofs();
		auto itG = G.begin();
		auto itFinner = inert_load.begin();
		
		while (itM < itM_end) {
			auto itM_col = itM->begin();
			auto itM_col_end = itM_col + ndofs();
			auto itG_col = itG->begin();
			while (itM_col < itM_col_end) {
				*itM_col = 0.0;
				*itG_col = 0.0;
				++itM_col;
				++itG_col;
			}
			++itM;
			++itG;
			
			
			*itFinner = 0.0;
			++itFinner;
		}
	}
	// M like block definition
	block_t block_M0(blocksize)		// M0 block
	       ,block_Mstat(blocksize)	// Mstat block
		   ,block_MstatNT(blocksize)// Mstat*NT block
	       ,block_M(blocksize)		// M block
	       ,block_dMdt(blocksize)	// dMdt block
		   ,block_G(blocksize)		// G block
		   ,block_sub2(blocksize);	// sub M like block
	math::vector_t<double,2> block_sub = math::zeros<double>(blocksize,blocksize); 	// sub M like block

	// pointers on transformation matrices
	math::vector_t<double,2> *pR1_left  // pointers on R1a 
							,*pR1_right	// or R1b
							,*pN_left
							,*pN_right
							,*pdNdt_left
							,*pdNdt_right;
	
	// loop over n blocks of M like matrices
	size_t n_blocks = ndofs()/blocksize;  // 4
	for (size_t i = 0; i < n_blocks; ++i) {
		switch (i)
		{
			case 0:
				pR1_left   = &R1a;
				pN_left    = &La;
				pdNdt_left = &dtadt_dLdta;
				break;
			case 1:
				pR1_left   = &R1a;
				pN_left    = &Ba;
				pdNdt_left = &dtadt_dBdta;
				break;
			case 2:
				pR1_left   = &R1b;
				pN_left    = &Lb;
				pdNdt_left = &dtbdt_dLdtb;
				break;
			case 3:
				pR1_left   = &R1b;
				pN_left    = &Bb;
				pdNdt_left = &dtbdt_dBdtb;
				break;
			default:
				assert(false);
		}
		
		for (size_t j = i; j < n_blocks; ++j) {
			
			// block declaretion
			math::matrix_block_set(M0   	,i,j,blocksize,blocksize,block_M0);
			math::matrix_block_set(Mstat	,i,j,blocksize,blocksize,block_Mstat);
			math::matrix_block_set(MstatNT	,i,j,blocksize,blocksize,block_MstatNT);
			math::matrix_block_set(M    	,i,j,blocksize,blocksize,block_M);
			math::matrix_block_set(dMdt 	,i,j,blocksize,blocksize,block_dMdt);
			
			// block calculation
			switch (j)
			{
				case 0:
					pR1_right   = &R1a;
					pN_right    = &La;
					pdNdt_right = &dtadt_dLdta;
					break;
				case 1:
					pR1_right   = &R1a;
					pN_right    = &Ba;
					pdNdt_right = &dtadt_dBdta;
					break;
				case 2:
					pR1_right   = &R1b;
					pN_right    = &Lb;
					pdNdt_right = &dtbdt_dLdtb;
					break;
				case 3:
					pR1_right   = &R1b;
					pN_right    = &Bb;
					pdNdt_right = &dtbdt_dBdtb;
					break;
				default:
					assert(false);
		}
			// Mstat = R1   * M0    * R1T
			math::dotT(block_M0, 	*pR1_right, block_sub);
			math::dot(*pR1_left, 	block_sub, block_Mstat);
			// M     = N    * Mstat * NT
			math::dotT(block_Mstat, *pN_right, 	block_MstatNT);	
			math::dot(*pN_left,  	block_MstatNT, block_M);
			// dMdt  = dNdt * Mstat * NT +
			math::dot(*pdNdt_left,  block_MstatNT, block_dMdt);
			//          + N * Mstat * dNTdt
			math::dotT(block_Mstat, *pdNdt_right, 	block_sub);
			math::dot(*pN_left,  	block_sub, 	block_dMdt);   	
			
			
			// if not diagonal ij element - get ji element by transpose ij
			if (i == j);
			math::matrix_block_set(Mstat,j,i,blocksize,blocksize,block_sub2);
			math::transpose(block_Mstat,block_sub2);
			
			math::matrix_block_set(Mstat,j,i,blocksize,blocksize,block_Mstat);
			math::dotT(block_sub2,*pN_right,block_Mstat);

			math::matrix_block_set(M,j,i,blocksize,blocksize,block_sub2);
			math::transpose(block_M,block_sub2);

			math::matrix_block_set(dMdt,j,i,blocksize,blocksize,block_sub2);
			math::transpose(block_dMdt,block_sub2);

		}
	}

	/* Calculate pseudo velocities vector */
	math::vector<double> pseudo_v(ndofs());
	math::vector<double> pseudo_v0(ndofs());
	auto itpseudo_v = pseudo_v0.begin();
	math::Slice pseudo_va(itpseudo_v+dof::uax, itpseudo_v+dof::uaz+1);
	math::Slice pseudo_vb(itpseudo_v+dof::ubx, itpseudo_v+dof::ubz+1);
	math::Slice pseudo_wa(itpseudo_v+dof::tax, itpseudo_v+dof::taz+1);
	math::Slice pseudo_wb(itpseudo_v+dof::tbx, itpseudo_v+dof::tbz+1);

	math::dot(duadt,La,pseudo_va);
	math::dot(dtadt,Ba,pseudo_wa);
	math::dot(dubdt,Lb,pseudo_vb);
	math::dot(dtbdt,Bb,pseudo_wb);

	math::dot(Mstat,pseudo_v0,pseudo_v);

	itpseudo_v = pseudo_v.begin();
	pseudo_va.new_slice(itpseudo_v+dof::uax, itpseudo_v+dof::uaz+1);
	pseudo_vb.new_slice(itpseudo_v+dof::ubx, itpseudo_v+dof::ubz+1);
	pseudo_wa.new_slice(itpseudo_v+dof::tax, itpseudo_v+dof::taz+1);
	pseudo_wb.new_slice(itpseudo_v+dof::tbx, itpseudo_v+dof::tbz+1);


	/* Calculate componsents of Finert */	
	math::dot(M,d2qdt2,inert_load);	// M*d2qdt2

	math::dot(dMdt,dqdt,inert_load); // M*d2qdt2 + dMdt*dqdt
	
	/* Calculate componsents of dTdq */
	// clculate dLTdv and dBTdv and place reults to dLdv and dBdv
	math::rotation_tensor_transpose_diff(ta,dLTdta,La,Ba);
	math::rotation_tensor_transpose_diff(tb,dLTdtb,Lb,Bb);
	math::zhilin_tensor_transpose_diff(ta,dBTdta,Ba,dLdta);
	math::zhilin_tensor_transpose_diff(tb,dBTdtb,Bb,dLdtb);
	
	// !!!!!!!
	// change preudo_v   ->  -preudo_v
	pseudo_v *= -1.0;
	// !!!!!!!
	
	auto itfinert = inert_load.begin();
	math::Slice inert_load_ta(itfinert+dof::tax,itfinert+dof::taz+1);
	math::Slice inert_load_tb(itfinert+dof::tbx,itfinert+dof::tbz+1);

	/* dTdta */
	math::matrix_block_set(G,1,0,blocksize,blocksize,block_G);
	// inert_load_ta = (dLdta * pseudo_va) * duadt  + [...]
	math::dot(dLTdta,duadt,block_G);
	math::dot(block_G,pseudo_va,inert_load_ta);
	
	math::matrix_block_set(G,1,1,blocksize,blocksize,block_G);
	// [...] = (dBdta * pseudo_wa) * dtadt  
	math::dot(dBTdta,dtadt,block_G);
	math::dot(block_G,pseudo_wa,inert_load_ta);

	/* dTdtb */
	math::matrix_block_set(G,3,2,blocksize,blocksize,block_G);
	// inert_load_tb = (dLdtb * pseudo_vb) * dubdt  + [...]
	math::dot(dLTdtb,dubdt,block_G);
	math::dot(block_G,pseudo_vb,inert_load_tb);
	
	math::matrix_block_set(G,3,3,blocksize,blocksize,block_G);
	// [...] = (dBdtb * pseudo_wb) * dtbdt  
	math::dot(dBTdtb,dtbdt,block_G);
	math::dot(block_G,pseudo_wb,inert_load_tb);


	/* Calculate G */
	// !!!!!!!
	// change preudo_v   ->  -preudo_v
	pseudo_v *= -1.0;
	// !!!!!!!
	math::matrix_block_set(Wsub,1,0,blocksize,blocksize,block_sub2);
	math::dot(dLdta,pseudo_va,block_sub2);
	
	math::matrix_block_set(Wsub,1,1,blocksize,blocksize,block_sub2);
	math::dot(dBdta,pseudo_wa,block_sub2);

	math::matrix_block_set(Wsub,3,2,blocksize,blocksize,block_sub2);
	math::dot(dLdtb,pseudo_vb,block_sub2);
	
	math::matrix_block_set(Wsub,3,3,blocksize,blocksize,block_sub2);
	math::dot(dBdtb,pseudo_wb,block_sub2);

	Wsub += math::dot(G,MstatNT);
	G = dMdt - Wsub + math::transpose(Wsub);

	// std::cout 
	// 		<< "|M| = " << math::norm(M)
	// 		<< "|dMdt| = " << math::norm(dMdt)
	// 		<< "|G| = " << math::norm(G)
	// 		<< "\nM = \n" << M
	// 		<< "\ndMdt = \n" << dMdt
	// 		<< "\nG = \n" << G 
	// 		<< '\n';
	
	// std::cout << ", |M| = " << math::norm(M) << ", |dMdt| = " << math::norm(dMdt) << '\n';
	
	/*
	abc * d = ab(c*d)
	(abc * d)^T = d * (abc)^T = d * cda = (d*c)da =(d*c)(ad)^T 

	(abc*d)*e = ((abc*d)*e)^T = e*(abc*d)^T 
	*/ 


}


void ElemBEAMLD::calc_parameters(const math::vector<Node>& nodes_info) {
	// R0
	auto t10 = nodes_info[nodes[1]].xyz - nodes_info[nodes[0]].xyz;
	double L = math::norm(t10); // length
	t10 /= L;
	auto t30 = prms_orientVec();
	auto t20 = math::cross(t30,t10);
	
	
	parameters = {
		// orient vector
		t30[0],t30[1],t30[2]

		// length
		,L

		// R0
		,t10[0],t20[0],t30[0]
		,t10[1],t20[1],t30[1]
		,t10[2],t20[2],t30[2]
		// ,t10[0],t10[1],t10[2]
		// ,t20[0],t20[1],t20[2]
		// ,t30[0],t30[1],t30[2]
		,0,0
	};
}

void ElemBEAMLD::calc_parameters(const math::vector<Node>& nodes_info
								,const ElemReleases& releases) {
	
	calc_parameters(nodes_info);
	// parameters.resize(parameters.size() + 2);
	
	double shifts[2]{};
	for (auto& node_release : releases.nodes) {
		shifts[node_release.node_localid] = node_release.shift + 1;
	} 
	parameters[ElemBEAMLD::prms::release_shift1] = shifts[0];
	parameters[ElemBEAMLD::prms::release_shift2] = shifts[1];
}



size_t ElemBEAMLD::nnodes() const {
	return 2;
}
size_t ElemBEAMLD::ndofs_node() const {
	return 6;
}

double ElemBEAMLD::prms_L() const {
	return parameters[ElemBEAMLD::prms::L];
}

block_t ElemBEAMLD::prms_R0() {
	return BaseElement::prms_R0<ElemBEAMLD>(parameters.begin());
}
const_block_t ElemBEAMLD::prms_R0() const {
	return BaseElement::prms_R0<ElemBEAMLD>(parameters.begin());
}

math::Slice<math::vector<double>::const_iterator
	       ,math::vector<double>::const_iterator> ElemBEAMLD::prms_basis0_ex0() const {
	auto it = parameters.begin();
	math::Slice ex0(it + ElemBEAMLD::prms::R01x, it + ElemBEAMLD::prms::R03x+3,3);
	return ex0;
}

math::Slice<math::vector<double>::iterator
	       ,math::vector<double>::iterator> ElemBEAMLD::prms_basis0_ex0() {
	auto it = parameters.begin();
	math::Slice ex0(it + ElemBEAMLD::prms::R01x, it + ElemBEAMLD::prms::R03x+3,3);
	return ex0;
}


math::Slice<math::vector<double>::const_iterator
	       ,math::vector<double>::const_iterator> ElemBEAMLD::prms_orientVec() const {
	auto it = parameters.begin();
	math::Slice v(it + ElemBEAMLD::prms::vx, it + ElemBEAMLD::prms::vz+1);
	return v;
}

math::Slice<math::vector<double>::iterator
	       ,math::vector<double>::iterator> ElemBEAMLD::prms_orientVec() {
	auto it = parameters.begin();
	math::Slice v(it + ElemBEAMLD::prms::vx, it + ElemBEAMLD::prms::vz+1);
	return v;
}


size_t ElemBEAMLD::prms_shift1() const {
    return static_cast<size_t>(parameters[ElemBEAMLD::prms::release_shift1]);
}

size_t ElemBEAMLD::prms_shift2() const {
    return static_cast<size_t>(parameters[ElemBEAMLD::prms::release_shift2]);
}





void ElemBEAMLD2::tangentStiffness_innerLoad(math::vector_t<double,2>& K
											 , math::vector_t<double,1>& inner_load
										,const math::vector_t<double,1>& property
	 								    ,const math::vector_t<double,1>& material
										,const math::vector_t<double,1>& q ) const
{
	auto itq = q.begin();
	math::Slice ua(itq+dof::uax,itq+dof::tax);
	math::Slice ub(itq+dof::ubx,itq+dof::tbx);
	math::Slice ta(itq+dof::tax,itq+dof::ubx);
	math::Slice tb(itq+dof::tbx,itq+dof::tbz+1);
	
	/*
    static const math::vector_t<double,2> H = {
        { 0.5,  0,   0,   0,   0,   0, -0.5,  0,   0,   0,   0,   0 },
        {  0,  0.5,  0,   0,   0,   0,   0, -0.5,  0,   0,   0,   0 },
        {  0,   0,  0.5,  0,   0,   0,   0,   0, -0.5,  0,   0,   0 },
        {  0,   0,   0,  1.0,  0,   0,   0,   0,   0,   0,   0,   0 },
        {  0,   0,   0,   0,  1.0,  0,   0,   0,   0,   0,   0,   0 },
        {  0,   0,   0,   0,   0,  1.0,  0,   0,   0,   0,   0,   0 },
        {-0.5,  0,   0,   0,   0,   0,  0.5,  0,   0,   0,   0,   0 },
        {  0, -0.5,  0,   0,   0,   0,   0,  0.5,  0,   0,   0,   0 },
        {  0,   0, -0.5,  0,   0,   0,   0,   0,  0.5,  0,   0,   0 },
        {  0,   0,   0,   0,   0,   0,   0,   0,   0,  1.0,  0,   0 },
        {  0,   0,   0,   0,   0,   0,   0,   0,   0,   0,  1.0,  0 },
        {  0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,  1.0}
    }; // H == HT
	*/

	math::vector_t<double,2> La = math::rotation_tensor(ta);
	math::vector_t<double,2> Lb = math::rotation_tensor(tb);
	math::vector_t<double,2> Ba = math::zhilin_tensor(ta);
	math::vector_t<double,2> Bb = math::zhilin_tensor(tb);
	
	
	math::vector<double> theta1  = math::vector_invariant(Lb,La); // vector invariant of (R0b * RoaT) tensor

    double theta1_abs = 0.5 * std::asin(0.5*math::norm(theta1));
	
	theta1 *= -0.25 / math::detail::rotation_tensor_helper2(2.*theta1_abs);
	

    // middle element point total rotation tensor
	math::vector_t<double,2> Rm = math::dot(math::rotation_tensor(theta1),La);
	
	// block of rotation tensor from 0 to actual position
	const_block_t R0 = prms_R0();
    math::vector_t<double,2> Rall_block = math::dot(Rm,R0); 
	math::vector_t<double,2> RallT_block = math::transpose(Rall_block); 

	auto ex0 = prms_basis0_ex0();
	auto ex_ = math::dot(R0,ex0)*(parameters[prms::L]/2.);
	math::vector_t<double,2> ex_xBa = math::zeros<double>(Node::DIM,Node::DIM);
	math::vector_t<double,2> ex_xBb = math::zeros<double>(Node::DIM,Node::DIM);
	
	math::cross(ex_,Ba,ex_xBa);
	math::cross(ex_,Bb,ex_xBb);

	// Small displacement vector (phisical items same as in global dispalcement vector q)
	math::vector<double> z(q.size());
	auto itz = z.begin();
	math::Slice xia(  itz+dof::uax,itz+dof::tax);
	math::Slice xib(  itz+dof::ubx,itz+dof::tbx);
	math::Slice betaa(itz+dof::tax,itz+dof::ubx);
	math::Slice betab(itz+dof::tbx,itz+dof::tbz+1);
	
	double L = parameters[ElemBEAM::prms::L];


	xia = (ua-ub)/2. + (math::dot(Rm,ex0) - ex0)*(L/2.);
	xib -= xia; // xib was zero vector
	
	betaa -= theta1; // betaa was zero vector
	betab += theta1; // betab was zero vector
	
	// stiffness matrix
	auto K0 = ElemBEAM::stiffnessStatic(property,material,parameters); // inital (for linear element)
	// std::cout << "K10 = \n" << K0 << std::endl;
	// block_t Kblock = matrix_block(K, 0, 0, ndofs(), ndofs());
	// math::fill(Kblock.begin(), Kblock.end(), 0.0);
	auto K2 = math::zeros<double>(K0); // in 2nd state 
	
	// calculate K2 = Rall * K0 * RallT
	static const size_t blocksize = 3;
	
	block_t K0block(blocksize);
	block_t K2block(blocksize); // upper triangular blocks
	block_t Kblock_sub(blocksize); // lower triangular blocks
								// // which equal to transposed upper triangular blocks
	
	
	// loop over n blocks of K2 matrix
	size_t n_blocks = ndofs()/blocksize;
	for (size_t i = 0; i < n_blocks; ++i) {
		for (size_t j = i; j < n_blocks; ++j) {
			
			math::matrix_block_set(K0,i,j,blocksize,blocksize,K0block);
			math::matrix_block_set(K2,i,j,blocksize,blocksize,K2block);
			
			// block calculation
			K2block = math::dot(Rall_block, math::dot(K0block, RallT_block) ); // Kblock_ij = Rall_block * K0block_ij * RallT_block
			
			if (i == j) continue;
			math::matrix_block_set(K2,j,i,blocksize,blocksize,Kblock_sub);
			Kblock_sub = math::transpose(K2block);
		}
	}
	/* Calculate betaa_ = Ba * betaa 
				 betab_ = Bb * betab */
	ex_ = betaa;
	math::dot(Ba,ex_,betaa);
	ex_ = betab;
	math::dot(Bb,ex_,betab);
	
	/* Calculate inner_load */
	math::dot(K2,z,inner_load);
	inner_load *= -1.0;
	
	// calculate HTK2 = HT * K2 and K1 = HTK2*H
	/*	HTK2 = (here K == K2) 
	 [  1/2(K1_1 - K3_1), 1/2(K1_2 - K3_2), 1/2(K1_3 - K3_3), 1/2(K1_4 - K3_4)
		K2_1, 				K2_2, 			K2_3, 				K2_4
		1/2(K3_1 - K1_1), 1/2(K3_2 - K1_2), 1/2(K3_3 - K1_3), 1/2(K3_4 - K1_4)
		K4_1, 				K4_2, 			K4_3, 				K4_4			]
	*/
	
	K = K2;
	block_t HTK2block1 = std::move(K0block), HTK2block2 = std::move(K2block);
	block_t tmpblock1 = std::move(Kblock_sub);
	block_t tmpblock2(blocksize);
	for (size_t k = 0; k < n_blocks; ++k) {
		math::matrix_block_set(K2,0,k,blocksize,blocksize,tmpblock1);
		math::matrix_block_set(K2,2,k,blocksize,blocksize,tmpblock2);
		math::matrix_block_set(K,0,k,blocksize,blocksize,HTK2block1);
		HTK2block1 = 0.5*(tmpblock1-tmpblock2);
		
		math::matrix_block_set(K,2,k,blocksize,blocksize,HTK2block2);
		HTK2block2 = -1*HTK2block1;
	}

	
	// stiffness matrix K1
	
	block_t Kblock1 = std::move(HTK2block1), Kblock2 = std::move(HTK2block2);
	
	math::matrix_block_set(K,0,1,blocksize,blocksize,Kblock1);
	math::matrix_block_set(K,1,0,blocksize,blocksize,Kblock2);
	Kblock2 = math::transpose(Kblock1);
	
	math::matrix_block_set(K,2,1,blocksize,blocksize,Kblock1);
	math::matrix_block_set(K,1,2,blocksize,blocksize,Kblock2);
	Kblock2 = math::transpose(Kblock1);
	
	math::matrix_block_set(K,0,3,blocksize,blocksize,Kblock1);
	math::matrix_block_set(K,3,0,blocksize,blocksize,Kblock2);
	Kblock2 = math::transpose(Kblock1);
	
	math::matrix_block_set(K,2,3,blocksize,blocksize,Kblock1);
	math::matrix_block_set(K,3,2,blocksize,blocksize,Kblock2);
	Kblock2 = math::transpose(Kblock1);
	
	math::matrix_block_set(K,0,0,blocksize,blocksize,Kblock1);
	math::matrix_block_set(K,2,2,blocksize,blocksize,Kblock2);
	Kblock1 = (Kblock1 + Kblock2)/2.;
	Kblock2 = Kblock1;
	
	math::matrix_block_set(K,0,2,blocksize,blocksize,Kblock1);
	math::matrix_block_set(K,2,0,blocksize,blocksize,Kblock2);
	Kblock1 = (Kblock1 + Kblock2)/2.;
	Kblock2 = Kblock1;


	// std::cout << "HTK2H - math::dot(H,K2,H) = \n" << K - math::dot(math::dot(H,K2),H) << std::endl;
	
	
	
	// // 2nd part of the stiffness matrix
	auto itloadInt = inner_load.begin();
	math::Slice Fa(itloadInt+dof::uax, itloadInt+dof::tax);
	math::Slice Fb(itloadInt+dof::ubx, itloadInt+dof::tbx);
	auto S1 = math::skew_tensor(Fa*0.25);
	auto S2 = math::skew_tensor(Fb*0.25);
	
	for (size_t k = 0; k < 2; ++k) {
		math::matrix_block_set(K,0,2*k+1,blocksize,blocksize,Kblock1);
		math::matrix_block_set(K,2*k+1,0,blocksize,blocksize,Kblock2);
		Kblock1 -= S1;
		Kblock2 += S1;
		
		math::matrix_block_set(K,2,2*k+1,blocksize,blocksize,Kblock1);
		math::matrix_block_set(K,2*k+1,2,blocksize,blocksize,Kblock2);
		Kblock1 -= S2;
		Kblock2 += S2;
	}
	
}

} // namespace fem 