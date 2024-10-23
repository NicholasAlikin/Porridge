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

void ElemBEAM2D::calc_parameters(const math::vector<Node>& elem_nodes) {
	parameters = math::vector<double>{math::norm(elem_nodes[0].xyz - elem_nodes[1].xyz)};
}

size_t ElemBEAM2D::nnodes() const {
	return 2;
}
size_t ElemBEAM2D::ndofs_node() const {
	return 2;
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
	double ky = property[ElemBEAM::prop::ky];
	double kz = property[ElemBEAM::prop::kz];
	
    double L = parameters[ElemBEAM::prms::L];
    double L2 = L*L;
    // double L3 = L2*L;

    
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

void ElemBEAM::calc_parameters(const math::vector<Node>& elem_nodes) {
	parameters = math::vector<double>{math::norm(elem_nodes[0].xyz - elem_nodes[1].xyz)};
}

size_t ElemBEAM::nnodes() const {
	return 2;
}
size_t ElemBEAM::ndofs_node() const {
	return 6;
}


void ElemBEAMLD::tangentStiffness_internalLoad(math::vector_t<double,2>& K
											 , math::vector_t<double,1>& internal_load
										,const math::vector_t<double,1>& property
	 								    ,const math::vector_t<double,1>& material
                                        ,const math::vector_t<double,2>& R0
										,const math::vector_t<double,2>& basis0
                                        ,const math::vector_t<double,3>& Rsum
										,const math::vector_t<double,1>& q ) const
{
    /*
    Rsum[0] == R0a, Rsum[1] == R0b - total rotation tensors
    */
	auto itq = q.begin();
	math::Slice ua(itq+ElemBEAMLD::dof::uax,itq+ElemBEAMLD::dof::tax);
	math::Slice ub(itq+ElemBEAMLD::dof::ubx,itq+ElemBEAMLD::dof::tbx);
	math::Slice ta(itq+ElemBEAMLD::dof::tax,itq+ElemBEAMLD::dof::ubx);
	math::Slice tb(itq+ElemBEAMLD::dof::tbx,itq+ElemBEAMLD::dof::tbz+1);
	
	
    // static const math::vector_t<double,2> H = {
    //     { 0.5,  0,   0,   0,   0,   0, -0.5,  0,   0,   0,   0,   0 },
    //     {  0,  0.5,  0,   0,   0,   0,   0, -0.5,  0,   0,   0,   0 },
    //     {  0,   0,  0.5,  0,   0,   0,   0,   0, -0.5,  0,   0,   0 },
    //     {  0,   0,   0,  1.0,  0,   0,   0,   0,   0,   0,   0,   0 },
    //     {  0,   0,   0,   0,  1.0,  0,   0,   0,   0,   0,   0,   0 },
    //     {  0,   0,   0,   0,   0,  1.0,  0,   0,   0,   0,   0,   0 },
    //     {-0.5,  0,   0,   0,   0,   0,  0.5,  0,   0,   0,   0,   0 },
    //     {  0, -0.5,  0,   0,   0,   0,   0,  0.5,  0,   0,   0,   0 },
    //     {  0,   0, -0.5,  0,   0,   0,   0,   0,  0.5,  0,   0,   0 },
    //     {  0,   0,   0,   0,   0,   0,   0,   0,   0,  1.0,  0,   0 },
    //     {  0,   0,   0,   0,   0,   0,   0,   0,   0,   0,  1.0,  0 },
    //     {  0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,  1.0}
    // }; // H == HT

    math::vector<double> theta1  = math::vector_invariant(Rsum[1],Rsum[0]); // vector invariant of (R0b * RoaT) tensor

    double theta1_abs = 0.5 * std::asin(0.5*math::norm(theta1));
	if (theta1_abs > 0)
	    theta1 *= (-0.5*theta1_abs/std::sin(theta1_abs));
	else
		theta1 *= 0;

    // middle element point total rotation tensor
	math::vector_t<double,2> Rm = math::dot(math::rotation_tensor(theta1_abs),Rsum[0]); 
	
	// block of rotation tensor from 0 to actual position
    math::vector_t<double,2> Rall_block = math::dot(Rm,R0); 
	math::vector_t<double,2> RallT_block = math::transpose(Rall_block); 

	
	// Small displacement vector (phisical items same as in global dispalcement vector q)
	math::vector<double> z(q.size());
	auto itz = z.begin();
	math::Slice xia(  itz+ElemBEAMLD::dof::uax,itz+ElemBEAMLD::dof::tax);
	math::Slice xib(  itz+ElemBEAMLD::dof::ubx,itz+ElemBEAMLD::dof::tbx);
	math::Slice betaa(itz+ElemBEAMLD::dof::tax,itz+ElemBEAMLD::dof::ubx);
	math::Slice betab(itz+ElemBEAMLD::dof::tbx,itz+ElemBEAMLD::dof::tbz+1);
	
	double L = parameters[ElemBEAM::prms::L];
	xia = (ua-ub)/2. + (math::dot(Rm,basis0[0]) - basis0[0])*(L/2.);
	xib = -1.*xia;
	
	betaa = ta - theta1;
	betab = tb + theta1;
	std::cout << "z = " << z << std::endl;
	
	// stiffness matrix
	auto K0 = ElemBEAM::stiffnessStatic(property,material,parameters); // inital (for linear element)
	std::cout << "K10 = \n" << K0 << std::endl;
	auto K2 = math::zeros<double>(K0); // in 2nd state 
	
	// calculate K2 = Rall * K0 * RallT
	static const size_t blocksize = 3;
	
	block_t K0block(blocksize);
	block_t K2block(blocksize); // upper triangular blocks
	block_t Kblock_sub(blocksize); // lower triangular blocks
																						// // which equal to transposed upper triangular blocks
	
	
	// loop over n blocks of K matrix
	size_t n = ndofs()/blocksize;
	for (size_t i = 0; i < n; ++i) {
		for (size_t j = i; j < n; ++j) {
			
			math::matrix_block_set(K0,i,j,blocksize,blocksize,K0block);
			math::matrix_block_set(K2,i,j,blocksize,blocksize,K2block);
			
			// block calculation
			K2block = math::dot(Rall_block, math::dot(K0block, RallT_block) ); // Kblock_ij = Rall_block * K0block_ij * RallT_block
			
			if (i == j) continue;
			math::matrix_block_set(K2,j,i,blocksize,blocksize,Kblock_sub);
			Kblock_sub = math::transpose(K2block);
		}
	}
	
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
	for (size_t k = 0; k < n; ++k) {
		math::matrix_block_set(K2,0,k,blocksize,blocksize,tmpblock1);
		math::matrix_block_set(K2,2,k,blocksize,blocksize,tmpblock2);
		math::matrix_block_set(K,0,k,blocksize,blocksize,HTK2block1);
		HTK2block1 = 0.5*(tmpblock1-tmpblock2);
		
		math::matrix_block_set(K,2,k,blocksize,blocksize,HTK2block2);
		HTK2block2 = HTK2block1;
	}
	
	// vector of innectia forces
	internal_load = math::dot(K,z);
	
	// stiffness matrix K1
	// auto K = math::zeros<double>(K0);
	
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
	
	
	
	// 2nd part of the stiffness matrix
	math::Slice Fa(internal_load.begin()+ElemBEAMLD::dof::uax,internal_load.begin()+ElemBEAMLD::dof::tax);
	math::Slice Fb(internal_load.begin()+ElemBEAMLD::dof::ubx,internal_load.begin()+ElemBEAMLD::dof::tbx);
	auto S1 = 1./4. * math::skew_symmetric_tensor(Fa);
	auto S2 = 1./4. * math::skew_symmetric_tensor(Fb);
	
	for (size_t k = 0; k < 2; ++k) {
		math::matrix_block_set(K,0,2*k+1,blocksize,blocksize,Kblock1);
		math::matrix_block_set(K,2*k+1,0,blocksize,blocksize,Kblock2);
		Kblock1 += S1;
		Kblock2 -= S1;
		
		math::matrix_block_set(K,2,2*k+1,blocksize,blocksize,Kblock1);
		math::matrix_block_set(K,2*k+1,2,blocksize,blocksize,Kblock2);
		Kblock1 += S2;
		Kblock2 -= S2;
	}
	
}

void ElemBEAMLD::calc_parameters(const math::vector<Node>& elem_nodes) {
	parameters = math::vector<double>{math::norm(elem_nodes[0].xyz - elem_nodes[1].xyz)};
}

} // namespace fem 