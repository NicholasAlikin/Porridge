/*Finite-element model class*/

#pragma once

#include <unordered_map>

#include "fem_base.hpp"
#include "elements.hpp"

#include "dft.hpp"

#include "Eigen/Sparse"

namespace fem {
class Model;
class Assemble;

struct ModelTraits {
    static void parse_nodal_data(const Model& model
								,Assemble& assemble);
    static void store_load_vector(const Model& model
								, const Assemble& assemble
                                , math::vector<double>& load);
    static void store_load_vector(Model& model
								, const Assemble& assemble
                                , math::vector<double>& load
								, double time);
    static void store_load_vector(Model& model
								, const Assemble& assemble
                                , math::vector<double>& load
								, double freq
                                , const ::npath::DFT& dft);
    
	static void assemble_precomputing(const Model& model, Assemble& assemble);

    // Assemble global matrix and load vector methods
    template <math::VectorLike V>
	requires std::same_as<typename V::basic_value_type, double>
	&& (math::Vector<V>	|| math::Matrix<V>)
	static void assemble(const Model& model, const Assemble& assemble
                                    ,       V&       matrix
                                    ,       math::vector<double>&       load
                                    , const math::vector<double>&       q
				                    , const math::vector_t<double,3>&   Rsum
                                    , void (BaseElement::* element_matrix_load)(     /* calculate element local matrix and load vector */
                                                           math::vector_t<double,2>& /* element local matrix */
													,      math::vector_t<double,1>& /* element local load */ 
													,const math::vector_t<double,1>& /* element property */ 
													,const math::vector_t<double,1>& /* element material */ 
													,const math::vector_t<double,3>& /* all Rsum */ 
													,const math::vector_t<double,1>& /* element q */  ) const
                                    );
	// Assemble for large rotation using only Euler vector
	template <math::VectorLike V>
	requires std::same_as<typename V::basic_value_type, double>
	&& (math::Vector<V>	|| math::Matrix<V>)
	static void assemble(const Model& model, const Assemble& assemble
                                    ,       V&       matrix
                                    ,       math::vector<double>&       load
                                    , const math::vector<double>&       q
                                    , void (BaseElement::* element_matrix_load)(     /* calculate element local matrix and load vector */
                                                           math::vector_t<double,2>& /* element local matrix */
													,      math::vector_t<double,1>& /* element local load */ 
													,const math::vector_t<double,1>& /* element property */ 
													,const math::vector_t<double,1>& /* element material */
													,const math::vector_t<double,1>& /* element q */  ) const
                                    );

	template <math::VectorLike V>
	requires std::same_as<typename V::basic_value_type, double>
	&& (math::Vector<V>	|| math::Matrix<V>)
	static void assemble(const Model& model, const Assemble& assemble
                                    ,       V&       matrix
                                    ,       math::vector<double>&       load
                                    , const math::vector<double>&       q
				                    , const math::vector<double>&       dqdt
									, const math::vector_t<double,3>&   Rsum
                                    , void (BaseElement::* element_matrix_load)(     /* calculate element local matrix and load vector */
                                                           math::vector_t<double,2>& /* element local matrix */
													,      math::vector_t<double,1>& /* element local load */ 
													,const math::vector_t<double,1>& /* element property */ 
													,const math::vector_t<double,1>& /* element material */ 
													,const math::vector_t<double,3>& /* all Rsum */ 
													,const math::vector_t<double,1>& /* element q */
													,const math::vector_t<double,1>& /* element dqdt */  ) const
                                    );
	
	template <math::VectorLike V>
	requires std::same_as<typename V::basic_value_type, double>
	&& (math::Vector<V>	|| math::Matrix<V>)
	static void assemble(const Model& model, const Assemble& assemble
                                    ,       V&       matrix
                                    ,       math::vector<double>&       load
                                    , const math::vector<double>&       q
				                    , const math::vector<double>&       dqdt
									, const math::vector<double>&       d2qdt2
									, const math::vector_t<double,3>&   Rsum
                                    , void (BaseElement::* element_matrix_load)(     /* calculate element local matrix and load vector */
                                                           math::vector_t<double,2>& /* element local matrix */
													,      math::vector_t<double,1>& /* element local load */ 
													,const math::vector_t<double,1>& /* element property */ 
													,const math::vector_t<double,1>& /* element material */ 
													,const math::vector_t<double,3>& /* all Rsum */ 
													,const math::vector_t<double,1>& /* element q */
													,const math::vector_t<double,1>& /* element dqdt */
													,const math::vector_t<double,1>& /* element d2qdt2 */  ) const
                                    );

	template <math::VectorLike V>
	requires std::same_as<typename V::basic_value_type, double>
	&& (math::Vector<V>	|| math::Matrix<V>)
	static void assemble(const Model& model, const Assemble& assemble
                                    ,       V&       matrix
									,       V&       matrix2
                                    ,       math::vector<double>&       load
                                    , const math::vector<double>&       q
				                    , const math::vector<double>&       dqdt
									, const math::vector<double>&       d2qdt2
									, const math::vector_t<double,3>&   Rsum
                                    , void (BaseElement::* element_matrix_load)(     /* calculate element local matrix and load vector */
                                                           math::vector_t<double,2>& /* element local matrix */
													,	   math::vector_t<double,2>& /* element local matrix2 */
													,      math::vector_t<double,1>& /* element local load */ 
													,const math::vector_t<double,1>& /* element property */ 
													,const math::vector_t<double,1>& /* element material */ 
													,const math::vector_t<double,3>& /* all Rsum */ 
													,const math::vector_t<double,1>& /* element q */
													,const math::vector_t<double,1>& /* element dqdt */
													,const math::vector_t<double,1>& /* element d2qdt2 */  ) const
                                    );
	
	template <math::VectorLike V>
	requires std::same_as<typename V::basic_value_type, double>
	&& (math::Vector<V>	|| math::Matrix<V>)
	static void assemble(const Model& model, const Assemble& assemble
                                    ,       V&       matrix
									,       V&       matrix2
                                    ,       math::vector<double>&       load
                                    , const math::vector<double>&       q
				                    , const math::vector<double>&       dqdt
									, const math::vector<double>&       d2qdt2
                                    , void (BaseElement::* element_matrix_load)(     /* calculate element local matrix and load vector */
                                                           math::vector_t<double,2>& /* element local matrix */
													,	   math::vector_t<double,2>& /* element local matrix2 */
													,      math::vector_t<double,1>& /* element local load */ 
													,const math::vector_t<double,1>& /* element property */ 
													,const math::vector_t<double,1>& /* element material */
													,const math::vector_t<double,1>& /* element q */
													,const math::vector_t<double,1>& /* element dqdt */
													,const math::vector_t<double,1>& /* element d2qdt2 */  ) const
                                    );
	

    static void assemble(const Model& model, const Assemble& assemble
                                    ,       math::vector_t<double,2>&   matrix
                                    ,       math::vector_t<double,1>&   load
                                    , const math::vector_t<double,1>&   q
                                    , const math::vector_t<double,1>&   dqdt
                                    , const math::vector_t<double,1>&   d2qdt2
                                    , const ::npath::DFT&               dft
                                    , double                            freq
                                    ,       math::vector_t<double,3>&   buffer_dft_matrix
                                    ,       math::vector_t<double,1>&   buffer_dft_vector
                                    , void (BaseElement::* element_matrix_load)(     /* calculate element local matrix and load vector */
                                                           math::vector_t<double,2>& /* element local matrix */
													,      math::vector_t<double,1>& /* element local load */ 
													,const math::vector_t<double,1>& /* element property */ 
													,const math::vector_t<double,1>& /* element material */
													,const math::vector_t<double,1>& /* element q */
													,const math::vector_t<double,1>& /* element dqdt */
													,const math::vector_t<double,1>& /* element d2qdt2 */  
                                                    ,const ::npath::DFT&             /* dft */ 
                                                    ,double                          /* freq */ 
                                                    ,      math::vector_t<double,3>& /* buffer_dft_matrix */
                                                    ,      math::vector_t<double,1>& /* buffer_dft_vector */) const
                                    );

    /* Assemble extendent Jacobi matrix */
    static void assemble(const Model& model, const Assemble& assemble
                                ,       math::vector_t<double,2>&   matrix      /* global system Jacobi matrix */
                                ,       math::vector_t<double,1>&   load        /* global internal system load vector */
                                , const math::vector_t<double,1>&   u           /* time domain displacement */
                                , const math::vector_t<double,1>&   dudt        /* time domain velocity */
                                , const math::vector_t<double,1>&   d2udt2      /* time domain acceleration */
                                , const math::vector_const_slice<double>&   q           /* frequency domain displacement */
                                , const ::npath::DFT&               dft         /* DFT transformer */
                                , double                            freq        /* current frequency */
                                ,       math::vector_t<double,3>&   buffer_dft_matrix
                                ,       math::vector_t<double,1>&   buffer_dft_vector
                                , void (BaseElement::* element_matrix_load)(     /* calculate element local matrix and load vector */
                                                        math::vector_t<double,2>&   /* element local matrix */
                                                ,      math::vector_t<double,1>&    /* element local load */
                                                ,      math::vector_t<double,1>&    /* element extendent matrix column */ 
                                                ,const math::vector_t<double,1>&    /* element property */ 
                                                ,const math::vector_t<double,1>&    /* element material */
                                                ,const math::vector_t<double,1>&    /* element u */
                                                ,const math::vector_t<double,1>&    /* element dudt */
                                                ,const math::vector_t<double,1>&    /* element d2udt2 */
                                                ,const math::vector_const_slice<double>& /* element q */  
                                                ,const ::npath::DFT&                /* dft */ 
                                                ,double                             /* freq */ 
                                                ,      math::vector_t<double,3>&    /* buffer_dft_matrix */
                                                ,      math::vector_t<double,1>&    /* buffer_dft_vector */) const
                                );
    /* Assemble load only */
    static void assemble(const Model& model, const Assemble& assemble
                                ,       math::vector_t<double,1>&   load        /* global internal system load vector */
                                , const math::vector_t<double,1>&   u           /* time domain displacement */
                                , const math::vector_t<double,1>&   dudt        /* time domain velocity */
                                , const math::vector_t<double,1>&   d2udt2      /* time domain acceleration */
                                , const ::npath::DFT&               dft         /* DFT transformer */
                                , double                            freq        /* current frequency */
                                ,       math::vector_t<double,1>&   buffer_dft_vector
                                , void (BaseElement::* element_matrix_load)(     /* calculate element local matrix and load vector */
                                                        math::vector_t<double,1>&    /* element local load */
                                                ,const math::vector_t<double,1>&    /* element property */ 
                                                ,const math::vector_t<double,1>&    /* element material */
                                                ,const math::vector_t<double,1>&    /* element u */
                                                ,const math::vector_t<double,1>&    /* element dudt */
                                                ,const math::vector_t<double,1>&    /* element d2udt2 */
                                                ,const ::npath::DFT&                /* dft */ 
                                                ,double                             /* freq */
                                                ,      math::vector_t<double,1>&    /* buffer_dft_vector */) const
                                );
    

    /* Assemble Jacobi as Eigen::SparseMatrix */
    void assemble(const Model& model, const Assemble& assemble
                                ,       Eigen::SparseMatrix<double>&   matrix
                                ,       math::vector_t<double,1>&   load
                                , const math::vector_t<double,1>&   q
                                , const math::vector_t<double,1>&   dqdt
                                , const math::vector_t<double,1>&   d2qdt2
                                , const ::npath::DFT&               dft
                                , double                            freq
                                ,       math::vector_t<double,3>&   buffer_dft_matrix
                                ,       math::vector_t<double,1>&   buffer_dft_vector
                                , void (BaseElement::* element_matrix_load)(     /* calculate element local matrix and load vector */
                                                        math::vector_t<double,2>& /* element local matrix */
                                                ,      math::vector_t<double,1>& /* element local load */ 
                                                ,const math::vector_t<double,1>& /* element property */ 
                                                ,const math::vector_t<double,1>& /* element material */
                                                ,const math::vector_t<double,1>& /* element q */
                                                ,const math::vector_t<double,1>& /* element dqdt */
                                                ,const math::vector_t<double,1>& /* element d2qdt2 */  
                                                ,const ::npath::DFT&             /* dft */ 
                                                ,double                          /* freq */ 
                                                ,      math::vector_t<double,3>& /* buffer_dft_matrix */
                                                ,      math::vector_t<double,1>& /* buffer_dft_vector */) const
                                );

    
    /* Assemble extendent Jacobi as Eigen::SparseMatrix */
    static void assemble(const Model& model, const Assemble& assemble
                            ,       Eigen::SparseMatrix<double>&   matrix      /* global system Jacobi matrix */
                            ,       math::vector_t<double,1>&   load        /* global internal system load vector */
                            , const math::vector_t<double,1>&   u           /* time domain displacement */
                            , const math::vector_t<double,1>&   dudt        /* time domain velocity */
                            , const math::vector_t<double,1>&   d2udt2      /* time domain acceleration */
                            , const math::vector_const_slice<double>&   q           /* frequency domain displacement */
                            , const ::npath::DFT&               dft         /* DFT transformer */
                            , double                            freq        /* current frequency */
                            ,       math::vector_t<double,3>&   buffer_dft_matrix
                            ,       math::vector_t<double,1>&   buffer_dft_vector
                            , void (BaseElement::* element_matrix_load)(     /* calculate element local matrix and load vector */
                                                    math::vector_t<double,2>&   /* element local matrix */
                                            ,      math::vector_t<double,1>&    /* element local load */
                                            ,      math::vector_t<double,1>&    /* element extendent matrix column */ 
                                            ,const math::vector_t<double,1>&    /* element property */ 
                                            ,const math::vector_t<double,1>&    /* element material */
                                            ,const math::vector_t<double,1>&    /* element u */
                                            ,const math::vector_t<double,1>&    /* element dudt */
                                            ,const math::vector_t<double,1>&    /* element d2udt2 */
                                            ,const math::vector_const_slice<double>& /* element q */  
                                            ,const ::npath::DFT&                /* dft */ 
                                            ,double                             /* freq */ 
                                            ,      math::vector_t<double,3>&    /* buffer_dft_matrix */
                                            ,      math::vector_t<double,1>&    /* buffer_dft_vector */) const
                            );

private:
	template <math::Matrix V>
	static void place_element_into_symm_matrix(V& matrix, double value, size_t row, size_t col, const Assemble& assemble);

	template <math::Matrix V>
	static void place_element_into_matrix(V& matrix, double value, size_t row, size_t col, const Assemble& assemble);

	template <math::Vector V>
	static void place_element_into_symm_matrix(V& matrix, double value, size_t row, size_t col, const Assemble& assemble);


    static void place_element_into_matrix(Eigen::SparseMatrix<double>& matrix
                                , const math::vector_slice<double>& local_row
                                , size_t row, size_t col, size_t step);


	static void store_element_state_vectors(typename math::vector<double>::iterator elem_displacement
								  		  , typename math::vector<size_t>::const_iterator elemgdof_row
								  		  , typename math::vector<size_t>::const_iterator elemgdof_row_end
								  		  , const    math::vector<double>& displacement);
	static void store_element_state_vectors(typename math::vector<double>::iterator elem_displacement
										  , typename math::vector<double>::iterator elem_velocity
								  		  , typename math::vector<size_t>::const_iterator elemgdof_row
								  		  , typename math::vector<size_t>::const_iterator elemgdof_row_end
								  		  , const    math::vector<double>& displacement
										  , const    math::vector<double>& velocity);
	static void store_element_state_vectors(typename math::vector<double>::iterator elem_displacement
										  , typename math::vector<double>::iterator elem_velocity
										  , typename math::vector<double>::iterator elem_acceleration
								  		  , typename math::vector<size_t>::const_iterator elemgdof_row
								  		  , typename math::vector<size_t>::const_iterator elemgdof_row_end
								  		  , const    math::vector<double>& displacement
										  , const    math::vector<double>& velocity
										  , const    math::vector<double>& acceleration);
    static void store_element_state_vectors(typename math::vector<double>::iterator elem_displacement
										  , typename math::vector<double>::iterator elem_velocity
										  , typename math::vector<double>::iterator elem_acceleration
								  		  , typename math::vector<size_t>::const_iterator elemgdof_row
								  		  , typename math::vector<size_t>::const_iterator elemgdof_row_end
								  		  , const    math::vector<double>& displacement
										  , const    math::vector<double>& velocity
										  , const    math::vector<double>& acceleration
                                          , const    ::npath::DFT&         dft
                                          , size_t elem_ndofs);

    static void store_element_state_vectors(typename math::vector<double>::iterator elem_displacement_time
										  , typename math::vector<double>::iterator elem_velocity_time
										  , typename math::vector<double>::iterator elem_acceleration_time
                                          , typename math::vector<double>::iterator elem_acceleration_freq
								  		  , typename math::vector<size_t>::const_iterator elemgdof_row
								  		  , typename math::vector<size_t>::const_iterator elemgdof_row_end
								  		  , const    math::vector<double>& displacement
										  , const    math::vector<double>& velocity
										  , const    math::vector<double>& acceleration
                                          , const    math::vector_const_slice<double>& displ_freq
                                          , const    ::npath::DFT&         dft
                                          , size_t elem_ndofs);
};

/*Finite element model.
Containts nodes, elements,
materials and properties data*/
class Model {
public:
    friend class ModelTraits;
    friend struct AnalysisTraits;
    
	math::vector<Node> nodes_info;					// nodes coords
	math::vector<NodeLoad> loads_info;				// nodes with static load
	math::vector<NodeLoadVar> loads_var_info;		// nodes with variable load
	math::vector<NodeConstraint> constraints_info;	// constraint nodes
	math::vector<BaseElement*> elements;			// elements
	math::vector_t<double,2> materials;				// material models
	math::vector_t<double,2> properties;			// properties
	math::vector<ElemReleases> elements_releases; 	// elements with their releases

private:
	/*Calculated data*/
	math::vector<NodeAllReleases> nodes_releases; 	// nodes with all existing releases in them

public:
	Model() = default;
	Model(const Model& other);
	Model(Model&& other);
	Model(const math::vector<Node>& nodes_info					
		 ,const math::vector<NodeLoad>& loads_info				
		 ,const math::vector<NodeConstraint>& constraints_info	
		 ,const math::vector<BaseElement*>& elements			
		 ,const math::vector_t<double,2>& materials				
		 ,const math::vector_t<double,2>& properties			
		 ,const math::vector<ElemReleases>& elements_releases);
	Model(const math::vector<Node>& nodes_info					
		 ,const math::vector<NodeLoadVar>& loads_info				
		 ,const math::vector<NodeConstraint>& constraints_info	
		 ,const math::vector<BaseElement*>& elements			
		 ,const math::vector_t<double,2>& materials				
		 ,const math::vector_t<double,2>& properties			
		 ,const math::vector<ElemReleases>& elements_releases);

	/*Precomputing methods*/
    static math::vector<NodeAllReleases> releases_precomputing(math::vector<ElemReleases>& elem_releases
											           , const math::vector<BaseElement*>& elements);

public:
	/* Constant flags */
	static constexpr size_t DOF_IS_CONSTRAINED = 0;
};


/*Assemble of global model matrix (stiffness)*/
class Assemble {
public:
    friend class ModelTraits;
    friend struct AnalysisTraits;

	size_t ndofs;									// number of model dofs
	math::vector_t<size_t,2> nodes_dofs;			// global dofs stored by nodes
	math::vector_t<size_t,2> elems_dofs;			// global dofs stored by elements
	
	math::vector<size_t> colhs;						// heights of global model matrix columns
	math::vector<size_t> diags;						// positions of global model matrix diagonal elements
	size_t band_width;								// width of global model matrix stored band-like

private:
	size_t max_elem_dofs;							// max number of dofs among all elements
    /* buffers */
    /* element locals */
	mutable math::vector_t<double,2> elem_matrix;			/* matrix */
	mutable math::vector_t<double,2> elem_matrix2;			/* matrix second
                                                            Uses to assemble 2 matrices simulteneously */
	mutable math::vector_t<double,1> elem_load;				/* load vector */
    mutable math::vector_t<double,1> elem_load2;				/* load vector 2.
                                                            Used to store extendent matrix column */
	mutable math::vector_t<double,1> elem_displ;			/* displacement vector */
    mutable math::vector_t<double,1> elem_displ2;			/* also displacement vector.
                                                            Uses to store displacement vector
                                                            in time and frequency domains */
	mutable math::vector_t<double,1> elem_vel;			    /* velocity vector */
	mutable math::vector_t<double,1> elem_accel;			/* acceleration vector */
public:
	
	Assemble() = default;
	Assemble(const Assemble& other);
	Assemble(Assemble&& other);

    /* Attribute getters */
    size_t get_max_elem_dofs();
	/*Precomputing methods*/
    void elem_matrix_size(size_t sz);
    void elem_load_size(size_t sz);
    void elem_state_size(size_t sz);
    void frequency_analyses_elem_state_size(size_t freq_size, size_t time_size);
};

/*Assemble of nonlinear system:
	o assamble global matrix (stiffness,mass)
	o assemble global load vector (internal, inectia)
Global matrix might be store like 1d (column band-like form) or 2d vector (ordinaly matrix-like form).
*/
template <math::VectorLike V>
requires std::same_as<typename V::basic_value_type, double>
&& (math::Vector<V>	|| math::Matrix<V>)
void ModelTraits::assemble(const Model& model, const Assemble& assemble		 /* Model and Assemble objects */
                                    ,       V&       matrix							 /* Global matrix */
                                    ,       math::vector<double>&       load		 /* global load vector */
                                    , const math::vector<double>&       q			 /* Global displacement vector */
				                    , const math::vector_t<double,3>&   Rsum		 /* Global nodes total rotational tensors */
                                    , void (BaseElement::* element_matrix_load)(     /* Function to calculate element matrix and load vector */
                                                           math::vector_t<double,2>& /*  element local matrix */
													,      math::vector_t<double,1>& /*  element local load */ 
													,const math::vector_t<double,1>& /*  element property */ 
													,const math::vector_t<double,1>& /*  element material */ 
													,const math::vector_t<double,3>& /*  all Rsum */ 
													,const math::vector_t<double,1>& /*  element q */  ) const
                                    )
{
	/*Iterators*/	
	// decltype(assemble.elem_displ.begin()) 		q_elem_dof,q_elem_dof_end;				/* Element displacement vector */
	decltype(assemble.elem_load.begin())  		load_dof;								/* Element load vector */
 
	decltype(assemble.elem_matrix.begin()) 		matloc_row;				/* Element matrix row values */
	decltype(assemble.elem_matrix[0].begin()) 	matloc_col;				/* Element matrix col values */
 
	auto 										elem = model.elements.begin()			/* Element*/
											  , elem_end = model.elements.end();	
	// decltype(model.elements[0]->nodes.begin()) 	elemnode, elemnode_end; 				/* Element node */
	
	auto 										elemgdofs = assemble.elems_dofs.begin();/* Element dofs */
	decltype(assemble.elems_dofs[0].begin()) 	elemgdof_row, elemgdof_row_end 			/* Element dof for loop over element matrix rows */
											  , elemgdof_col, elemgdof_col_end;  		/*             for loop over element matrix cols */
	// decltype(assemble.nodes_dofs[0].begin()) 	nodegdof, nodegdof_end; 				/* Element global dof id*/
	
	
	/* Element matrix rows counter */
	size_t matloc_row_num;
	
	/* Loop over all elements */
	while (elem != elem_end) {

		/* store element displacements */
		ModelTraits::store_element_state_vectors(assemble.elem_displ.begin()	// where store to
												,elemgdofs->begin()				// with dofs id
												,elemgdofs->end()
												,q);							// store from there

		/* calculate element matrix and load vector */
		((*elem)->*element_matrix_load)( assemble.elem_matrix 					// where to store matrix
										,assemble.elem_load						// where to store load vector
										,model.properties[(*elem)->propID]		// element property
										,model.materials[ (*elem)->matlID]		// element material
										,Rsum									// total rotational tensors
										,assemble.elem_displ);					// element displacement vector
		
		/* store element matrix and load vector to global matrix and load vector*/
		// loop over matrix rows
		for (matloc_row = assemble.elem_matrix.begin()    			// iterator on elem matrix row
		    ,elemgdof_row = elemgdofs->begin()   					// iterator on global dofs, corresponding row global dofs of the elem
			,elemgdof_row_end = elemgdofs->end() 					// same
			,matloc_row_num = 0					 					// elem matrix count - to go throw the upper triangular part of the matrix only
			,load_dof = assemble.elem_load.begin()			  		// iterator on elem load vector
		   				;elemgdof_row < elemgdof_row_end 			// loop over all matrix rows
									;++matloc_row
									,++elemgdof_row
									,++matloc_row_num
									,++load_dof)
		{

			// if dof is constrained - go to next matrix row and vector component
			if (*elemgdof_row == 0) continue;     					// assamble only matrix rows, which correspond not constrained global dofs

			// loop over column of the current matrix row
			for (matloc_col = matloc_row->begin() + matloc_row_num  // iterator on column element - only upper triangular part
			    ,elemgdof_col = elemgdofs->begin() + matloc_row_num // iterator on global dofs, corresponding col global dofs of the elem
				,elemgdof_col_end = elemgdofs->end() 				// same
							;elemgdof_col < elemgdof_col_end 		// loop over row elements from diagonal to the end
										;++matloc_col
										,++elemgdof_col)
			{
				// if dof is constrained - go to next component
				if (*elemgdof_col == Model::DOF_IS_CONSTRAINED) continue;  					// assamble only matrix columns, which correspond not constrained global dofs
				//
				ModelTraits::place_element_into_symm_matrix(matrix,*matloc_col,*elemgdof_row,*elemgdof_col,assemble);
				
			}

			// store global load vector
			load[*elemgdof_row-1] += *load_dof;
		}
		
		
		++elem; ++elemgdofs;
	}
	
}

/*Assemble of nonlinear system:
	o assamble global matrix (stiffness,mass)
	o assemble global load vector (internal, inectia)
Global matrix might be store like 1d (column band-like form) or 2d vector (ordinaly matrix-like form).
*/
template <math::VectorLike V>
requires std::same_as<typename V::basic_value_type, double>
&& (math::Vector<V>	|| math::Matrix<V>)
void ModelTraits::assemble(const Model& model, const Assemble& assemble		 /* Model and Assemble objects */
                                    ,       V&       matrix							 /* Global matrix */
                                    ,       math::vector<double>&       load		 /* global load vector */
                                    , const math::vector<double>&       q			 /* Global displacement vector */
                                    , void (BaseElement::* element_matrix_load)(     /* Function to calculate element matrix and load vector */
                                                           math::vector_t<double,2>& /*  element local matrix */
													,      math::vector_t<double,1>& /*  element local load */ 
													,const math::vector_t<double,1>& /*  element property */ 
													,const math::vector_t<double,1>& /*  element material */
													,const math::vector_t<double,1>& /*  element q */  ) const
                                    )
{
	/*Iterators*/	
	// decltype(assemble.elem_displ.begin()) 		q_elem_dof,q_elem_dof_end;				/* Element displacement vector */
	decltype(assemble.elem_load.begin())  		load_dof;								/* Element load vector */
 
	decltype(assemble.elem_matrix.begin()) 		matloc_row;				/* Element matrix row values */
	decltype(assemble.elem_matrix[0].begin()) 	matloc_col;				/* Element matrix col values */
 
	auto 										elem = model.elements.begin()			/* Element*/
											  , elem_end = model.elements.end();	
	// decltype(model.elements[0]->nodes.begin()) 	elemnode, elemnode_end; 				/* Element node */
	
	auto 										elemgdofs = assemble.elems_dofs.begin();/* Element dofs */
	decltype(assemble.elems_dofs[0].begin()) 	elemgdof_row, elemgdof_row_end 			/* Element dof for loop over element matrix rows */
											  , elemgdof_col, elemgdof_col_end;  		/*             for loop over element matrix cols */
	// decltype(assemble.nodes_dofs[0].begin()) 	nodegdof, nodegdof_end; 				/* Element global dof id*/
	
	
	/* Element matrix rows counter */
	size_t matloc_row_num;
	
	/* Loop over all elements */
	while (elem != elem_end) {

		/* store element displacements */
		ModelTraits::store_element_state_vectors(assemble.elem_displ.begin()	// where store to
												,elemgdofs->begin()				// with dofs id
												,elemgdofs->end()
												,q);							// store from there

		/* calculate element matrix and load vector */
		((*elem)->*element_matrix_load)( assemble.elem_matrix 					// where to store matrix
										,assemble.elem_load						// where to store load vector
										,model.properties[(*elem)->propID]		// element property
										,model.materials[ (*elem)->matlID]		// element material
										,assemble.elem_displ);					// element displacement vector
		
		/* store element matrix and load vector to global matrix and load vector*/
		// loop over matrix rows
		for (matloc_row = assemble.elem_matrix.begin()    			// iterator on elem matrix row
		    ,elemgdof_row = elemgdofs->begin()   					// iterator on global dofs, corresponding row global dofs of the elem
			,elemgdof_row_end = elemgdofs->end() 					// same
			,matloc_row_num = 0					 					// elem matrix count - to go throw the upper triangular part of the matrix only
			,load_dof = assemble.elem_load.begin()			  		// iterator on elem load vector
		   				;elemgdof_row < elemgdof_row_end 			// loop over all matrix rows
									;++matloc_row
									,++elemgdof_row
									,++matloc_row_num
									,++load_dof)
		{

			// if dof is constrained - go to next matrix row and vector component
			if (*elemgdof_row == 0) continue;     					// assamble only matrix rows, which correspond not constrained global dofs

			// loop over column of the current matrix row
			for (matloc_col = matloc_row->begin() + matloc_row_num  // iterator on column element - only upper triangular part
			    ,elemgdof_col = elemgdofs->begin() + matloc_row_num // iterator on global dofs, corresponding col global dofs of the elem
				,elemgdof_col_end = elemgdofs->end() 				// same
							;elemgdof_col < elemgdof_col_end 		// loop over row elements from diagonal to the end
										;++matloc_col
										,++elemgdof_col)
			{
				// if dof is constrained - go to next component
				if (*elemgdof_col == Model::DOF_IS_CONSTRAINED) continue;  					// assamble only matrix columns, which correspond not constrained global dofs
				//
				ModelTraits::place_element_into_symm_matrix(matrix,*matloc_col,*elemgdof_row,*elemgdof_col,assemble);
				
			}

			// store global load vector
			load[*elemgdof_row-1] += *load_dof;
		}
		
		
		++elem; ++elemgdofs;
	}
	// throw 1;
}

template <math::VectorLike V>
requires std::same_as<typename V::basic_value_type, double>
&& (math::Vector<V>	|| math::Matrix<V>)
void ModelTraits::assemble(const Model& model, const Assemble& assemble
                                ,       V&       matrix
                                ,       math::vector<double>&       load
                                , const math::vector<double>&       q
			                    , const math::vector<double>&       dqdt
								, const math::vector_t<double,3>&   Rsum
                                , void (BaseElement::* element_matrix_load)(     /* calculate element local matrix and load vector */
                                                       math::vector_t<double,2>& /* element local matrix */
												,      math::vector_t<double,1>& /* element local load */ 
												,const math::vector_t<double,1>& /* element property */ 
												,const math::vector_t<double,1>& /* element material */ 
												,const math::vector_t<double,3>& /* all Rsum */ 
												,const math::vector_t<double,1>& /* element q */
												,const math::vector_t<double,1>& /* element dqdt */  ) const
                                )
{
	/*Iterators*/	
	// decltype(assemble.elem_displ.begin()) 		q_elem_dof,q_elem_dof_end;				/* Element displacement vector */
	decltype(assemble.elem_load.begin())  		load_dof;								/* Element load vector */
 
	decltype(assemble.elem_matrix.begin()) 		matloc_row;				/* Element matrix row values */
	decltype(assemble.elem_matrix[0].begin()) 	matloc_col;				/* Element matrix col values */
 
	auto 										elem = model.elements.begin()			/* Element*/
											  , elem_end = model.elements.end();	
	// decltype(model.elements[0]->nodes.begin()) 	elemnode, elemnode_end; 				/* Element node */
	
	auto 										elemgdofs = assemble.elems_dofs.begin();/* Element dofs */
	decltype(assemble.elems_dofs[0].begin()) 	elemgdof_row, elemgdof_row_end 			/* Element dof for loop over element matrix rows */
											  , elemgdof_col, elemgdof_col_end;  		/*             for loop over element matrix cols */
	// decltype(assemble.nodes_dofs[0].begin()) 	nodegdof, nodegdof_end; 				/* Element global dof id*/
	
	
	/* Element matrix rows counter */
	size_t matloc_row_num;
	
	/* Loop over all elements */
	while (elem != elem_end) {

		/* store element displacements */
		ModelTraits::store_element_state_vectors(assemble.elem_displ.begin()	// where store to
												,assemble.elem_vel.begin()
												,elemgdofs->begin()				// with dofs id
												,elemgdofs->end()
												,q
												,dqdt);							// store from there

		/* calculate element matrix and load vector */
		((*elem)->*element_matrix_load)( assemble.elem_matrix 					// where to store matrix
										,assemble.elem_load						// where to store load vector
										,model.properties[(*elem)->propID]		// element property
										,model.materials[ (*elem)->matlID]		// element material
										,Rsum									// total rotational tensors
										,assemble.elem_displ					// element displacement vector
										,assemble.elem_vel);
		
		/* store element matrix and load vector to global matrix and load vector*/
		// loop over matrix rows
		for (matloc_row = assemble.elem_matrix.begin()    			// iterator on elem matrix row
		    ,elemgdof_row = elemgdofs->begin()   					// iterator on global dofs, corresponding row global dofs of the elem
			,elemgdof_row_end = elemgdofs->end() 					// same
			,matloc_row_num = 0					 					// elem matrix count - to go throw the upper triangular part of the matrix only
			,load_dof = assemble.elem_load.begin()			  		// iterator on elem load vector
		   				;elemgdof_row < elemgdof_row_end 			// loop over all matrix rows
									;++matloc_row
									,++elemgdof_row
									,++matloc_row_num
									,++load_dof)
		{

			// if dof is constrained - go to next matrix row and vector component
			if (*elemgdof_row == 0) continue;     					// assamble only matrix rows, which correspond not constrained global dofs

			// loop over column of the current matrix row
			for (matloc_col = matloc_row->begin() + matloc_row_num  // iterator on column element - only upper triangular part
			    ,elemgdof_col = elemgdofs->begin() + matloc_row_num // iterator on global dofs, corresponding col global dofs of the elem
				,elemgdof_col_end = elemgdofs->end() 				// same
							;elemgdof_col < elemgdof_col_end 		// loop over row elements from diagonal to the end
										;++matloc_col
										,++elemgdof_col)
			{
				// if dof is constrained - go to next component
				if (*elemgdof_col == Model::DOF_IS_CONSTRAINED) continue;  					// assamble only matrix columns, which correspond not constrained global dofs
				//
				ModelTraits::place_element_into_symm_matrix(matrix,*matloc_col,*elemgdof_row,*elemgdof_col,assemble);
				
			}

			// store global load vector
			load[*elemgdof_row-1] += *load_dof;
		}
		
		
		++elem; ++elemgdofs;
	}
	
}


template <math::VectorLike V>
requires std::same_as<typename V::basic_value_type, double>
&& (math::Vector<V>	|| math::Matrix<V>)
void ModelTraits::assemble(const Model& model, const Assemble& assemble
                                ,       V&       matrix
                                ,       math::vector<double>&       load
                                , const math::vector<double>&       q
			                    , const math::vector<double>&       dqdt
								, const math::vector<double>&       d2qdt2
								, const math::vector_t<double,3>&   Rsum
                                , void (BaseElement::* element_matrix_load)(     /* calculate element local matrix and load vector */
                                                       math::vector_t<double,2>& /* element local matrix */
												,      math::vector_t<double,1>& /* element local load */ 
												,const math::vector_t<double,1>& /* element property */ 
												,const math::vector_t<double,1>& /* element material */ 
												,const math::vector_t<double,3>& /* all Rsum */ 
												,const math::vector_t<double,1>& /* element q */
												,const math::vector_t<double,1>& /* element dqdt */
												,const math::vector_t<double,1>& /* element d2qdt2 */  ) const
                                )
{
	/*Iterators*/	
	// decltype(assemble.elem_displ.begin()) 		q_elem_dof,q_elem_dof_end;				/* Element displacement vector */
	decltype(assemble.elem_load.begin())  		load_dof;								/* Element load vector */
 
	decltype(assemble.elem_matrix.begin()) 		matloc_row;				/* Element matrix row values */
	decltype(assemble.elem_matrix[0].begin()) 	matloc_col;				/* Element matrix col values */
 
	auto 										elem = model.elements.begin()			/* Element*/
											  , elem_end = model.elements.end();	
	// decltype(model.elements[0]->nodes.begin()) 	elemnode, elemnode_end; 				/* Element node */
	
	auto 										elemgdofs = assemble.elems_dofs.begin();/* Element dofs */
	decltype(assemble.elems_dofs[0].begin()) 	elemgdof_row, elemgdof_row_end 			/* Element dof for loop over element matrix rows */
											  , elemgdof_col, elemgdof_col_end;  		/*             for loop over element matrix cols */
	// decltype(assemble.nodes_dofs[0].begin()) 	nodegdof, nodegdof_end; 				/* Element global dof id*/
	
	
	/* Element matrix rows counter */
	size_t matloc_row_num;
	/* Loop over all elements */
	while (elem != elem_end) {

		/* store element displacements */
		ModelTraits::store_element_state_vectors(assemble.elem_displ.begin()	// where store to
												,assemble.elem_vel.begin()
												,assemble.elem_accel.begin()
												,elemgdofs->begin()				// with dofs id
												,elemgdofs->end()
												,q								// store from there
												,dqdt
												,d2qdt2);
		/* calculate element matrix and load vector */
		((*elem)->*element_matrix_load)( assemble.elem_matrix 					// where to store matrix
										,assemble.elem_load						// where to store load vector
										,model.properties[(*elem)->propID]		// element property
										,model.materials[ (*elem)->matlID]		// element material
										,Rsum									// total rotational tensors
										,assemble.elem_displ					// element displacement vector
										,assemble.elem_vel
										,assemble.elem_accel);
		/* store element matrix and load vector to global matrix and load vector*/
		// loop over matrix rows
		for (matloc_row = assemble.elem_matrix.begin()    			// iterator on elem matrix row
		    ,elemgdof_row = elemgdofs->begin()   					// iterator on global dofs, corresponding row global dofs of the elem
			,elemgdof_row_end = elemgdofs->end() 					// same
			,matloc_row_num = 0					 					// elem matrix count - to go throw the upper triangular part of the matrix only
			,load_dof = assemble.elem_load.begin()			  		// iterator on elem load vector
		   				;elemgdof_row < elemgdof_row_end 			// loop over all matrix rows
									;++matloc_row
									,++elemgdof_row
									,++matloc_row_num
									,++load_dof)
		{

			// if dof is constrained - go to next matrix row and vector component
			if (*elemgdof_row == 0) continue;     					// assamble only matrix rows, which correspond not constrained global dofs

			// loop over column of the current matrix row
			for (matloc_col = matloc_row->begin() + matloc_row_num  // iterator on column element - only upper triangular part
			    ,elemgdof_col = elemgdofs->begin() + matloc_row_num // iterator on global dofs, corresponding col global dofs of the elem
				,elemgdof_col_end = elemgdofs->end() 				// same
							;elemgdof_col < elemgdof_col_end 		// loop over row elements from diagonal to the end
										;++matloc_col
										,++elemgdof_col)
			{
				// if dof is constrained - go to next component
				if (*elemgdof_col == Model::DOF_IS_CONSTRAINED) continue;  					// assamble only matrix columns, which correspond not constrained global dofs
				//
				ModelTraits::place_element_into_symm_matrix(matrix,*matloc_col,*elemgdof_row,*elemgdof_col,assemble);
				
			}

			// store global load vector
			load[*elemgdof_row-1] += *load_dof;
		}
		
		
		++elem; ++elemgdofs;
	}
	
}

template <math::VectorLike V>
requires std::same_as<typename V::basic_value_type, double>
&& (math::Vector<V>	|| math::Matrix<V>)
void ModelTraits::assemble(const Model& model, const Assemble& assemble
                                ,       V&       matrix
								,       V&       matrix2
                                ,       math::vector<double>&       load
                                , const math::vector<double>&       q
			                    , const math::vector<double>&       dqdt
								, const math::vector<double>&       d2qdt2
								, const math::vector_t<double,3>&   Rsum
                                , void (BaseElement::* element_matrix_load)(     /* calculate element local matrix and load vector */
                                                       math::vector_t<double,2>& /* element local matrix */
												,	   math::vector_t<double,2>& /* element local matrix2 */
												,      math::vector_t<double,1>& /* element local load */ 
												,const math::vector_t<double,1>& /* element property */ 
												,const math::vector_t<double,1>& /* element material */ 
												,const math::vector_t<double,3>& /* all Rsum */ 
												,const math::vector_t<double,1>& /* element q */
												,const math::vector_t<double,1>& /* element dqdt */
												,const math::vector_t<double,1>& /* element d2qdt2 */  ) const
                                )
{
	/*Iterators*/	
	// decltype(assemble.elem_displ.begin()) 		q_elem_dof,q_elem_dof_end;				/* Element displacement vector */
	decltype(assemble.elem_load.begin())  		load_dof;								/* Element load vector */
 
	decltype(assemble.elem_matrix.begin()) 		matloc_row,matloc2_row;				/* Element matrix row values */
	decltype(assemble.elem_matrix[0].begin()) 	matloc_col,matloc2_col;				/* Element matrix col values */
 
	auto 										elem = model.elements.begin()			/* Element*/
											  , elem_end = model.elements.end();	
	// decltype(model.elements[0]->nodes.begin()) 	elemnode, elemnode_end; 				/* Element node */
	
	auto 										elemgdofs = assemble.elems_dofs.begin();/* Element dofs */
	decltype(assemble.elems_dofs[0].begin()) 	elemgdof_row, elemgdof_row_end 			/* Element dof for loop over element matrix rows */
											  , elemgdof_col, elemgdof_col_end;  		/*             for loop over element matrix cols */
	// decltype(assemble.nodes_dofs[0].begin()) 	nodegdof, nodegdof_end; 				/* Element global dof id*/
	
	
	/* Loop over all elements */
	while (elem != elem_end) {

		/* store element displacements */
		ModelTraits::store_element_state_vectors(assemble.elem_displ.begin()	// where store to
												,assemble.elem_vel.begin()
												,assemble.elem_accel.begin()
												,elemgdofs->begin()				// with dofs id
												,elemgdofs->end()
												,q								// store from there
												,dqdt
												,d2qdt2);
		/* calculate element matrix and load vector */
		((*elem)->*element_matrix_load)( assemble.elem_matrix 					// where to store matrix
										,assemble.elem_matrix2 					// where to store matrix2
										,assemble.elem_load						// where to store load vector
										,model.properties[(*elem)->propID]		// element property
										,model.materials[ (*elem)->matlID]		// element material
										,Rsum									// total rotational tensors
										,assemble.elem_displ					// element displacement vector
										,assemble.elem_vel
										,assemble.elem_accel);
		/* store element matrix and load vector to global matrix and load vector*/
		// loop over matrix rows
		for (matloc_row = assemble.elem_matrix.begin()    			// iterator on elem matrix row
			,matloc2_row = assemble.elem_matrix2.begin()
		    ,elemgdof_row = elemgdofs->begin()   					// iterator on global dofs, corresponding row global dofs of the elem
			,elemgdof_row_end = elemgdofs->end() 					// same
			,load_dof = assemble.elem_load.begin()			  		// iterator on elem load vector
		   				;elemgdof_row < elemgdof_row_end 			// loop over all matrix rows
									;++matloc_row
									,++matloc2_row
									,++elemgdof_row
									,++load_dof)
		{

			// if dof is constrained - go to next matrix row and vector component
			if (*elemgdof_row == 0) continue;     					// assamble only matrix rows, which correspond not constrained global dofs

			// loop over column of the current matrix row
			// loop over all element - assume that matrix and martix2 are not symmetric
			for (matloc_col = matloc_row->begin()  					// iterator on column element
				,matloc2_col = matloc2_row->begin()
			    ,elemgdof_col = elemgdofs->begin()					// iterator on global dofs, corresponding col global dofs of the elem
				,elemgdof_col_end = elemgdofs->end() 				// same
							;elemgdof_col < elemgdof_col_end 		// loop over row elements from diagonal to the end
										;++matloc_col
										,++matloc2_col
										,++elemgdof_col)
			{
				// if dof is constrained - go to next component
				if (*elemgdof_col == Model::DOF_IS_CONSTRAINED) continue;  					// assamble only matrix columns, which correspond not constrained global dofs
				
				ModelTraits::place_element_into_matrix(matrix ,*matloc_col ,*elemgdof_row,*elemgdof_col,assemble);
				ModelTraits::place_element_into_matrix(matrix2,*matloc2_col,*elemgdof_row,*elemgdof_col,assemble);
				
			}

			// store global load vector
			load[*elemgdof_row-1] += *load_dof;
		}
		
		
		++elem; ++elemgdofs;
	}
	
}

template <math::VectorLike V>
requires std::same_as<typename V::basic_value_type, double>
&& (math::Vector<V>	|| math::Matrix<V>)
void ModelTraits::assemble(const Model& model, const Assemble& assemble
                                ,       V&       matrix
								,       V&       matrix2
                                ,       math::vector<double>&       load
                                , const math::vector<double>&       q
			                    , const math::vector<double>&       dqdt
								, const math::vector<double>&       d2qdt2
                                , void (BaseElement::* element_matrix_load)(     /* calculate element local matrix and load vector */
                                                       math::vector_t<double,2>& /* element local matrix */
												,	   math::vector_t<double,2>& /* element local matrix2 */
												,      math::vector_t<double,1>& /* element local load */ 
												,const math::vector_t<double,1>& /* element property */ 
												,const math::vector_t<double,1>& /* element material */ 
												,const math::vector_t<double,1>& /* element q */
												,const math::vector_t<double,1>& /* element dqdt */
												,const math::vector_t<double,1>& /* element d2qdt2 */  ) const
                                )
{
	/*Iterators*/	
	// decltype(assemble.elem_displ.begin()) 		q_elem_dof,q_elem_dof_end;				/* Element displacement vector */
	decltype(assemble.elem_load.begin())  		load_dof;								/* Element load vector */
 
	decltype(assemble.elem_matrix.begin()) 		matloc_row,matloc2_row;				/* Element matrix row values */
	decltype(assemble.elem_matrix[0].begin()) 	matloc_col,matloc2_col;				/* Element matrix col values */
 
	auto 										elem = model.elements.begin()			/* Element*/
											  , elem_end = model.elements.end();	
	// decltype(model.elements[0]->nodes.begin()) 	elemnode, elemnode_end; 				/* Element node */
	
	auto 										elemgdofs = assemble.elems_dofs.begin();/* Element dofs */
	decltype(assemble.elems_dofs[0].begin()) 	elemgdof_row, elemgdof_row_end 			/* Element dof for loop over element matrix rows */
											  , elemgdof_col, elemgdof_col_end;  		/*             for loop over element matrix cols */
	// decltype(assemble.nodes_dofs[0].begin()) 	nodegdof, nodegdof_end; 				/* Element global dof id*/
	
	
	/* Loop over all elements */
	while (elem != elem_end) {

		/* store element displacements */
		ModelTraits::store_element_state_vectors(assemble.elem_displ.begin()	// where store to
												,assemble.elem_vel.begin()
												,assemble.elem_accel.begin()
												,elemgdofs->begin()				// with dofs id
												,elemgdofs->end()
												,q								// store from there
												,dqdt
												,d2qdt2);
		/* calculate element matrix and load vector */
		((*elem)->*element_matrix_load)( assemble.elem_matrix 					// where to store matrix
										,assemble.elem_matrix2 					// where to store matrix2
										,assemble.elem_load						// where to store load vector
										,model.properties[(*elem)->propID]		// element property
										,model.materials[ (*elem)->matlID]		// element material
										,assemble.elem_displ					// element displacement vector
										,assemble.elem_vel
										,assemble.elem_accel);
		/* store element matrix and load vector to global matrix and load vector*/
		// loop over matrix rows
		for (matloc_row = assemble.elem_matrix.begin()    			// iterator on elem matrix row
			,matloc2_row = assemble.elem_matrix2.begin()
		    ,elemgdof_row = elemgdofs->begin()   					// iterator on global dofs, corresponding row global dofs of the elem
			,elemgdof_row_end = elemgdofs->end() 					// same
			,load_dof = assemble.elem_load.begin()			  		// iterator on elem load vector
		   				;elemgdof_row < elemgdof_row_end 			// loop over all matrix rows
									;++matloc_row
									,++matloc2_row
									,++elemgdof_row
									,++load_dof)
		{

			// if dof is constrained - go to next matrix row and vector component
			if (*elemgdof_row == 0) continue;     					// assamble only matrix rows, which correspond not constrained global dofs

			// loop over column of the current matrix row
			// loop over all element - assume that matrix and martix2 are not symmetric
			for (matloc_col = matloc_row->begin()  					// iterator on column element
				,matloc2_col = matloc2_row->begin()
			    ,elemgdof_col = elemgdofs->begin()					// iterator on global dofs, corresponding col global dofs of the elem
				,elemgdof_col_end = elemgdofs->end() 				// same
							;elemgdof_col < elemgdof_col_end 		// loop over row elements from diagonal to the end
										;++matloc_col
										,++matloc2_col
										,++elemgdof_col)
			{
				// if dof is constrained - go to next component
				if (*elemgdof_col == Model::DOF_IS_CONSTRAINED) continue;  					// assamble only matrix columns, which correspond not constrained global dofs
				
				ModelTraits::place_element_into_matrix(matrix ,*matloc_col ,*elemgdof_row,*elemgdof_col,assemble);
				ModelTraits::place_element_into_matrix(matrix2,*matloc2_col,*elemgdof_row,*elemgdof_col,assemble);
				
			}

			// store global load vector
			load[*elemgdof_row-1] += *load_dof;
		}
		
		
		++elem; ++elemgdofs;
	}
	
}



template <math::Matrix V>
void ModelTraits::place_element_into_symm_matrix(V& matrix, double value, size_t row, size_t col, const Assemble& assemble) {
	matrix[row-1][col-1] += value;
	if (row != col) {
		matrix[col-1][row-1] += value;
	}
}

template <math::Matrix V>
void ModelTraits::place_element_into_matrix(V& matrix, double value, size_t row, size_t col, const Assemble& assemble) {
	matrix[row-1][col-1] += value;
}

template <math::Vector V>
void ModelTraits::place_element_into_symm_matrix(V& matrix, double value, size_t row, size_t col, const Assemble& assemble) {
	matrix[(assemble.diags[col-1] + col) - row] += value;
}

} // namespace fem