// Incompressible Flow Solver

// Copyright (C) 2005  Amir R. Baserinia

//----------------------------------------------------------------------------
// Date:   08 Apr 2004
// Update: 25 Feb 2005
//----------------------------------------------------------------------------

#ifndef CLASS_NIFS_C_CELL_H
#define CLASS_NIFS_C_CELL_H

#include <vector>
#include "mesh_obj.h"
#include "point_2d.h"
#include "vector_block.h"
#include "matrix_block.h"
#include "vector_2d.h"
#include "matrix_2d.h"

namespace NIFS{    

//--- INTERFACE ---
	
class c_cell: public c_mesh_obj {
public:
	c_cell();
	virtual ~c_cell() {}
	double get_volume() { return m_volume; }    
	const c_point_2d& get_node() { return m_node_point; }
	
	void calc_volume_and_cent(); // area and centroid
	void calc_gradient(); // general gradient function
	void calc_hessian();
	void calc_coefs( void );
	void calc_anis_spacing( void );
	
	c_matrix_2d get_anis_spacing_memu( void );
		
	void update_geometry( void );
	
	vector_block< double > get_state( void )  { return m_state; }
	vector_block< c_vector_2d > get_grad( void ) { return m_grad; } 
	vector_block< c_matrix_2d > get_hess( void ) { return m_hess; } 
	
	void set_state( const vector_block< double >& state );
	void set_state( const vector_block< double >& state, unsigned index_beg, unsigned index_end);
	void set_force_state( const vector_block< double >& state );
	void set_force_state( const vector_block< double >& state, unsigned index_beg, unsigned index_end);

	const vector_block< double >& get_coef_rhs( void ) const { return m_coef_rhs; }
	const matrix_block< double >& get_coef_slf( void ) const { return m_coef_slf; }
	const matrix_block< double >& get_coef_nbr( int i ) const { return m_coef_nbr[i]; }
	
	vector_block< double > get_residual();
	
	double get_df( void ) { return 0.5*( m_coef_slf(1)(1) + m_coef_slf(2)(2) ); }
	void print_result( void );
	
	bool is_valid_move( c_mesh_obj* vert , c_vector_2d vec );
	double get_charac_size( void ); // characteristic size
	vector_block< double > get_err_iso_coef( void );
  double get_aspect_ratio();	
	
private:
	double m_volume; // volume of cell ( area here )
	c_point_2d m_node_point; // cell node
	vector_block< double >		m_state;
	vector_block< double >		m_state_old;
	vector_block< double >		m_resid;
	vector_block< c_vector_2d >	m_grad;
	vector_block< c_matrix_2d >	m_hess;

	vector_block< double >		m_coef_rhs;
	matrix_block< double >		m_coef_slf;
	std::vector< matrix_block< double > > m_coef_nbr;
	
	vector_block< c_matrix_2d >	m_ans_spc; // anisotropic spacing function

	double sqr( double x ) { return x*x; }	
};

} //--- namespace NBCFD ---

#endif // CLASS_NIFS_C_CELL_H

