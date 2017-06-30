// Incompressible Flow Solver

// Copyright (C) 2005  Amir R. Baserinia

//------------------------------------------------------------------------------
// Date:   08 Apr 2004
// Update: 25 Feb 2005
//------------------------------------------------------------------------------

#ifndef CLASS_NIFS_C_FACE_H
#define CLASS_NIFS_C_FACE_H

#include "mesh_obj.h"
#include "point_2d.h"
#include "vector_2d.h"
#include "matrix_2d.h"
#include "vector_block.h"
#include "matrix_block.h"

#define	FACE_INTERNAL 0
#define FACE_BOUNDARY 1

namespace NIFS{    

class c_face : public c_mesh_obj {
	public:
		c_face();
		virtual ~c_face() {}

		const vector_block< double >& get_state( void ) { return m_state; }
		const vector_block< double >& get_error( void ) { return m_error; }
		const vector_block< c_vector_2d >& get_grad( void ) { return m_grad; } 
		const vector_block< c_matrix_2d >& get_hess( void ) { return m_hess; } 
		
		const vector_block< double >& get_weight( void ) { return m_weight; }
		const matrix_block< double >& get_coef_1( void ) { return m_coef_1; }
		const matrix_block< double >& get_coef_2( void ) { return m_coef_2; }
		const vector_block< double >& get_coef_c( void ) { return m_coef_c; }
		const vector_block< double >& get_flow( void )   { return m_flow; }

		c_point_2d get_fmp( void ) { return m_fmp; }
		double get_area() { return m_area; }
		const c_vector_2d& get_unv() { return m_unv; }
		const c_vector_2d& get_utv() { return m_utv; }
		double get_mfr( void ) { return m_mfr; }
		virtual unsigned get_type( void ) = 0;
		virtual double project( const c_vector_2d& pnt, c_vector_2d& prj ) = 0;
		
		virtual void update_geometry( void );
		virtual void update_solution( void );
		
		c_vector_2d get_collapse_pos( void );
		double get_size( void );
		double get_df() { return m_df;}

	protected:
		void calc_area();				// face area
		void calc_utv();				// face unit tangent vector
		void calc_unv();				// face unit normal vector
		void calc_fmp();				// calc face mid-point
		void calc_alpha();				// correction due to mesh anisotropy
		void calc_r1();					// vector from node 1 to mid-point
		virtual void calc_r2() = 0;		// vector from node 2 to mid-point
		virtual void calc_rc() = 0;		// vector from node 1&2 m-p to face m-p 
		virtual void calc_uv12() = 0;	// unit vector from node 1 to 2
		virtual void calc_s12() = 0;	// distance between node 1 and 2
		virtual void calc_coef_1() = 0;
		virtual void calc_coef_2() = 0;
		virtual void calc_coef_c() = 0;
		virtual void calc_mfr() = 0;	// mass flow rate
		virtual void calc_df() = 0;		// Rhie_Chow pressure correction 
		virtual void calc_state() = 0;
		virtual void calc_error() = 0;	// flow error

	protected:	
		double m_area;		// face area per unit depth
		double m_s12;		// distance between two neighbouring nodes	
		double m_alpha;		// non-orthogonality factor	\alpha = n_hat \dot s_hat 
		double m_df;		// Rhie-Chow pressure correction factor
		double m_mfr;		// mass flow rate

		c_point_2d m_fmp;	// face midpoint

		c_vector_2d m_unv;	// unit normal vector
		c_vector_2d m_utv;	// unit tangential vector
		c_vector_2d m_uv12;	// unit vector from node 1 to node 2 
		c_vector_2d m_r1;	// Vector joining the face mid-point to node 1
		c_vector_2d m_r2;	// Vector joining the face mid-point to node 2
		c_vector_2d m_rc;	// Vector joining the face mid-point to node 1&2 mid-point

		vector_block< double >		m_weight;
		vector_block< double >		m_state;
		vector_block< c_vector_2d >	m_grad;
		vector_block< c_matrix_2d >	m_hess;
		matrix_block< double >		m_coef_1;
		matrix_block< double >		m_coef_2;
		vector_block< double >		m_coef_c;
		vector_block< double >		m_flow;
		vector_block< double >		m_error;

	// | mass  flow |   |        | | p1 |   |        | | p2 |   | p_c |
	// | mom-u flow | = | coef_1 | | u1 | + | coef_2 | | u2 | + | u_c |
	// | mom-v flow |   |        | | v1 |   |        | | v2 |   | v_c |
	// | heat  flow |   |        | | T1 |   |        | | T2 |   | T_c |
};

} //--- namespace NIFS

#endif //--- CLASS_NIFS_C_FACE_H ---

