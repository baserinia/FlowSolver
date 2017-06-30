// Incompressible Flow Solver

// Copyright (C) 2005  Amir R. Baserinia

//------------------------------------------------------------------------------
// Created on:	09 Apr 2004
// Last update:	27 Jan 2005
//------------------------------------------------------------------------------

#ifndef CLASS_NIFS_C_MESH_H
#define CLASS_NIFS_C_MESH_H

#include <string>
#include <list>
#include "cell.h"
#include "face.h"
#include "vertex.h"
#include "parser.h"
#include "info.h"
#include "matrix_2d.h"

namespace NIFS{    

//--- INTERFACE ---

class c_mesh {
	public:
		c_mesh() {}
		c_mesh( c_mesh& mesh ); // copy constructor
		virtual ~c_mesh() { release_mesh(); }
			 
   	int load_mesh( std::string file_name );
 		void load_initial_sol( std::string file_name );    	
   	int save_mesh();

		void reset_vert_index( void );
		void reset_cell_index( void );
		void reset_face_index( void );
	
		void push_vert( c_vertex* vert )	{ m_verts.push_back( vert ); }
		void push_face( c_face*   face )	{ m_faces.push_back( face ); }
		void push_cell( c_cell*   cell )	{ m_cells.push_back( cell ); }
	
		int get_cell_num() { return m_cells.size(); }
		int get_vert_num() { return m_verts.size(); }
		int get_face_num() { return m_faces.size(); }
	
		void release_verts( void );
		void release_faces( void );
		void release_cells( void );
		void release_mesh( void );
	
		c_cell* get_first_cell( void );
		c_cell* get_next_cell( void );
		
		void update_solution( void );
		void update_geometry( void );
		
		vector_block< double > get_residual_rms( void );
		vector_block< double > get_residual_max( void );

		int get_coef_num() { return m_coef_num; }
		void print_info( void );
		// remove these
		void   load_size( void );
		double get_size( c_vector_2d vec );
		c_matrix_2d get_size_matrix( c_vector_2d vec );
    double get_length_ratio( c_vector_2d , c_vector_2d );
    double get_length_ratio_move( c_vector_2d , c_vector_2d );				
		
		void operator=( c_mesh& mesh ); 	// overloaded assignment operator	
		c_cell* search( c_vector_2d vec );
		void reset_current_cell( void ) { m_current_cell = get_first_cell(); } 
		vector_block< double > get_avg_err( void ) { return m_avg_err; }

	protected:
		double sqr( double x ) { return x * x; }
		std::list< c_vertex* >	m_verts; 
		std::list< c_face* >	m_faces; 
		std::list< c_cell* >	m_cells; 
		std::list< c_vertex* >::iterator	m_vert_iter; 
		std::list< c_face* >::iterator		m_face_iter; 
		std::list< c_cell* >::iterator		m_cell_iter; 
		c_cell* m_current_cell;
		c_info m_info;
		c_parser parser;
		int m_coef_num;
		double E[256][256];
		double A11[256][256];
		double A12[256][256];
		double A22[256][256];
		vector_block<double> m_avg_err;
};

} // namespace NIFS

#endif // CLASS_C_MESH_H
