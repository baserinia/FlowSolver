// Incompressible Flow Solver

// Copyright (C) 2005  Amir R. Baserinia

//----------------------------------------------------------------------------
// Created on:	08 Apr 2004
// Last update:	28 Mar 2005
//-----------------------------------------------------------------------------

#ifndef CLASS_NIFS_C_MESH_OBJ_H
#define CLASS_NIFS_C_MESH_OBJ_H

#include "list_nbr.h"
#include "info.h"

namespace NIFS{    

//--- INTERFACE ---

class c_mesh_obj
{
	public:
		c_mesh_obj() {}
		virtual ~c_mesh_obj() {}
	
		unsigned int get_vert_num() { return m_verts_nbr.size(); }
		unsigned int get_face_num() { return m_faces_nbr.size(); }        
		unsigned int get_cell_num() { return m_cells_nbr.size(); }    

		bool is_vert_nbr( c_mesh_obj* vert ) { return m_verts_nbr.has_obj( vert ); }
		bool is_face_nbr( c_mesh_obj* face ) { return m_faces_nbr.has_obj( face ); }
		bool is_cell_nbr( c_mesh_obj* cell ) { return m_cells_nbr.has_obj( cell ); }

		void push_vert( c_mesh_obj* vert ) { m_verts_nbr.push_back( vert ); }
		void push_face( c_mesh_obj* face ) { m_faces_nbr.push_back( face ); }    
		void push_cell( c_mesh_obj* cell ) { m_cells_nbr.push_back( cell ); }    

		void remove_vert( c_mesh_obj* vert ) { m_verts_nbr.remove( vert ); }
		void remove_face( c_mesh_obj* face ) { m_faces_nbr.remove( face ); }
		void remove_cell( c_mesh_obj* cell ) { m_cells_nbr.remove( cell ); }

		void clear_verts() { m_verts_nbr.clear( ); }
		void clear_faces() { m_faces_nbr.clear( ); }
		void clear_cells() { m_cells_nbr.clear( ); }
		void clear_all() { clear_verts(); clear_faces(); clear_cells(); }

		c_mesh_obj* get_vert( unsigned index ) { return m_verts_nbr[ index ]; }
		c_mesh_obj* get_face( unsigned index ) { return m_faces_nbr[ index ]; }
		c_mesh_obj* get_cell( unsigned index ) { return m_cells_nbr[ index ]; }

		void set_vert( unsigned index , c_mesh_obj* obj ) { m_verts_nbr.set( index , obj ); }
		void set_face( unsigned index , c_mesh_obj* obj ) { m_faces_nbr.set( index , obj ); }
		void set_cell( unsigned index , c_mesh_obj* obj ) { m_cells_nbr.set( index , obj ); }

		
		void set_index( unsigned int index ) { m_index = index; }
		void set_code( unsigned int code ) { m_code = code; }
		unsigned int get_index( void ) { return m_index; }
		unsigned int get_code( void ) { return m_code; }
		const c_info& get_info() { return m_info; }

	private:
		list_nbr<c_mesh_obj> m_verts_nbr; // Neighbouring Vertices
		list_nbr<c_mesh_obj> m_faces_nbr; // Neighbouring Faces
		list_nbr<c_mesh_obj> m_cells_nbr; // Neighbouring Cells
		unsigned int m_index;
		unsigned int m_code;
		c_info m_info;
};

} // namespace NIFS

#endif //--- CLASS_NIFS_C_MESH_OBJ_H

