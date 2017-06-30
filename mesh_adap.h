//	Copyright (C) 2006  Amir R. Baserinia
//
// 	This file is part of SMAIF (Simple Mesh Adaptor for Incompressible Flows)
//
//	SMAIF is free software; you can redistribute it and/or modify
//	it under the terms of the GNU General Public License as published by
//	the Free Software Foundation; either version 2 of the License, or
//	(at your option) any later version.
//
//	SMAIF is distributed in the hope that it will be useful,
//	but WITHOUT ANY WARRANTY; without even the implied warranty of
//	MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//	GNU General Public License for more details.
//
//	You should have received a copy of the GNU General Public License
//	along with SMAIF; if not, write to the Free Software
//	Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA

//------------------------------------------------------------------------------
// Created on:	17 Sep 2004
// Last update:	17 Oct 2005
//------------------------------------------------------------------------------

#ifndef CLASS_NIFS_C_MESH_ADAP_H
#define CLASS_NIFS_C_MESH_ADAP_H

#include "mesh.h"

namespace NIFS{    

//--- INTERFACE ---

class c_mesh_adap : public c_mesh  
{
public:
	c_mesh_adap() {}
	c_mesh_adap( c_mesh& mesh ) : c_mesh::c_mesh( mesh ) {} 
	~c_mesh_adap() {}
	int split_face( c_face* face );
	int collapse_face( c_face* face );
	int collapse_face2( c_face* face );
	std::list< c_face* >::iterator collapse_face( std::list<c_face*>::iterator face_iter );
	void remove_vert( c_vertex* );
	void set_bg_mesh( c_mesh* bg_mesh ) { m_bg_mesh = bg_mesh; }
	void refine( void );
	void modify( void );
	int swap_face( c_face* face );
	int	 delaunay( c_face* face );
	int	 delaunay2( c_face* face );
	void delaunay( c_vertex* );
	bool delaunay( void );
//	void smooth( void );
	double smooth( void );	
	double smooth_vertex( c_vertex* );
  int merge_tri_to_quad( c_face* face );
  int split_quad_to_tri( c_cell* cell );		
	void adapt( void );	
	void switch_quad( void );
//	void switch_quad2( void );
//  void switch_quad3( void );
//  void switch_quad4( void );
	int switch_quad_nbr( c_cell* cell );
  void cleanup_tri( void );
  void pushout_tri( void );  
  int enhance_topolog( void );
  int enhance_topolog( c_face* face );
  int enhance_topolog( c_vertex* vert );

private:
	double pi() { return 3.14159; }
	c_mesh* m_bg_mesh;	// back ground mesh
};

} // namespace NBCFD

#endif // CLASS_NIFS_C_MESH_ADAP_H
