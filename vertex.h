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
// Created on:	08 Apr 2004
// Last update:	07 Feb 2006
//------------------------------------------------------------------------------

#ifndef CLASS_NIFS_C_VERTEX_H
#define CLASS_NIFS_C_VERTEX_H

#include "mesh_obj.h"
#include "point_2d.h"
#include "vector_2d.h"

namespace NIFS{    

//--- INTERFACE ---

class c_vertex: public c_mesh_obj
{
public:
	c_vertex() { m_active = 0; }
	c_vertex( const c_point_2d& point ) { m_point = point; m_active = 0; }
	c_vertex( double x, double y ) { m_point.set_pos( x, y ); m_active = 0; }
	virtual ~c_vertex() {}
	const c_point_2d& get_point() const { return m_point; }
	void set_point( const c_point_2d& point ) { m_point = point; }
	void set_pos( const c_vector_2d& vec ) { m_point.set_pos( vec ); }
	virtual unsigned get_type( void ); // pure virtual -> later
	c_vector_2d get_boundary_tangent( void );
	c_vector_2d move_vertex_bnd( c_vector_2d vec );
	void move_vertex_rel( c_vector_2d vec ) {}
	void move_vertex_to( c_vector_2d vec ) {}
	
	void set_active( int active ) { m_active = active; }
	int get_active( void ) { return m_active; }

private:
	int m_type;				// type of the vertex
	c_point_2d m_point; 	// point object of vertex
	int m_active;
};

} //--- namespace NIFS ---

#endif //--- CLASS_NIFS_C_VERTEX_H ---

