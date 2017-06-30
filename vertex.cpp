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

//-----------------------------------------------------------------------------
// Date:   15 Feb 2004
// Update: 15 Feb 2005
//-----------------------------------------------------------------------------

#include <iostream>
#include "vertex.h"
#include "face.h"

//--- IMPLEMENTATION ---


//	This function returns the type of the boundary vertex
//	The returned value is:
//		0: internal 
//		1: boundary
//		2: corner
unsigned
NIFS::c_vertex::get_type( void )
{
	c_face* face_ptr;
	unsigned bc_code_1 = 0;
	unsigned bc_code_2 = 0;
	
	for ( unsigned i = 0 ; i < get_face_num( ) ; ++i )
	{
		face_ptr = static_cast< c_face* >( get_face( i ) );
		if ( face_ptr->get_type() != 0 ) // non-internal face
			if ( bc_code_1 == 0 )
				bc_code_1 = face_ptr->get_code();
			else
			{
				bc_code_2 = face_ptr->get_code();
				break;
			}
	}

	if ( bc_code_1 == 0 )	
		return 0;	// internal vertex
	else if ( bc_code_1 == bc_code_2 )	
		return 1;	// non-corner boundary vertex
	else 
		return 2;	// corner boundary vertex
}


NIFS::c_vector_2d 
NIFS::c_vertex::get_boundary_tangent( void )
{
	c_vertex*	vert_ptr;
	c_vertex*	vert_ptr_1 = NULL;
	c_vertex*	vert_ptr_2 = NULL;
	
	if ( get_type() != 1 )		// internal or corner vertex
		return c_vector_2d( 0.0 , 0.0 );	
	else						// non-corner boundary vertex
	{
		for ( unsigned i = 0 ; i < get_vert_num( ) ; ++i )
		{
			vert_ptr = static_cast< c_vertex* >( get_vert( i ) );
			if ( vert_ptr->get_type() != 0 ) // nin-internal face
			{
				if ( vert_ptr_1 == NULL )
					vert_ptr_1 = vert_ptr;
				else
					vert_ptr_2 = vert_ptr;
			}
		}
		return ( vert_ptr_1->get_point().get_pos() - vert_ptr_2->get_point().get_pos() ).get_unit();
	}
}


NIFS::c_vector_2d
NIFS::c_vertex::move_vertex_bnd( c_vector_2d vec )
{
	c_face* face_ptr;
	c_face* face_ptr_1 = NULL;
	c_face* face_ptr_2 = NULL;
	c_vector_2d prj_1;
	c_vector_2d prj_2;

	for ( unsigned i = 0 ; i < get_face_num( ) ; ++i )
	{
		face_ptr = static_cast< c_face* >( get_face( i ) );
		if ( face_ptr->get_type() != 0 ) // non-internal face
		{
			if ( face_ptr_1 == NULL )
				face_ptr_1 = face_ptr;
			else
				face_ptr_2 = face_ptr;
		}
	}
	
	double t1 = face_ptr_1->project( vec , prj_1 );
	double t2 = face_ptr_2->project( vec , prj_2 );
		
	bool cnd1 = ( ( 0.0 < t1 ) && ( t1 < 1.0 ) );
	bool cnd2 = ( ( 0.0 < t2 ) && ( t2 < 1.0 ) );

	if ( cnd1 && !cnd2 )
		return prj_1;
	else if ( !cnd1 && cnd2 )
		return prj_2;
  else
    return vec;	
}
