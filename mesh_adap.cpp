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

//----------------------------------------------------------------------------
// Created on:	17 Sep 2004
// Last update:	13 Aug 2006
//-----------------------------------------------------------------------------

#include <vector>
#include <list>
#include <algorithm>
#include <cstdlib>
#include "mesh_adap.h"
#include "vector_2d.h"
#include "matrix_2d.h"
#include "point_2d.h"
#include <iostream>
#include "face_bnd.h"
#include "face_int.h"
#include "face_bnd_wall.h"
#include "face_bnd_wall_dir.h"
#include "face_bnd_wall_neu.h"
#include "face_bnd_wall_dir_pres.h"
#include "face_bnd_inflow.h"
#include "face_bnd_outflow.h"
#include "face_bnd_sym.h"
#include "my_except.h"

//--- IMPLEMENTATION ---


//---( Split a face into two sub-faces )------------------------------------------
// This method splits a face into two sub-faces by inserting a node at the 
// face mid-point and modifying the local topology.
//
// note: this method is general, but it will be used for triangular meshes only
int
NIFS::c_mesh_adap::split_face( c_face* face )
{
//  short face_type;
  unsigned i;
  unsigned j;

	for ( i = 0 ; i < face->get_cell_num(); ++i ) {
		// for triangular mesh
		if ( static_cast<c_cell*>( face->get_cell(i) )->get_face_num() != 3 ) 
			return -1; // at least one of the neighbours is non-triangular
	}


	// consruct the new vertex at the face mid-point and push it into vertex list
	c_vector_2d vec;
	if ( face->get_type() == 0 )
		vec = face->get_fmp().get_pos();
	else
		face->project( face->get_fmp().get_pos() , vec );
		
	c_vertex* vert = new c_vertex( vec );
	m_verts.push_back( vert );
	
  if ( face->get_type() == 0 ) {	// if face is internal face 
	// create a temporary vector of neighbouring cells
	unsigned cell_num =	static_cast< c_cell* >( face->get_cell( 0 ) )->get_face_num() +
	                    static_cast< c_cell* >( face->get_cell( 1 ) )->get_face_num() - 2;
	std::vector< c_cell* > cell_temp( cell_num );
	cell_temp[0] = static_cast< c_cell* >( face->get_cell( 0 ) );
	cell_temp[1] = static_cast< c_cell* >( face->get_cell( 1 ) );
	for ( i = 2 ; i < cell_num ; ++i) {
	
	  cell_temp[i] = new c_cell( );      // initialize cell and push it into temp vector
	  m_cells.push_back( cell_temp[i] ); // push new cells into cell list
	}
	
	// create temporary faces
	unsigned face_num =	static_cast< c_cell* >( face->get_cell( 0 ) )->get_face_num() +
	                    static_cast< c_cell* >( face->get_cell( 1 ) )->get_face_num() - 2;
	std::vector< c_face* > face_temp( face_num );
	face_temp[0] = face; // the first face is the existing one
	for ( i = 1 ; i < face_num ; ++i ) {
	  face_temp[i] = new c_face_int();	// initialize face and it into temp ventor
	  m_faces.push_back( face_temp[i] );// push new faces into the face list
	}
	
	// faces of cavity associated with the removal of face
	std::vector< c_face* > face_bnd_temp( face_num );	
	j = 0;
	for ( i = 0 ; i < cell_temp[0]->get_face_num(); ++i ) {
	  if ( cell_temp[0]->get_face( i ) != face ) {
		face_bnd_temp[j] =  static_cast< c_face* >( cell_temp[0]->get_face( i ) ) ;
		++j;
	  }
	}
	for ( i = 0 ; i < cell_temp[1]->get_face_num(); ++i ) {
	  if ( cell_temp[1]->get_face( i ) != face ) {
		face_bnd_temp[j] =  static_cast< c_face* >( cell_temp[1]->get_face( i ) ) ;
		++j;
	  }
	}

	// vertices of cavity associated with the removal of face
	unsigned vert_num = face_num;
	std::vector< c_vertex* > vert_temp( vert_num ); 

	vert_temp[0] = static_cast< c_vertex*>( face->get_vert( 0 ) );
	vert_temp[1] = static_cast< c_vertex*>( face->get_vert( 1 ) );
	j = 2;
	for ( i = 0 ; i < cell_temp[0]->get_vert_num(); ++i ) {
	  if ( ( cell_temp[0]->get_vert( i ) != vert_temp[0] ) &&
		   ( cell_temp[0]->get_vert( i ) != vert_temp[1] ) ) {
		vert_temp[j] =  static_cast< c_vertex* >( cell_temp[0]->get_vert( i ) ) ;
		++j;
	  }
	}
	for ( i = 0 ; i < cell_temp[1]->get_vert_num(); ++i ) {
    if ( ( cell_temp[1]->get_vert( i ) != vert_temp[0] ) &&
         ( cell_temp[1]->get_vert( i ) != vert_temp[1] ) ) {
		vert_temp[j] =  static_cast< c_vertex* >( cell_temp[1]->get_vert( i ) ) ;
		++j;
	  }
	}

	for ( i = 0 ; i < vert_num ; ++i ) {
	  // all the vertices of the cavity
	  if ( vert_temp[i]->is_face_nbr( face ) ) {
		// the end-vertices of the face being removed
		if ( vert_temp[i] == face->get_vert( 0 ) )
		  vert_temp[i]->remove_vert( face->get_vert(1) );
		if ( vert_temp[i] == face->get_vert( 1 ) )
		  vert_temp[i]->remove_vert( face->get_vert( 0 ) );
		vert_temp[i]->remove_face( face );
		vert_temp[i]->remove_cell( face->get_cell( 0 ) );
		vert_temp[i]->remove_cell( face->get_cell( 1 ) );
	  } 
	  else {
		// vertices other than the end-vertices of the face
		if ( vert_temp[i]->is_cell_nbr( face->get_cell( 0 ) ) )
		  vert_temp[i]->remove_cell( face->get_cell(0) );
		if ( vert_temp[i]->is_cell_nbr( face->get_cell( 1 ) ) )
		  vert_temp[i]->remove_cell( face->get_cell(1) );
	  }
	}
	// clear all the connectivities associated with the face
//	std::cout << face->get_cell_num() << std::endl;
	face->clear_all( );

	for ( i = 0 ; i < face_num ; ++i ) {
	  // sweep all the boundary faces of the cavity
	  if ( ( face_bnd_temp[i]->get_cell( 0 ) == cell_temp[0] ) ||
		   ( face_bnd_temp[i]->get_cell( 0 ) == cell_temp[1] ) ) {
		// if the neighbour 0 is inside the cavity
		if ( face_bnd_temp[i]->get_cell_num() > 1 ) 
		  // if the face is interior to the domain
		  static_cast< c_cell* >( face_bnd_temp[i]->get_cell( 1 ) )->remove_cell( face_bnd_temp[i]->get_cell( 0 ) );
		// set the new neibour of the cavity boundary face
		face_bnd_temp[i]->set_cell( 0 , cell_temp[i] );
	  }
	  if ( ( face_bnd_temp[i]->get_cell( 1 ) == cell_temp[0] ) ||
		   ( face_bnd_temp[i]->get_cell( 1 ) == cell_temp[1] ) ) {
		// if the neibour 1 is inside the cavity
		if ( face_bnd_temp[i]->get_cell_num() > 1 ) 
		  // if the face is interior to the domain 
		  static_cast< c_cell* >( face_bnd_temp[i]->get_cell( 0 ) )->remove_cell( face_bnd_temp[i]->get_cell( 1 ) );
		// set the new neibour of the cavity boundary face
		face_bnd_temp[i]->set_cell( 1 , cell_temp[i] );
	  }
	}

//	std::cout << "step" << std::endl;
	// note: face_num == cell_num
	// I assumed "face_temp[i]" is adjuscent to "cell_temp[i]"
	for ( i = 0 ; i < face_num ; ++i ) {
	  // clear all the neighbours of the new cells
//	std::cout << cell_temp[i] << std::endl;			  
	  cell_temp[i]->clear_all( );

	  // the new inserted vertex and the end vertices of face_temp[i] are the new
	  // neighbours of the control volume temp_cell[i]
	  cell_temp[i]->push_vert( static_cast< c_vertex* >( face_bnd_temp[i]->get_vert( 0 ) ) );
	  cell_temp[i]->push_vert( static_cast< c_vertex* >( face_bnd_temp[i]->get_vert( 1 ) ) );
	  cell_temp[i]->push_vert( vert );	// vert is common among all the new cells

	  // face_temp[i] is also a neighbour of the new cell
	  cell_temp[i]->push_face( face_bnd_temp[i] );
//		std::cout << face_num << std::endl;		  
	}

	// note: face_num == vert_num
	// I assumed "vert_temp[i]" is the end vertex of "face_temp[i]"	
	for ( i = 0 ; i < vert_num ; ++i ) {
	  face_temp[i]->clear_verts( );
	  
	  face_temp[i]->push_vert( vert );
	  face_temp[i]->push_vert( vert_temp[i] );
	  
	  for ( j = 0 ; j < face_num ; ++j ) 
		if ( vert_temp[i]->is_face_nbr( face_bnd_temp[j] ) ) {
		  face_temp[i]->push_cell( cell_temp[j] );
		  cell_temp[j]->push_face( face_temp[i] );
		}
	}

	for ( i = 0 ; i < face_num ; ++i ) {
	  static_cast<c_cell*>( face_temp[i]->get_cell( 0 ) )->push_cell( face_temp[i]->get_cell( 1 ) );
	  static_cast<c_cell*>( face_temp[i]->get_cell( 1 ) )->push_cell( face_temp[i]->get_cell( 0 ) );
	  
	  if (  face_bnd_temp[i]->get_cell_num() > 1 ) { 
		  static_cast<c_cell*>( face_bnd_temp[i]->get_cell(0) )->push_cell( face_bnd_temp[i]->get_cell(1) );
		  static_cast<c_cell*>( face_bnd_temp[i]->get_cell(1) )->push_cell( face_bnd_temp[i]->get_cell(0) );
	  }
	}

	for ( i = 0 ; i < face_num ; ++i ) {
	  vert_temp[i]->push_vert( vert );
	  vert->push_vert( vert_temp[i] );
	  vert_temp[i]->push_face( face_temp[i] );
	  vert->push_face( face_temp[i] );
	  vert_temp[i]->push_cell( face_temp[i]->get_cell(0) );
	  vert_temp[i]->push_cell( face_temp[i]->get_cell(1) );
	  vert->push_cell( cell_temp[i] );
	}
	

	
	for ( i=0; i<cell_num; i++)
	{
	  cell_temp[i]->update_geometry();
//	  std::cout << cell_temp[i]->get_volume() << std::endl;
	}	  	
	
	for ( i=0; i<cell_num; i++)
	  face_bnd_temp[i]->update_geometry();
	
	for ( i=0; i<face_num; i++){
	  face_temp[i]->update_geometry();

	}

  } 
  
  else {  // BOUNDARY FACE
	
	// create temporary cells
	unsigned cell_num =	static_cast<c_cell*>( face->get_cell( 0 ) )->get_face_num()-1;
	std::vector<c_cell*> cell_temp( cell_num );
	cell_temp[0] = static_cast<c_cell*>( face->get_cell( 0 ) );
	for ( i = 1 ; i < cell_num ; ++i ) {
	  cell_temp[i] = new c_cell( );
	  m_cells.push_back( cell_temp[i] );
	}


	// create temporary faces
	unsigned face_num =	static_cast<c_cell*>( face->get_cell( 0 ) )->get_vert_num();
	std::vector<c_vertex*> vert_temp( face_num );	// cavity boundary vertices
	for ( i = 0 ; i < face_num ; ++i ) {
	  vert_temp[i] = static_cast<c_vertex*>( cell_temp[0]->get_vert( i ) );
	}

	std::vector<c_face*> face_temp( face_num );
//	face_temp[0] = face;
	for ( i = 0 ; i < face_num ; ++i ) {
	  if ( ( vert_temp[i] != face->get_vert( 0 ) ) &&
	  	   ( vert_temp[i] != face->get_vert( 1 ) ) )
//	if ( vert_temp[i]->get_type() == 0  )
	  {
	    face_temp[i] = new c_face_int();
	    m_faces.push_back( face_temp[i] );
	  }
	  else if ( vert_temp[i] == face->get_vert( 0 ) )
	  	face_temp[i] = face;
	  else {
	    switch ( m_info.get_bnd_cond()->get_bc_type( face->get_code() ) ) {     
		case 1:
			face_temp[i] = new c_face_bnd_inflow();		// in-flow
			break;
		case 2:
			face_temp[i] = new c_face_bnd_outflow();	// out-flow
			break;
		case 3:
			face_temp[i] = new c_face_bnd_sym();		// symmetry
			break;
		case 4:
			face_temp[i] = new c_face_bnd_wall_dir( );	// wall dirichlet
			break;
		case 5:
			face_temp[i] = new c_face_bnd_wall_neu( );	// wall neumann
			break;
		case 6:
			face_temp[i] = new c_face_bnd_wall_dir_pres();	// wall dirichlet & pres
			break;
        default:
			throw NIFS::my_except( "unidentified boundary condition type" );
		} // switch
		face_temp[i]->set_code( face->get_code() );
		m_faces.push_back( face_temp[i] );
	  }
	}
		
	std::vector<c_face*> face_bnd_temp( cell_num );	// cavity boundary faces
	j = 0;
	for (i=0; i<cell_temp[0]->get_face_num(); i++) {
	  if ( cell_temp[0]->get_face( i ) != face ) {
		face_bnd_temp[j] =  static_cast<c_face*>( cell_temp[0]->get_face( i ) ) ;
		j++;
	  }
	}

	
	for (i=0; i<face_num; i++) {
	  if ( vert_temp[i]->is_face_nbr( face ) ) {
		if ( vert_temp[i] == face->get_vert(0) )
		  vert_temp[i]->remove_vert( face->get_vert(1) );
		if ( vert_temp[i] == face->get_vert(1) )
		  vert_temp[i]->remove_vert( face->get_vert(0) );
		vert_temp[i]->remove_face( face );
		vert_temp[i]->remove_cell( face->get_cell(0) );
	  } else {
		if ( vert_temp[i]->is_cell_nbr( face->get_cell(0) ) )
		  vert_temp[i]->remove_cell( face->get_cell(0) );
	  }
	}
	face->clear_all();
	
	for (i=0; i<cell_num; i++) {
	  if ( face_bnd_temp[i]->get_cell(0) == cell_temp[0] )  {
		if ( face_bnd_temp[i]->get_cell_num() > 1) 
		  static_cast<c_cell*>( face_bnd_temp[i]->get_cell(1) )->remove_cell( face_bnd_temp[i]->get_cell(0) );
		face_bnd_temp[i]->set_cell( 0, cell_temp[i] );
	  }
	  if ( ( face_bnd_temp[i]->get_cell_num() > 1) &&
		   ( face_bnd_temp[i]->get_cell(1) == cell_temp[0] ) ) {
		if ( face_bnd_temp[i]->get_cell_num() > 1) 
		  static_cast<c_cell*>( face_bnd_temp[i]->get_cell(0) )->remove_cell( face_bnd_temp[i]->get_cell(1) );
		face_bnd_temp[i]->set_cell( 1, cell_temp[i] );
	  }
	}
	
	for (i=0; i<cell_num; i++) {
	  cell_temp[i]->clear_all();
	  
	  cell_temp[i]->push_vert( static_cast<c_vertex*>( face_bnd_temp[i]->get_vert(0) ) );
	  cell_temp[i]->push_vert( static_cast<c_vertex*>( face_bnd_temp[i]->get_vert(1) ) );
	  cell_temp[i]->push_vert( vert );
	  
	  cell_temp[i]->push_face( face_bnd_temp[i] );
	}
	
	for (i=0; i<face_num; i++) {
	  
	  face_temp[i]->clear_verts();
	  face_temp[i]->push_vert( vert );
	  face_temp[i]->push_vert( vert_temp[i] );
	  
	  for (j=0; j<cell_num; j++) 
		if ( vert_temp[i]->is_face_nbr( face_bnd_temp[j] ) ) {
		  face_temp[i]->push_cell( cell_temp[j] );
		  cell_temp[j]->push_face( face_temp[i] );
		}
	}

	for (i=0; i<face_num; i++) 
	  if ( face_temp[i]->get_cell_num() > 1) {
		static_cast<c_cell*>( face_temp[i]->get_cell(0) )->push_cell( face_temp[i]->get_cell(1) );
		static_cast<c_cell*>( face_temp[i]->get_cell(1) )->push_cell( face_temp[i]->get_cell(0) );
	  }
	
	for (i=0; i<cell_num; i++) 
	  if (  face_bnd_temp[i]->get_cell_num() > 1 ) { 
		static_cast<c_cell*>( face_bnd_temp[i]->get_cell(0) )->push_cell( face_bnd_temp[i]->get_cell(1) );
		static_cast<c_cell*>( face_bnd_temp[i]->get_cell(1) )->push_cell( face_bnd_temp[i]->get_cell(0) );
	  }

	for (i=0; i<face_num; i++) {
	  vert_temp[i]->push_vert( vert );
	  vert->push_vert( vert_temp[i] );
	  vert_temp[i]->push_face( face_temp[i] );
	  vert->push_face( face_temp[i] );
	  vert_temp[i]->push_cell( face_temp[i]->get_cell(0) );
	  if ( face_temp[i]->get_cell_num() > 1)
		vert_temp[i]->push_cell( face_temp[i]->get_cell(1) );
	}
	vert->push_cell( cell_temp[0] );
	vert->push_cell( cell_temp[1] );
	
	for ( i=0; i<cell_num; i++)
	  cell_temp[i]->update_geometry();
	
	for ( i=0; i<cell_num; i++)
	  face_bnd_temp[i]->update_geometry();
	
	for ( i=0; i<face_num; i++)
	  face_temp[i]->update_geometry();

  }

	return 0;
	
} // end of "split_face" method




//-------------------------------------
//	Collapse Face						  
//	This method removes a face and brings the end points together
//	At the moment, it is very slow. The design of the mesh lists
//	must be modified to improve the performance of this method.
//
//int
//NIFS::c_mesh_adap::collapse_face( c_face* face )
std::list<NIFS::c_face*>::iterator 
NIFS::c_mesh_adap::collapse_face( std::list<c_face*>::iterator   face_iter)
{

	// special treatment for trangles
	unsigned i;
	unsigned j;
	c_face* face = *face_iter;
	c_cell* cell_ptr; //cell_ptr
	c_vertex* vert_ptr_0 = static_cast<c_vertex*>( face->get_vert(0) );
	c_vertex* vert_ptr_1 = static_cast<c_vertex*>( face->get_vert(1) );
	c_face* face_ptr_0;
	c_face* face_ptr_1;
	c_cell* cell_ptr_0;
	c_cell* cell_ptr_1;
	std::list<NIFS::c_face*>::iterator face_iter_ret = face_iter;
	face_iter_ret++;

//	c_vector_2d new_pos = face->get_collapse_pos();
	
	if ( ( vert_ptr_0->get_type() == 2 ) && 
		 ( vert_ptr_1->get_type() == 2 ) ) 
//		return -1; // both face end-points are corner vertices
		return (face_iter_ret);


	if ( ( vert_ptr_0->get_type() != 0 ) && 
		 ( vert_ptr_1->get_type() != 0 ) &&
		 ( face->get_cell_num() == 2  ) )
//		return -2; // face joins 2 boundary points and cannot be removed
		return (face_iter_ret);		

	for ( i = 0 ; i < face->get_cell_num(); ++i ) {
		// for triangular mesh
		if ( static_cast<c_cell*>( face->get_cell(i) )->get_face_num() != 3 ) 
//			return -3; // at least one of the neighbours is non-triangular
			return (face_iter_ret);			
	}

	// vert_0 is retained and vert_1 is removed
	if ( vert_ptr_1->get_type() > vert_ptr_0->get_type() ) {
		c_vertex* vert_temp = vert_ptr_1;
		vert_ptr_1 = vert_ptr_0;
		vert_ptr_0 = vert_temp;
	}

	c_vector_2d new_pos = face->get_collapse_pos();
//	c_vector_2d new_pos = vert_ptr_0->get_point().get_pos();

	for ( i = 0 ; i < vert_ptr_0->get_cell_num() ; ++i ) {
		cell_ptr = static_cast<c_cell*>( vert_ptr_0->get_cell(i) );
		if ( ( cell_ptr != face->get_cell( 0 ) ) &&
			 ( cell_ptr != face->get_cell( 1 ) ) &&
			 ( cell_ptr->is_valid_move( vert_ptr_0, new_pos ) == false )
		   )
//			return -4;
			return (face_iter_ret);
	}

	for ( i = 0 ; i < vert_ptr_1->get_cell_num() ; ++i ) {
		cell_ptr = static_cast<c_cell*>( vert_ptr_1->get_cell(i) );
		if ( ( cell_ptr != face->get_cell( 0 ) ) &&
			 ( cell_ptr != face->get_cell( 1 ) ) &&
			 ( cell_ptr->is_valid_move( vert_ptr_1, new_pos ) == false )
		   )
//			return -4;
			return (face_iter_ret);
	}


	for ( i = 0 ; i < face->get_cell_num(); ++i ) {
		cell_ptr_0 = NULL;
		cell_ptr_1 = NULL;
		face_ptr_0 = NULL;
		face_ptr_1 = NULL;

		cell_ptr = static_cast<c_cell*>( face->get_cell(i) );
		for ( j = 0 ; j < cell_ptr->get_face_num(); ++j ) {
			if ( face != static_cast<c_face*>( cell_ptr->get_face(j) ) ) {
				if ( static_cast<c_face*>( cell_ptr->get_face(j) )->is_vert_nbr( vert_ptr_0 ) ) {
					face_ptr_0 = static_cast<c_face*>( cell_ptr->get_face(j) );
				} else {
					face_ptr_1 = static_cast<c_face*>( cell_ptr->get_face(j) );
				}
			}
		}
			
		if ( face_ptr_0->get_cell_num() > 1) {
			if ( face_ptr_0->get_cell(0) == cell_ptr ) {
				cell_ptr_0 = static_cast<c_cell*>( face_ptr_0->get_cell(1) );
			} else {
				cell_ptr_0 = static_cast<c_cell*>( face_ptr_0->get_cell(0) );
			}
		}

		if ( face_ptr_1->get_cell_num() > 1) {
			if ( face_ptr_1->get_cell(0) == cell_ptr ) {
				cell_ptr_1 = static_cast<c_cell*>( face_ptr_1->get_cell(1) );
			} else {
				cell_ptr_1 = static_cast<c_cell*>( face_ptr_1->get_cell(0) );
			}
		}

		if ( cell_ptr_0 != NULL) {
			cell_ptr_0->remove_cell( cell_ptr );
			if ( cell_ptr_1 != NULL ) 
				cell_ptr_0->push_cell( cell_ptr_1 );
		}

		if ( cell_ptr_1 != NULL) {
			cell_ptr_1->remove_cell( cell_ptr );
			cell_ptr_1->remove_face( face_ptr_1 );
			cell_ptr_1->push_face( face_ptr_0 );
			if ( cell_ptr_0 != NULL ) 
				cell_ptr_1->push_cell( cell_ptr_0 );
		}
			
		face_ptr_0->remove_cell( cell_ptr );
		if ( cell_ptr_1 != NULL )
			face_ptr_0->push_cell( cell_ptr_1 );
		vert_ptr_0->remove_cell( cell_ptr );
		vert_ptr_1->remove_cell( cell_ptr );
		vert_ptr_1->remove_face( face_ptr_1 );
		if ( vert_ptr_1 == static_cast<c_vertex*>( face_ptr_1->get_vert(0) ) ) {
			vert_ptr_1->remove_vert( static_cast<c_vertex*>( face_ptr_1->get_vert(1) ) );
			static_cast<c_vertex*>( face_ptr_1->get_vert(1) )->remove_vert( vert_ptr_1 );
			static_cast<c_vertex*>( face_ptr_1->get_vert(1) )->remove_face( face_ptr_1 );
			static_cast<c_vertex*>( face_ptr_1->get_vert(1) )->remove_cell( cell_ptr );
		} else {
			vert_ptr_1->remove_vert( static_cast<c_vertex*>( face_ptr_1->get_vert(0) ) );
			static_cast<c_vertex*>( face_ptr_1->get_vert(0) )->remove_vert( vert_ptr_1 );
			static_cast<c_vertex*>( face_ptr_1->get_vert(0) )->remove_face( face_ptr_1 );
			static_cast<c_vertex*>( face_ptr_1->get_vert(0) )->remove_cell( cell_ptr );
		}

		m_faces.remove( face_ptr_1 );
		delete face_ptr_1;
		face_ptr_1 = NULL;

		m_cells.remove( cell_ptr );
		delete cell_ptr;
		cell_ptr = NULL;
	}

		
	vert_ptr_1->remove_vert( vert_ptr_0 );
	vert_ptr_1->remove_face( face );
	vert_ptr_0->remove_vert( vert_ptr_1 );
	vert_ptr_0->remove_face( face );

	vert_ptr_0->set_pos( new_pos );
	
	
	for ( i = 0 ; i < vert_ptr_1->get_vert_num() ; ++i ) {
		vert_ptr_0->push_vert( static_cast<c_vertex*>( vert_ptr_1->get_vert(i) ) );
		static_cast<c_vertex*>( vert_ptr_1->get_vert(i) )->remove_vert( vert_ptr_1 );
		static_cast<c_vertex*>( vert_ptr_1->get_vert(i) )->push_vert( vert_ptr_0 );
	}
	
	for ( i = 0 ; i < vert_ptr_1->get_cell_num() ; ++i ) {
		cell_ptr = static_cast<c_cell*>( vert_ptr_1->get_cell(i) );
		vert_ptr_0->push_cell( static_cast<c_cell*>( vert_ptr_1->get_cell(i) ) );
		static_cast<c_cell*>( vert_ptr_1->get_cell(i) )->remove_vert( vert_ptr_1 );
		static_cast<c_cell*>( vert_ptr_1->get_cell(i) )->push_vert( vert_ptr_0 );
	}
	
	for ( i = 0 ; i < vert_ptr_1->get_face_num() ; ++i ) {
		vert_ptr_0->push_face( static_cast<c_face*>( vert_ptr_1->get_face(i) ) );
		static_cast<c_face*>( vert_ptr_1->get_face(i) )->remove_vert( vert_ptr_1 );
		static_cast<c_face*>( vert_ptr_1->get_face(i) )->push_vert( vert_ptr_0 );
	}
	

  for ( i = 0 ; i < vert_ptr_0->get_cell_num() ; ++i ) {
    cell_ptr = static_cast<c_cell*>( vert_ptr_0->get_cell(i) );
    cell_ptr->update_geometry();
    for ( j = 0 ; j < cell_ptr->get_face_num() ; ++j ) {
    static_cast<c_face*>( cell_ptr->get_face(j) )->update_geometry() ;
    }
  }

	face->clear_all();
	for ( std::list<NIFS::c_face*>::iterator iter = m_faces.begin(); 
		  iter != m_faces.end(); 
		  ++iter ) {
		if ( (*iter) == face ) {
			face_iter_ret = m_faces.erase( face_iter );
			delete face;			
			break;
		}
	}

	vert_ptr_1->clear_all();
	m_verts.remove( vert_ptr_1 );
	delete vert_ptr_1;
	vert_ptr_1 = NULL;
	
	return face_iter_ret;
}


int
NIFS::c_mesh_adap::collapse_face2( c_face* face )
{

	// special treatment for trangles
	unsigned i;
	unsigned j;
	c_cell* cell_ptr; //cell_ptr
	c_vertex* vert_ptr_0 = static_cast<c_vertex*>( face->get_vert(0) );
	c_vertex* vert_ptr_1 = static_cast<c_vertex*>( face->get_vert(1) );
	c_face* face_ptr_0;
	c_face* face_ptr_1;
	c_cell* cell_ptr_0;
	c_cell* cell_ptr_1;

//	c_vector_2d new_pos = face->get_collapse_pos();
	
	if ( ( vert_ptr_0->get_type() == 2 ) && 
		 ( vert_ptr_1->get_type() == 2 ) ) 
		return -1; // both face end-points are corner vertices



	if ( ( vert_ptr_0->get_type() != 0 ) && 
		 ( vert_ptr_1->get_type() != 0 ) &&
		 ( face->get_cell_num() == 2  ) )
		return -2; // face joins 2 boundary points and cannot be removed


	for ( i = 0 ; i < face->get_cell_num(); ++i ) {
		// for triangular mesh
		if ( static_cast<c_cell*>( face->get_cell(i) )->get_face_num() != 3 ) 
			return -3; // at least one of the neighbours is non-triangular

	}

	// vert_0 is retained and vert_1 is removed
	if ( vert_ptr_1->get_type() > vert_ptr_0->get_type() ) {
		c_vertex* vert_temp = vert_ptr_1;
		vert_ptr_1 = vert_ptr_0;
		vert_ptr_0 = vert_temp;
	}

	c_vector_2d new_pos = face->get_collapse_pos();
//	c_vector_2d new_pos = vert_ptr_0->get_point().get_pos();

	for ( i = 0 ; i < vert_ptr_0->get_cell_num() ; ++i ) {
		cell_ptr = static_cast<c_cell*>( vert_ptr_0->get_cell(i) );
		if ( ( cell_ptr != face->get_cell( 0 ) ) &&
			 ( cell_ptr != face->get_cell( 1 ) ) &&
			 ( cell_ptr->is_valid_move( vert_ptr_0, new_pos ) == false )
		   )
			return -4;
	}

	for ( i = 0 ; i < vert_ptr_1->get_cell_num() ; ++i ) {
		cell_ptr = static_cast<c_cell*>( vert_ptr_1->get_cell(i) );
		if ( ( cell_ptr != face->get_cell( 0 ) ) &&
			 ( cell_ptr != face->get_cell( 1 ) ) &&
			 ( cell_ptr->is_valid_move( vert_ptr_1, new_pos ) == false )
		   )
			return -4;
	}


	for ( i = 0 ; i < face->get_cell_num(); ++i ) {
		cell_ptr_0 = NULL;
		cell_ptr_1 = NULL;
		face_ptr_0 = NULL;
		face_ptr_1 = NULL;

		cell_ptr = static_cast<c_cell*>( face->get_cell(i) );
		for ( j = 0 ; j < cell_ptr->get_face_num(); ++j ) {
			if ( face != static_cast<c_face*>( cell_ptr->get_face(j) ) ) {
				if ( static_cast<c_face*>( cell_ptr->get_face(j) )->is_vert_nbr( vert_ptr_0 ) ) {
					face_ptr_0 = static_cast<c_face*>( cell_ptr->get_face(j) );
				} else {
					face_ptr_1 = static_cast<c_face*>( cell_ptr->get_face(j) );
				}
			}
		}
			
		if ( face_ptr_0->get_cell_num() > 1) {
			if ( face_ptr_0->get_cell(0) == cell_ptr ) {
				cell_ptr_0 = static_cast<c_cell*>( face_ptr_0->get_cell(1) );
			} else {
				cell_ptr_0 = static_cast<c_cell*>( face_ptr_0->get_cell(0) );
			}
		}

		if ( face_ptr_1->get_cell_num() > 1) {
			if ( face_ptr_1->get_cell(0) == cell_ptr ) {
				cell_ptr_1 = static_cast<c_cell*>( face_ptr_1->get_cell(1) );
			} else {
				cell_ptr_1 = static_cast<c_cell*>( face_ptr_1->get_cell(0) );
			}
		}

		if ( cell_ptr_0 != NULL) {
			cell_ptr_0->remove_cell( cell_ptr );
			if ( cell_ptr_1 != NULL ) 
				cell_ptr_0->push_cell( cell_ptr_1 );
		}

		if ( cell_ptr_1 != NULL) {
			cell_ptr_1->remove_cell( cell_ptr );
			cell_ptr_1->remove_face( face_ptr_1 );
			cell_ptr_1->push_face( face_ptr_0 );
			if ( cell_ptr_0 != NULL ) 
				cell_ptr_1->push_cell( cell_ptr_0 );
		}
			
		face_ptr_0->remove_cell( cell_ptr );
		if ( cell_ptr_1 != NULL )
			face_ptr_0->push_cell( cell_ptr_1 );
		vert_ptr_0->remove_cell( cell_ptr );
		vert_ptr_1->remove_cell( cell_ptr );
		vert_ptr_1->remove_face( face_ptr_1 );
		if ( vert_ptr_1 == static_cast<c_vertex*>( face_ptr_1->get_vert(0) ) ) {
			vert_ptr_1->remove_vert( static_cast<c_vertex*>( face_ptr_1->get_vert(1) ) );
			static_cast<c_vertex*>( face_ptr_1->get_vert(1) )->remove_vert( vert_ptr_1 );
			static_cast<c_vertex*>( face_ptr_1->get_vert(1) )->remove_face( face_ptr_1 );
			static_cast<c_vertex*>( face_ptr_1->get_vert(1) )->remove_cell( cell_ptr );
		} else {
			vert_ptr_1->remove_vert( static_cast<c_vertex*>( face_ptr_1->get_vert(0) ) );
			static_cast<c_vertex*>( face_ptr_1->get_vert(0) )->remove_vert( vert_ptr_1 );
			static_cast<c_vertex*>( face_ptr_1->get_vert(0) )->remove_face( face_ptr_1 );
			static_cast<c_vertex*>( face_ptr_1->get_vert(0) )->remove_cell( cell_ptr );
		}

		m_faces.remove( face_ptr_1 );
		delete face_ptr_1;
		face_ptr_1 = NULL;

		m_cells.remove( cell_ptr );
		delete cell_ptr;
		cell_ptr = NULL;
	}

		
	vert_ptr_1->remove_vert( vert_ptr_0 );
	vert_ptr_1->remove_face( face );
	vert_ptr_0->remove_vert( vert_ptr_1 );
	vert_ptr_0->remove_face( face );

	vert_ptr_0->set_pos( new_pos );
	
	
	for ( i = 0 ; i < vert_ptr_1->get_vert_num() ; ++i ) {
		vert_ptr_0->push_vert( static_cast<c_vertex*>( vert_ptr_1->get_vert(i) ) );
		static_cast<c_vertex*>( vert_ptr_1->get_vert(i) )->remove_vert( vert_ptr_1 );
		static_cast<c_vertex*>( vert_ptr_1->get_vert(i) )->push_vert( vert_ptr_0 );
	}
	
	for ( i = 0 ; i < vert_ptr_1->get_cell_num() ; ++i ) {
		cell_ptr = static_cast<c_cell*>( vert_ptr_1->get_cell(i) );
		vert_ptr_0->push_cell( static_cast<c_cell*>( vert_ptr_1->get_cell(i) ) );
		static_cast<c_cell*>( vert_ptr_1->get_cell(i) )->remove_vert( vert_ptr_1 );
		static_cast<c_cell*>( vert_ptr_1->get_cell(i) )->push_vert( vert_ptr_0 );
	}
	
	for ( i = 0 ; i < vert_ptr_1->get_face_num() ; ++i ) {
		vert_ptr_0->push_face( static_cast<c_face*>( vert_ptr_1->get_face(i) ) );
		static_cast<c_face*>( vert_ptr_1->get_face(i) )->remove_vert( vert_ptr_1 );
		static_cast<c_face*>( vert_ptr_1->get_face(i) )->push_vert( vert_ptr_0 );
	}
	

  for ( i = 0 ; i < vert_ptr_0->get_cell_num() ; ++i ) {
    cell_ptr = static_cast<c_cell*>( vert_ptr_0->get_cell(i) );
    cell_ptr->update_geometry();
    for ( j = 0 ; j < cell_ptr->get_face_num() ; ++j ) {
    static_cast<c_face*>( cell_ptr->get_face(j) )->update_geometry() ;
    }
  }

	face->clear_all();
  m_faces.remove( face );
  delete face;
  face = NULL;  
  
	vert_ptr_1->clear_all();
	m_verts.remove( vert_ptr_1 );
	delete vert_ptr_1;
	vert_ptr_1 = NULL;
	
	return 0;
}


// another version of the same function
// less efficient but more flexible
int
NIFS::c_mesh_adap::collapse_face( c_face* face )
{

	// special treatment for trangles
	unsigned i;
	unsigned j;
	c_cell* cell_ptr = NULL; //cell_ptr
	c_vertex* vert_ptr_0 = static_cast<c_vertex*>( face->get_vert(0) );
	c_vertex* vert_ptr_1 = static_cast<c_vertex*>( face->get_vert(1) );
	c_vertex* vert_ptr = NULL;
	c_face* face_ptr_0 = NULL;
	c_face* face_ptr_1 = NULL;
	c_cell* cell_ptr_0 = NULL;
	c_cell* cell_ptr_1 = NULL;

//	c_vector_2d new_pos = face->get_collapse_pos();
	if ( ( vert_ptr_0->get_type() == 2 ) && 
		 ( vert_ptr_1->get_type() == 2 ) ) 
		return -1; // both face end-points are corner vertices

	if ( ( vert_ptr_0->get_type() != 0 ) && 
		   ( vert_ptr_1->get_type() != 0 ) &&
		   ( face->get_cell_num() == 2  ) )
		return -2; // face joins 2 boundary points and cannot be removed

	for ( i = 0 ; i < face->get_cell_num(); ++i ) {
		// for triangular mesh
		if ( static_cast<c_cell*>( face->get_cell(i) )->get_face_num() != 3 ) 
			return -3; // at least one of the neighbours is non-triangular
	}

	// vert_0 is retained and vert_1 is removed
	if ( vert_ptr_1->get_type() > vert_ptr_0->get_type() ) {
		vert_ptr   = vert_ptr_1;
		vert_ptr_1 = vert_ptr_0;
		vert_ptr_0 = vert_ptr;
	}


  c_vector_2d new_pos = vert_ptr_0->get_point().get_pos();

	for ( i = 0 ; i < vert_ptr_1->get_cell_num() ; ++i ) {
		cell_ptr = static_cast<c_cell*>( vert_ptr_1->get_cell(i) );
		if ( ( cell_ptr != face->get_cell( 0 ) ) &&
         ( cell_ptr != face->get_cell( 1 ) ) &&
         ( cell_ptr->is_valid_move( vert_ptr_1, new_pos ) == false ) )
			return -4;
	}

	for ( i = 0 ; i < face->get_cell_num() ; ++i ) {
		cell_ptr_0 = NULL;
		cell_ptr_1 = NULL;
		face_ptr_0 = NULL;
		face_ptr_1 = NULL;

		cell_ptr = static_cast<c_cell*>( face->get_cell(i) );
		for ( j = 0 ; j < cell_ptr->get_face_num(); ++j ) {
			if ( face != static_cast<c_face*>( cell_ptr->get_face(j) ) ) {
				if ( static_cast<c_face*>( cell_ptr->get_face(j) )->is_vert_nbr( vert_ptr_0 ) ) {
					face_ptr_0 = static_cast<c_face*>( cell_ptr->get_face(j) );
				} else {
					face_ptr_1 = static_cast<c_face*>( cell_ptr->get_face(j) );
				}
			}
		}

		for ( j = 0 ; j < cell_ptr->get_vert_num(); ++j ) {
		  if ( face->is_vert_nbr( cell_ptr->get_vert( j ) ) == false )
		    vert_ptr = static_cast<c_vertex*>( cell_ptr->get_vert( j ) );
		}
			
		if ( face_ptr_0->get_cell_num() == 2) {
			if ( face_ptr_0->get_cell(0) == cell_ptr ) {
				cell_ptr_0 = static_cast<c_cell*>( face_ptr_0->get_cell(1) );
			} else {
				cell_ptr_0 = static_cast<c_cell*>( face_ptr_0->get_cell(0) );
			}
		}

		if ( face_ptr_1->get_cell_num() == 2) {
			if ( face_ptr_1->get_cell(0) == cell_ptr ) {
				cell_ptr_1 = static_cast<c_cell*>( face_ptr_1->get_cell(1) );
			} else {
				cell_ptr_1 = static_cast<c_cell*>( face_ptr_1->get_cell(0) );
			}
		}

		if ( cell_ptr_0 != NULL) {
			cell_ptr_0->remove_cell( cell_ptr );
			if ( cell_ptr_1 != NULL ) 
				cell_ptr_0->push_cell( cell_ptr_1 );
		}

		if ( cell_ptr_1 != NULL) {
			cell_ptr_1->remove_cell( cell_ptr );
			cell_ptr_1->remove_face( face_ptr_1 );
			cell_ptr_1->push_face( face_ptr_0 );
			cell_ptr_1->remove_vert( vert_ptr_1 );
			cell_ptr_1->push_vert( vert_ptr_0 );
			if ( cell_ptr_0 != NULL ) 
				cell_ptr_1->push_cell( cell_ptr_0 );
		}
			
		face_ptr_0->remove_cell( cell_ptr );
		if ( cell_ptr_1 != NULL ) {
			face_ptr_0->push_cell( cell_ptr_1 );
    }

    vert_ptr_0->remove_cell( cell_ptr );
    vert_ptr_1->remove_cell( cell_ptr );
    vert_ptr_1->remove_face( face_ptr_1 );
    vert_ptr_1->remove_vert( vert_ptr );

    vert_ptr->remove_cell( cell_ptr );
    vert_ptr->remove_face( face_ptr_1 );
    vert_ptr->remove_vert( vert_ptr_1 );
    
		m_faces.remove( face_ptr_1 );
		delete face_ptr_1;
		face_ptr_1 = NULL;

		m_cells.remove( cell_ptr );
		delete cell_ptr;
		cell_ptr = NULL;
	}
		
	vert_ptr_1->remove_vert( vert_ptr_0 );
	vert_ptr_1->remove_face( face );
	vert_ptr_0->remove_vert( vert_ptr_1 );
	vert_ptr_0->remove_face( face );

//	vert_ptr_0->set_pos( new_pos );
	
	for ( i = 0 ; i < vert_ptr_1->get_vert_num() ; ++i ) {
		vert_ptr_0->push_vert( static_cast<c_vertex*>( vert_ptr_1->get_vert(i) ) );
		static_cast<c_vertex*>( vert_ptr_1->get_vert(i) )->remove_vert( vert_ptr_1 );
		static_cast<c_vertex*>( vert_ptr_1->get_vert(i) )->push_vert( vert_ptr_0 );
	}
	
	for ( i = 0 ; i < vert_ptr_1->get_cell_num() ; ++i ) {
//		cell_ptr = static_cast<c_cell*>( vert_ptr_1->get_cell(i) );
		vert_ptr_0->push_cell( static_cast<c_cell*>( vert_ptr_1->get_cell(i) ) );
		static_cast<c_cell*>( vert_ptr_1->get_cell(i) )->remove_vert( vert_ptr_1 );
		static_cast<c_cell*>( vert_ptr_1->get_cell(i) )->push_vert( vert_ptr_0 );
	}
	
	for ( i = 0 ; i < vert_ptr_1->get_face_num() ; ++i ) {
		vert_ptr_0->push_face( dynamic_cast<c_face*>( vert_ptr_1->get_face(i) ) );
		static_cast<c_face*>( vert_ptr_1->get_face(i) )->remove_vert( vert_ptr_1 );
		static_cast<c_face*>( vert_ptr_1->get_face(i) )->push_vert( vert_ptr_0 );
	}


  for ( i = 0 ; i < vert_ptr_0->get_cell_num() ; ++i ) {
    cell_ptr = static_cast<c_cell*>( vert_ptr_0->get_cell(i) );
    cell_ptr->update_geometry();
    for ( j = 0 ; j < cell_ptr->get_face_num() ; ++j ) {
      static_cast<c_face*>( cell_ptr->get_face(j) )->update_geometry() ;
    }
  }

	face->clear_all();
	m_faces.remove( face );
	delete face;				
	face = NULL;

	vert_ptr_1->clear_all();
	m_verts.remove( vert_ptr_1 );
	delete vert_ptr_1;
	vert_ptr_1 = NULL;

	return 0;
}




//------------------------------------------------------------------------------
// removes a vertex
// note: vert_0 is retained and vert_1 is removed
//
void 
NIFS::c_mesh_adap::remove_vert( c_vertex* vert )
{
	// special treatment for trangles
	unsigned i;
	unsigned j;
	c_face* face = NULL;
	c_cell* cell_ptr; //cell_ptr
	c_vertex* vert_ptr_0;
	c_vertex* vert_ptr_1 = vert;
	c_face* face_ptr_0;
	c_face* face_ptr_1;
	c_cell* cell_ptr_0;
	c_cell* cell_ptr_1;

	
  if ( vert->get_type() == 2 ) 
    return; // corner vertex
  else if ( vert->get_type() == 1 ) { // boundary but not corner
    for ( i = 0 ; i < vert->get_face_num() ; ++i ) {
      if ( static_cast<c_face*>( vert->get_face(i) )->get_type() != 0 ) // boundary face
        face = static_cast<c_face*>( vert->get_face(i) ); 
    }
  } else 
    face = static_cast<c_face*>( vert->get_face( 0 ) );


  if ( face->get_vert( 0 ) == vert )
    vert_ptr_0 = static_cast<c_vertex*>( face->get_vert( 1 ) );
  else
    vert_ptr_0 = static_cast<c_vertex*>( face->get_vert( 0 ) );

  for ( i = 0 ; i < vert->get_cell_num(); ++i ) {
    if ( static_cast<c_cell*>( vert->get_cell(i) )->get_face_num() != 3 ) 
      return;	// at least one of the neighbours is not triangular
  }

  c_vector_2d new_pos = vert_ptr_0->get_point().get_pos();

  for ( i = 0 ; i < vert_ptr_1->get_cell_num() ; ++i ) {
    cell_ptr = static_cast<c_cell*>( vert_ptr_1->get_cell( i ) );
      if ( ( cell_ptr != face->get_cell( 0 ) ) &&
           ( cell_ptr != face->get_cell( 1 ) ) &&
           ( cell_ptr->is_valid_move( vert_ptr_1, new_pos ) == false )
         )
        return;  // inhibits folded mesh
  }

  for ( i = 0 ; i < face->get_cell_num() ; ++i ) {
    cell_ptr_0 = NULL;
    cell_ptr_1 = NULL;
    face_ptr_0 = NULL;
    face_ptr_1 = NULL;

    cell_ptr = static_cast<c_cell*>( face->get_cell(i) );
		for ( j = 0 ; j < cell_ptr->get_face_num(); ++j ) {
			if ( face != static_cast<c_face*>( cell_ptr->get_face(j) ) ) {
				if ( static_cast<c_face*>( cell_ptr->get_face(j) )->is_vert_nbr( vert_ptr_0 ) ) {
					face_ptr_0 = static_cast<c_face*>( cell_ptr->get_face(j) );
				} else {
					face_ptr_1 = static_cast<c_face*>( cell_ptr->get_face(j) );
				}
			}
		}
			
		if ( face_ptr_0->get_cell_num() > 1) {
			if ( face_ptr_0->get_cell(0) == cell_ptr ) {
				cell_ptr_0 = static_cast<c_cell*>( face_ptr_0->get_cell(1) );
			} else {
				cell_ptr_0 = static_cast<c_cell*>( face_ptr_0->get_cell(0) );
			}
		}

		if ( face_ptr_1->get_cell_num() > 1) {
			if ( face_ptr_1->get_cell(0) == cell_ptr ) {
				cell_ptr_1 = static_cast<c_cell*>( face_ptr_1->get_cell(1) );
			} else {
				cell_ptr_1 = static_cast<c_cell*>( face_ptr_1->get_cell(0) );
			}
		}

		if ( cell_ptr_0 != NULL) {
			cell_ptr_0->remove_cell( cell_ptr );
			if ( cell_ptr_1 != NULL ) 
				cell_ptr_0->push_cell( cell_ptr_1 );
		}

		if ( cell_ptr_1 != NULL) {
			cell_ptr_1->remove_cell( cell_ptr );
			cell_ptr_1->remove_face( face_ptr_1 );
			cell_ptr_1->push_face( face_ptr_0 );
			if ( cell_ptr_0 != NULL ) 
				cell_ptr_1->push_cell( cell_ptr_0 );
		}
			
		face_ptr_0->remove_cell( cell_ptr );
		if ( cell_ptr_1 != NULL )
			face_ptr_0->push_cell( cell_ptr_1 );
		vert_ptr_0->remove_cell( cell_ptr );
		vert_ptr_1->remove_cell( cell_ptr );
		vert_ptr_1->remove_face( face_ptr_1 );
		if ( vert_ptr_1 == static_cast<c_vertex*>( face_ptr_1->get_vert(0) ) ) {
			vert_ptr_1->remove_vert( static_cast<c_vertex*>( face_ptr_1->get_vert(1) ) );
			static_cast<c_vertex*>( face_ptr_1->get_vert(1) )->remove_vert( vert_ptr_1 );
			static_cast<c_vertex*>( face_ptr_1->get_vert(1) )->remove_face( face_ptr_1 );
			static_cast<c_vertex*>( face_ptr_1->get_vert(1) )->remove_cell( cell_ptr );
		} else {
			vert_ptr_1->remove_vert( static_cast<c_vertex*>( face_ptr_1->get_vert(0) ) );
			static_cast<c_vertex*>( face_ptr_1->get_vert(0) )->remove_vert( vert_ptr_1 );
			static_cast<c_vertex*>( face_ptr_1->get_vert(0) )->remove_face( face_ptr_1 );
			static_cast<c_vertex*>( face_ptr_1->get_vert(0) )->remove_cell( cell_ptr );
		}

		m_faces.remove( face_ptr_1 );
		delete face_ptr_1;
		face_ptr_1 = NULL;

		m_cells.remove( cell_ptr );
		delete cell_ptr;
		cell_ptr = NULL;
	}

		
	vert_ptr_1->remove_vert( vert_ptr_0 );
	vert_ptr_1->remove_face( face );
	vert_ptr_0->remove_vert( vert_ptr_1 );
	vert_ptr_0->remove_face( face );

	
	for ( i = 0 ; i < vert_ptr_1->get_vert_num() ; ++i ) {
		vert_ptr_0->push_vert( static_cast<c_vertex*>( vert_ptr_1->get_vert(i) ) );
		static_cast<c_vertex*>( vert_ptr_1->get_vert(i) )->remove_vert( vert_ptr_1 );
		static_cast<c_vertex*>( vert_ptr_1->get_vert(i) )->push_vert( vert_ptr_0 );
	}
	
	for ( i = 0 ; i < vert_ptr_1->get_cell_num() ; ++i ) {
		cell_ptr = static_cast<c_cell*>( vert_ptr_1->get_cell(i) );
		vert_ptr_0->push_cell( static_cast<c_cell*>( vert_ptr_1->get_cell(i) ) );
		static_cast<c_cell*>( vert_ptr_1->get_cell(i) )->remove_vert( vert_ptr_1 );
		static_cast<c_cell*>( vert_ptr_1->get_cell(i) )->push_vert( vert_ptr_0 );
	}
	
	for ( i = 0 ; i < vert_ptr_1->get_face_num() ; ++i ) {
		vert_ptr_0->push_face( static_cast<c_face*>( vert_ptr_1->get_face(i) ) );
		static_cast<c_face*>( vert_ptr_1->get_face(i) )->remove_vert( vert_ptr_1 );
		static_cast<c_face*>( vert_ptr_1->get_face(i) )->push_vert( vert_ptr_0 );
	}
	

  for ( i = 0 ; i < vert_ptr_0->get_cell_num() ; ++i ) {
    cell_ptr = static_cast<c_cell*>( vert_ptr_0->get_cell(i) );
    cell_ptr->update_geometry();
    for ( j = 0 ; j < cell_ptr->get_face_num() ; ++j ) {
    static_cast<c_face*>( cell_ptr->get_face(j) )->update_geometry() ;
    }
  }
	
  face->clear_all();
  m_faces.remove( face );
  delete face;			
  face = NULL;  
  
  vert_ptr_1->clear_all();
  m_verts.remove( vert_ptr_1 );
  delete vert_ptr_1;
  vert_ptr_1 = NULL;
	
  return;
}


// 	swaps a face
//	return valuse:
//		-1:	boundary face
//		-2: one or both neighbouring cells are not triangle
//		-3:	the associated quadrilateral is concave
//		 n: number of swaped faces ( n >= 0 )
int 
NIFS::c_mesh_adap::delaunay2( c_face* face )
{
	// ignore boundary faces
	if ( face->get_cell_num() != 2 ) // non-internal face
		return -1;

	// only applicable to triangular cells
	if ( ( static_cast<c_cell*>( face->get_cell( 0 ) )->get_face_num() != 3 ) ||
       ( static_cast<c_cell*>( face->get_cell( 1 ) )->get_face_num() != 3 ) ) 
		return -2;

	// create temporary cells
	int cell_num = 2;
	std::vector< c_cell* > cell_temp( cell_num );
	cell_temp[0] = static_cast<c_cell*>( face->get_cell( 0 ) );
	cell_temp[1] = static_cast<c_cell*>( face->get_cell( 1 ) );

	// cavity boundary vertices
	int vert_num = 4;
	std::vector< c_vertex* > bnd_vert_temp( vert_num );	
	bnd_vert_temp[0] = static_cast<c_vertex*>( face->get_vert( 0 ) );
	bnd_vert_temp[1] = static_cast<c_vertex*>( face->get_vert( 1 ) );
	for ( int i = 0 ; i < 3 ; ++i) {
		if ( ( cell_temp[0]->get_vert(i) != bnd_vert_temp[0] ) &&
			 ( cell_temp[0]->get_vert(i) != bnd_vert_temp[1] ) )
			bnd_vert_temp[2] = static_cast<c_vertex*>( cell_temp[0]->get_vert(i) );
		if ( ( cell_temp[1]->get_vert(i) != bnd_vert_temp[0] ) &&
			 ( cell_temp[1]->get_vert(i) != bnd_vert_temp[1] ) )
			bnd_vert_temp[3] = static_cast<c_vertex*>( cell_temp[1]->get_vert(i) );
	}

	// ignore concave quadrilaterals
	c_vector_2d vec_0 = bnd_vert_temp[1]->get_point().get_pos() - bnd_vert_temp[0]->get_point().get_pos();
	c_vector_2d vec_1 = bnd_vert_temp[2]->get_point().get_pos() - bnd_vert_temp[0]->get_point().get_pos();
	c_vector_2d vec_2 = bnd_vert_temp[3]->get_point().get_pos() - bnd_vert_temp[0]->get_point().get_pos();
	double theta_1 = vec_0.get_angle( vec_1 ) + vec_0.get_angle( vec_2 );
	if ( theta_1 > pi() )
		return -3;

	vec_0 = bnd_vert_temp[0]->get_point().get_pos() - bnd_vert_temp[1]->get_point().get_pos();
	vec_1 = bnd_vert_temp[2]->get_point().get_pos() - bnd_vert_temp[1]->get_point().get_pos();
	vec_2 = bnd_vert_temp[3]->get_point().get_pos() - bnd_vert_temp[1]->get_point().get_pos();
	double theta_2 = vec_0.get_angle( vec_1 ) + vec_0.get_angle( vec_2 );
	if ( theta_2 > pi() )
		return -3;

	if ( (theta_1 + theta_2) > pi() )
		return 0;
  
  swap_face( face );		

	return 1;
}

// anisotropic trinagulation
int 
NIFS::c_mesh_adap::delaunay( c_face* face )
{
	// ignore boundary faces
	if ( face->get_cell_num() != 2 ) // non-internal face
		return -1;

	// only applicable to triangular cells
	if ( ( static_cast<c_cell*>( face->get_cell( 0 ) )->get_face_num() != 3 ) ||
       ( static_cast<c_cell*>( face->get_cell( 1 ) )->get_face_num() != 3 ) ) 
		return -2;

	// create temporary cells
	int cell_num = 2;
	std::vector< c_cell* > cell_temp( cell_num );
	cell_temp[0] = static_cast<c_cell*>( face->get_cell( 0 ) );
	cell_temp[1] = static_cast<c_cell*>( face->get_cell( 1 ) );

	// cavity boundary vertices
	// vert_num = 4;
	std::vector< c_vertex* > bnd_vert_temp( 4 );	
	std::vector< c_matrix_2d > metric( 4 );	
	
	bnd_vert_temp[0] = static_cast<c_vertex*>( face->get_vert( 0 ) );
	metric[0] = m_bg_mesh->get_size_matrix( bnd_vert_temp[0]->get_point().get_pos() );
	bnd_vert_temp[1] = static_cast<c_vertex*>( face->get_vert( 1 ) );
	metric[1] = m_bg_mesh->get_size_matrix( bnd_vert_temp[1]->get_point().get_pos() );
	for ( int i = 0 ; i < 3 ; ++i) {
		if ( ( cell_temp[0]->get_vert(i) != bnd_vert_temp[0] ) &&
			 ( cell_temp[0]->get_vert(i) != bnd_vert_temp[1] ) ) {
			bnd_vert_temp[2] = static_cast<c_vertex*>( cell_temp[0]->get_vert(i) );
			metric[2] = m_bg_mesh->get_size_matrix( bnd_vert_temp[2]->get_point().get_pos() );
		}
		if ( ( cell_temp[1]->get_vert(i) != bnd_vert_temp[0] ) &&
			 ( cell_temp[1]->get_vert(i) != bnd_vert_temp[1] ) ) {
			bnd_vert_temp[3] = static_cast<c_vertex*>( cell_temp[1]->get_vert(i) );
			metric[3] = m_bg_mesh->get_size_matrix( bnd_vert_temp[3]->get_point().get_pos() );
		}
	}

	// ignore concave quadrilaterals
	std::vector< c_vector_2d > vect( 5 );
	vect[0] = bnd_vert_temp[1]->get_point().get_pos() - bnd_vert_temp[0]->get_point().get_pos();
	vect[1] = bnd_vert_temp[2]->get_point().get_pos() - bnd_vert_temp[0]->get_point().get_pos();
	vect[2] = bnd_vert_temp[3]->get_point().get_pos() - bnd_vert_temp[0]->get_point().get_pos();
	vect[3] = bnd_vert_temp[2]->get_point().get_pos() - bnd_vert_temp[1]->get_point().get_pos();
	vect[4] = bnd_vert_temp[3]->get_point().get_pos() - bnd_vert_temp[1]->get_point().get_pos();
	
	std::vector< double > ang( 6 );
	c_vector_2d v1,v2;
	v1 = metric[0] * vect[0];
	v2 = metric[0] * vect[1];
	ang[0] =  v1.get_angle( v2 );
	v1 = metric[0] * vect[0];
	v2 = metric[0] * vect[2];
	ang[1] =  v1.get_angle( v2 );
	v1 = metric[1] * (-vect[0]);
	v2 = metric[1] * vect[3];
	ang[2] =  v1.get_angle( v2 );
	v1 = metric[1] * (-vect[0]);
	v2 = metric[1] * vect[4];
	ang[3] = v1.get_angle( v2 );
	v1 = metric[2] * (-vect[1]);
	v2 = metric[2] * (-vect[3]);
	ang[4] = v1.get_angle( v2 );
	v1 = metric[3] * (-vect[2]);
	v2 = metric[3] * (-vect[4]);
	ang[5] = v1.get_angle( v2 );
	
	
	if ( ( ang[0] + ang[1] > pi() ) || ( ang[2] + ang[3] > pi() ) )
		return -3;


	if ( ang[0] + ang[1] + ang[2] + ang[3] > ang[4] + ang[5]  )
		return 0;
  
  swap_face( face );		

	return 1;
}


int 
NIFS::c_mesh_adap::swap_face( c_face* face )
{
	int i;

	// ignore boundary faces
	if ( face->get_cell_num() != 2 ) // non-internal face
		return -1;

	// only applicable to triangular cells
	if ( ( static_cast<c_cell*>( face->get_cell( 0 ) )->get_face_num() != 3 ) ||
       ( static_cast<c_cell*>( face->get_cell( 1 ) )->get_face_num() != 3 ) ) 
		return -2;

	// create temporary cells
	int cell_num = 2;
	std::vector< c_cell* > cell_temp( cell_num );
	cell_temp[0] = static_cast<c_cell*>( face->get_cell( 0 ) );
	cell_temp[1] = static_cast<c_cell*>( face->get_cell( 1 ) );


	// cavity boundary vertices
	int vert_num = 4;
	std::vector< c_vertex* > bnd_vert_temp( vert_num );	
	bnd_vert_temp[0] = static_cast<c_vertex*>( face->get_vert( 0 ) );
	bnd_vert_temp[1] = static_cast<c_vertex*>( face->get_vert( 1 ) );
	for ( i = 0 ; i < 3 ; ++i) {
		if ( ( cell_temp[0]->get_vert(i) != bnd_vert_temp[0] ) &&
			 ( cell_temp[0]->get_vert(i) != bnd_vert_temp[1] ) )
			bnd_vert_temp[2] = static_cast<c_vertex*>( cell_temp[0]->get_vert(i) );
		if ( ( cell_temp[1]->get_vert(i) != bnd_vert_temp[0] ) &&
			 ( cell_temp[1]->get_vert(i) != bnd_vert_temp[1] ) )
			bnd_vert_temp[3] = static_cast<c_vertex*>( cell_temp[1]->get_vert(i) );
	}
	
	c_face* face_ptr;
	// cavity boundary faces
	int face_num = 4;
	std::vector< c_face* > bnd_face_temp( face_num );	
	for ( i = 0 ; i < 3 ; ++i ) {
		if  ( cell_temp[0]->get_face(i) != face ) {
			face_ptr = static_cast<c_face*>( cell_temp[0]->get_face(i) );
			if ( face_ptr->is_vert_nbr( bnd_vert_temp[0] ) )
				bnd_face_temp[0] =  face_ptr;
			else
				bnd_face_temp[1] =  face_ptr;
		}

		if  ( cell_temp[1]->get_face(i) != face ) {
			face_ptr = static_cast<c_face*>( cell_temp[1]->get_face(i) );
			if ( face_ptr->is_vert_nbr( bnd_vert_temp[0] ) )
				bnd_face_temp[2] =  face_ptr;
			else
				bnd_face_temp[3] =  face_ptr;
		}
	}

	// ignore concave quadrilaterals
	c_vector_2d vec_0 = bnd_vert_temp[1]->get_point().get_pos() - bnd_vert_temp[0]->get_point().get_pos();
	c_vector_2d vec_1 = bnd_vert_temp[2]->get_point().get_pos() - bnd_vert_temp[0]->get_point().get_pos();
	c_vector_2d vec_2 = bnd_vert_temp[3]->get_point().get_pos() - bnd_vert_temp[0]->get_point().get_pos();
	double theta_1 = vec_0.get_angle( vec_1 ) + vec_0.get_angle( vec_2 );
	if ( theta_1 > pi() )
		return -3;

	vec_0 = bnd_vert_temp[0]->get_point().get_pos() - bnd_vert_temp[1]->get_point().get_pos();
	vec_1 = bnd_vert_temp[2]->get_point().get_pos() - bnd_vert_temp[1]->get_point().get_pos();
	vec_2 = bnd_vert_temp[3]->get_point().get_pos() - bnd_vert_temp[1]->get_point().get_pos();
	double theta_2 = vec_0.get_angle( vec_1 ) + vec_0.get_angle( vec_2 );
	if ( theta_2 > pi() )
		return -3;

// Delaunay condition
//	if ( (theta_1 + theta_2) > pi() )
//		return 0;
		
	face->clear_verts();
	face->push_vert( bnd_vert_temp[2] );
	face->push_vert( bnd_vert_temp[3] );

	bnd_vert_temp[0]->remove_cell( cell_temp[1] );
	bnd_vert_temp[0]->remove_face( face );
	bnd_vert_temp[0]->remove_vert( bnd_vert_temp[1] );

	bnd_vert_temp[1]->remove_cell( cell_temp[0] );
	bnd_vert_temp[1]->remove_face( face );
	bnd_vert_temp[1]->remove_vert( bnd_vert_temp[0] );

	bnd_vert_temp[2]->push_cell( cell_temp[1] );
	bnd_vert_temp[2]->push_face( face );
	bnd_vert_temp[2]->push_vert( bnd_vert_temp[3] );
	
	bnd_vert_temp[3]->push_cell( cell_temp[0] );
	bnd_vert_temp[3]->push_face( face );
	bnd_vert_temp[3]->push_vert( bnd_vert_temp[2] );

	cell_temp[0]->remove_vert( bnd_vert_temp[1] );
	cell_temp[0]->push_vert( bnd_vert_temp[3] );
	
	cell_temp[1]->remove_vert( bnd_vert_temp[0] );
	cell_temp[1]->push_vert( bnd_vert_temp[2] );
	
	cell_temp[0]->remove_face( bnd_face_temp[1] );
	cell_temp[0]->push_face( bnd_face_temp[2] );
	
	cell_temp[1]->remove_face( bnd_face_temp[2] );
	cell_temp[1]->push_face( bnd_face_temp[1] );
	
	if ( bnd_face_temp[2]->get_cell_num() > 1 ) {
		if ( bnd_face_temp[2]->get_cell( 0 ) == cell_temp[1] ){
			static_cast<c_cell*>( bnd_face_temp[2]->get_cell(1) )->remove_cell( cell_temp[1] );
			static_cast<c_cell*>( bnd_face_temp[2]->get_cell(1) )->push_cell( cell_temp[0] );
			cell_temp[1]->remove_cell( bnd_face_temp[2]->get_cell(1) );
			cell_temp[0]->push_cell( bnd_face_temp[2]->get_cell(1) );
		} else {
			static_cast<c_cell*>( bnd_face_temp[2]->get_cell(0) )->remove_cell( cell_temp[1] );
			static_cast<c_cell*>( bnd_face_temp[2]->get_cell(0) )->push_cell( cell_temp[0] );
			cell_temp[1]->remove_cell( bnd_face_temp[2]->get_cell(0) );
			cell_temp[0]->push_cell( bnd_face_temp[2]->get_cell(0) );
		}
	}

	if ( bnd_face_temp[1]->get_cell_num() > 1 ) {
		if ( bnd_face_temp[1]->get_cell( 0 ) == cell_temp[0] ){
			static_cast<c_cell*>( bnd_face_temp[1]->get_cell(1) )->remove_cell( cell_temp[0] );
			static_cast<c_cell*>( bnd_face_temp[1]->get_cell(1) )->push_cell( cell_temp[1] );
			cell_temp[0]->remove_cell( bnd_face_temp[1]->get_cell(1) );
			cell_temp[1]->push_cell( bnd_face_temp[1]->get_cell(1) );
		} else {
			static_cast<c_cell*>( bnd_face_temp[1]->get_cell(0) )->remove_cell( cell_temp[0] );
			static_cast<c_cell*>( bnd_face_temp[1]->get_cell(0) )->push_cell( cell_temp[1] );
			cell_temp[0]->remove_cell( bnd_face_temp[1]->get_cell(0) );
			cell_temp[1]->push_cell( bnd_face_temp[1]->get_cell(0) );
		}
	}

	bnd_face_temp[1]->remove_cell( cell_temp[0] );
	bnd_face_temp[1]->push_cell( cell_temp[1] ); 
	
	bnd_face_temp[2]->remove_cell( cell_temp[1] );
	bnd_face_temp[2]->push_cell( cell_temp[0] ); 

	for (i=0; i<cell_num; i++)
		cell_temp[i]->update_geometry();

	for (i=0; i<face_num; i++)
		bnd_face_temp[i]->update_geometry();

	face->update_geometry();

	return 1;
}



int
NIFS::c_mesh_adap::merge_tri_to_quad( c_face* face )
{
  if ( face == NULL )
    return -1;
  
	if ( face->get_cell_num() != 2 ) // non-internal face
		return -1;

	// only applicable to triangular cells
  if ( ( static_cast<c_cell*>( face->get_cell( 0 ) )->get_face_num() != 3 ) ||
       ( static_cast<c_cell*>( face->get_cell( 1 ) )->get_face_num() != 3 ) ) 
    return -2;

	// create temporary cells
	c_cell* cell_ptr_0;
	c_cell* cell_ptr_1;
	c_cell* cell_ptr;
	c_face* face_ptr;
	c_vertex* vert_ptr_0;
	c_vertex* vert_ptr_1;
	c_vertex* vert_ptr;
	unsigned i;
	
	cell_ptr_0 = static_cast< c_cell* >( face->get_cell( 0 ) );
	cell_ptr_1 = static_cast< c_cell* >( face->get_cell( 1 ) );

  // modify the local cell neighbourhood
  cell_ptr_0->remove_cell( cell_ptr_1 );
  for ( i = 0 ; i < cell_ptr_1->get_cell_num() ; ++i ) {
    cell_ptr = static_cast< c_cell* >( cell_ptr_1->get_cell( i ) );
    if ( cell_ptr != cell_ptr_0 ) {
      cell_ptr_0->push_cell( cell_ptr );
      cell_ptr->remove_cell( cell_ptr_1 );
      cell_ptr->push_cell( cell_ptr_0 );
    }
  }

  // remove the common face and add the faces of the neighbouring control volume
  cell_ptr_0->remove_face( face );
  for ( i = 0 ; i < cell_ptr_1->get_face_num() ; ++i ) {
    face_ptr = static_cast< c_face* >( cell_ptr_1->get_face( i ) );
    if ( face_ptr != face ) {
      cell_ptr_0->push_face( face_ptr );
      face_ptr->remove_cell( cell_ptr_1 );
      face_ptr->push_cell( cell_ptr_0 );
    }
  }

  // add the extra vertex of the neighbouring control volume
  for ( i = 0 ; i < cell_ptr_1->get_vert_num() ; ++i ) {
    vert_ptr = static_cast< c_vertex* >( cell_ptr_1->get_vert( i ) );
    if ( face->is_vert_nbr( vert_ptr ) == false ) {
      cell_ptr_0->push_vert( vert_ptr );
      vert_ptr->remove_cell( cell_ptr_1 );
      vert_ptr->push_cell( cell_ptr_0 );
    }
  }

  vert_ptr_0 = static_cast< c_vertex* >( face->get_vert( 0 ) );
  vert_ptr_1 = static_cast< c_vertex* >( face->get_vert( 1 ) );
  vert_ptr_0->remove_vert( vert_ptr_1 );
  vert_ptr_0->remove_face( face );
  vert_ptr_0->remove_cell( cell_ptr_1 );
  vert_ptr_1->remove_vert( vert_ptr_0 );
  vert_ptr_1->remove_face( face );
  vert_ptr_1->remove_cell( cell_ptr_1 );
  
  face->clear_all();
  m_faces.remove( face );
  delete face;			
  face = NULL;  
  
  cell_ptr_1->clear_all();
  m_cells.remove( cell_ptr_1 );
  delete cell_ptr_1;
  cell_ptr_1 = NULL;

  cell_ptr_0->update_geometry();
  
	for ( i = 0 ; i < cell_ptr_0->get_face_num() ; ++i ) 
		static_cast< c_face* >( cell_ptr_0->get_face( i ) )->update_geometry();
  
  return 0;
}	


// splits a quad into two triangles
// if the quad is not convex, the method fails. This constraint is to be 
// fixed later. 
int 
NIFS::c_mesh_adap::split_quad_to_tri( c_cell* cell )
{
	if ( cell->get_face_num() != 4 ) // non-quadrilateral cell
		return -1;

  std::vector< c_vertex* > verts( cell->get_vert_num() );
  std::vector< c_face* > faces( cell->get_face_num() );
  c_cell* new_cell = new c_cell(); // new triangular cell
  c_face* new_face = new c_face_int(); // new face splitting up the quad cell
  
  m_faces.push_back( new_face );
  m_cells.push_back( new_cell );
  
  c_vertex* vert_ptr; // temporary vertex
  c_face* face_ptr;   // temporary face
  
  faces[0] = static_cast<c_face*>( cell->get_face( 0 ) );
  verts[0] = static_cast<c_vertex*>( faces[0]->get_vert( 0 ) );
  verts[1] = static_cast<c_vertex*>( faces[0]->get_vert( 1 ) );
  unsigned i;
  unsigned j;
  for ( i = 2 ; i < cell->get_vert_num() ; ++i )
    for ( j = 0 ; j < cell->get_vert_num() ; ++j ) {
      vert_ptr = static_cast<c_vertex*>( cell->get_vert( j ) );
      if ( vert_ptr->is_vert_nbr( verts[i-1] ) && ( vert_ptr != verts[i-2] ) )
        verts[i] = vert_ptr;
    }

  for ( i = 1 ; i < cell->get_face_num() ; ++i )
    for ( j = 0 ; j < cell->get_face_num() ; ++j ) {
      face_ptr = static_cast<c_face*>( cell->get_face( j ) );
      if ( face_ptr->is_vert_nbr( verts[i] ) && ( face_ptr != faces[i-1] ) )
        faces[i] = face_ptr;
    }

  verts[0]->push_vert( verts[2] );
  verts[0]->push_face( new_face );
  verts[0]->push_cell( new_cell );
  verts[2]->push_vert( verts[0] );
  verts[2]->push_face( new_face );
  verts[2]->push_cell( new_cell );
  verts[3]->remove_cell( cell );
  verts[3]->push_cell( new_cell );
  
  faces[2]->remove_cell( cell );
  faces[2]->push_cell( new_cell );
  faces[3]->remove_cell( cell );
  faces[3]->push_cell( new_cell );
  new_face->push_vert( verts[0] );
  new_face->push_vert( verts[2] );
  new_face->push_cell( cell );
  new_face->push_cell( new_cell );
  
  cell->remove_vert( verts[3] );
  cell->remove_face( faces[2] );
  cell->remove_face( faces[3] );
  cell->push_face( new_face );
  
  new_cell->push_vert( verts[0] );
  new_cell->push_vert( verts[2] );
  new_cell->push_vert( verts[3] );
  new_cell->push_face( faces[2] );
  new_cell->push_face( faces[3] );
  new_cell->push_face( new_face );
  
  static_cast<c_cell*>( faces[2]->get_cell( 0 ) )->remove_cell( cell );
  static_cast<c_cell*>( faces[2]->get_cell( 0 ) )->push_cell( new_cell );
  static_cast<c_cell*>( faces[3]->get_cell( 0 ) )->remove_cell( cell );
  static_cast<c_cell*>( faces[3]->get_cell( 0 ) )->push_cell( new_cell );
  
  cell->remove_cell( static_cast<c_cell*>( faces[2]->get_cell( 0 ) ) );
  cell->remove_cell( static_cast<c_cell*>( faces[3]->get_cell( 0 ) ) );
  cell->push_cell( new_cell );
  
  new_cell->push_cell( static_cast<c_cell*>( faces[2]->get_cell( 0 ) ) );
  new_cell->push_cell( static_cast<c_cell*>( faces[3]->get_cell( 0 ) ) );
  new_cell->push_cell( cell );
  
  
  cell->update_geometry();
  new_cell->update_geometry();
  new_face->update_geometry(); 
  for ( i = 1 ; i < 4 ; ++i )
    faces[i]->update_geometry();

  delaunay(new_face); // delaunay condition
   
  return 0;
}


int 
NIFS::c_mesh_adap::enhance_topolog( c_face* face )
{
  c_cell* cell_ptr;
  c_vector_2d vec_0;
  c_vector_2d vec_1;
  std::vector< c_vertex* > vert_ptr( 4 );

	if ( face->get_cell_num() != 2 ) // non-internal face
		return -1;

 
  vert_ptr[0] = static_cast<c_vertex*>( face->get_vert( 0 ) );
  vert_ptr[1] = static_cast<c_vertex*>( face->get_vert( 1 ) );
  for ( unsigned i = 0 ; i < face->get_cell_num() ; ++i ) {
    cell_ptr = static_cast<c_cell*>( face->get_cell( 0 ) );
    if ( cell_ptr->get_face_num() != 3 )
      return -2;
    if ( cell_ptr->get_aspect_ratio() < 1.5 )
      return -3;

    for ( unsigned j = 0 ; j < cell_ptr->get_vert_num() ; ++j ) {
      if ( face->is_vert_nbr( cell_ptr->get_vert( j ) ) == false )
        vert_ptr[i+2] = static_cast<c_vertex*>( cell_ptr->get_vert( j ) );
    }
  }

  vec_0 = vert_ptr[2]->get_point().get_pos() - vert_ptr[0]->get_point().get_pos();
  vec_1 = vert_ptr[3]->get_point().get_pos() - vert_ptr[0]->get_point().get_pos();
  if ( fabs( vec_0.get_unit() * vec_1.get_unit() ) > 0.5 )
    return -4;
  vec_0 = vert_ptr[2]->get_point().get_pos() - vert_ptr[1]->get_point().get_pos();
  vec_1 = vert_ptr[3]->get_point().get_pos() - vert_ptr[1]->get_point().get_pos();
  if ( fabs( vec_0.get_unit() * vec_1.get_unit() ) > 0.5 )
    return -4;
  
    if ( ( vert_ptr[0]->get_vert_num() + vert_ptr[1]->get_vert_num() - 2 ) >
         ( vert_ptr[2]->get_vert_num() + vert_ptr[3]->get_vert_num() ) ) 
    return swap_face( face );

  return 0;
}


int 
NIFS::c_mesh_adap::enhance_topolog( c_vertex* vert )
{
  c_face* face_ptr = NULL;
  c_face* face_ptr_temp = NULL;
  c_face* face_ptr_0 = NULL;
  c_face* face_ptr_1 = NULL;
  c_cell* cell_ptr_0 = NULL;
  c_cell* cell_ptr_1 = NULL;
  
  for ( unsigned i = 0 ; i < vert->get_face_num() ; ++i ) {
    face_ptr = static_cast<c_face*>( vert->get_face( i ) );
    if ( face_ptr->get_cell_num() < 2 )
      continue;
    cell_ptr_0 = static_cast<c_cell*>( face_ptr->get_cell( 0 ) );
    cell_ptr_1 = static_cast<c_cell*>( face_ptr->get_cell( 1 ) );
    if ( ( cell_ptr_0->get_face_num() != 3 ) ||
         ( cell_ptr_1->get_face_num() != 3 ) )
      continue;
    if ( ( cell_ptr_0->get_aspect_ratio() < 1.5 ) ||
         ( cell_ptr_1->get_aspect_ratio() < 1.5 ) )
      continue; 
    for ( unsigned j = 0 ; j < 3 ; ++ j ) {
      face_ptr_temp = static_cast<c_face*>( cell_ptr_0->get_face( j ) );
      if ( ( face_ptr_temp != face_ptr ) &&
           ( vert->is_face_nbr( face_ptr_temp ) ) )
         face_ptr_0 = face_ptr_temp;
    }
    for ( unsigned j = 0 ; j < 3 ; ++ j ) {
      face_ptr_temp = static_cast<c_face*>( cell_ptr_1->get_face( j ) );
      if ( ( face_ptr_temp != face_ptr ) &&
           ( vert->is_face_nbr( face_ptr_temp ) ) )
         face_ptr_1 = face_ptr_temp;
    }
    if ( ( fabs( face_ptr->get_utv() * face_ptr_0->get_utv() ) > 0.5 ) &&
         ( fabs( face_ptr->get_utv() * face_ptr_1->get_utv() ) > 0.5 ) ) {
      merge_tri_to_quad( face_ptr_0 );
      merge_tri_to_quad( face_ptr_1 );
      return 1;
    }
  }
  
  return 0;
}


void 
NIFS::c_mesh_adap::delaunay( c_vertex* vert)
{
  c_face*   face_nbr;
  c_vertex* vert_nbr;
    
  for ( unsigned i = 0 ; i < vert->get_face_num() ; ++i ) {    
    face_nbr = static_cast<c_face*>( vert->get_face( i ) );
    if ( face_nbr->get_vert( 0 ) == vert )
      vert_nbr = static_cast<c_vertex*>( face_nbr->get_vert( 1 ) );
    else
      vert_nbr = static_cast<c_vertex*>( face_nbr->get_vert( 0 ) );
    if ( delaunay( face_nbr ) > 0 )
        delaunay( vert_nbr );
  }
  return;	
}


double 
NIFS::c_mesh_adap::smooth_vertex( c_vertex* vert )
{
  c_vertex* vert_ptr;
//  c_face*   face_ptr;
  c_cell*   cell_ptr;
  c_vector_2d vec( 0.0, 0.0 );
  c_vector_2d del_vec( 0.0, 0.0 );
  c_matrix_2d eig;
  double eig_1, eig_2;
//  double r;
  double avg;
  double d;
  double f;
  double cndno;
  int div;
  int div_old;
  unsigned n = 0;
  int cond = 1;

  if ( vert->get_type() == 2 )
    return 0.0; //corner
  
  avg = 0.0;  
  
  eig = m_bg_mesh->get_size_matrix( vert->get_point().get_pos() ).eig_val();
  eig_1 = std::abs( eig(0)(0) );
  eig_2 = std::abs( eig(1)(1) );
  cndno = ( eig_1 < eig_2 ) ? ( eig_2/eig_1) : ( eig_1/eig_2);
  avg = ( eig_1 < eig_2 ) ? ( 1.0/eig_2) : ( 1.0/eig_1);
  
  for ( unsigned k = 0 ; k < vert->get_cell_num() ; ++k ) {
  cell_ptr = static_cast<c_cell*>( vert->get_cell( k ) );
  for ( unsigned i = 0 ; i < cell_ptr->get_vert_num() ; ++i ) {
    vert_ptr = static_cast<c_vertex*>( cell_ptr->get_vert( i ) );

    cond = 1;
    if ( vert_ptr == vert )
      cond = 0;
      
      
    if ( vert->get_type() != 0 ) { // vert is boundary
      if ( vert_ptr->get_type() == 2 )
        cond = 2; 
		}
		
		if ( cond == 0 ) continue;
		
    d = m_bg_mesh->get_length_ratio( vert->get_point().get_pos() ,
		                           vert_ptr->get_point().get_pos() );
		if ( ( cond == 2 ) && ( d > 1.0 ) ) d = 1.0;
    if ( cell_ptr->get_vert_num() == 3 ) {
      if ( cndno > 4.0 )
        f = - 3.0 * std::pow( d , 6.0 ) * ( 1.0 - std::pow( d , 4.0 ) ) * std::exp( - std::pow( d , 10.0 ) );
      else 
        f = - ( 1.0 - std::pow( d , 2.0 ) ) * std::exp( - std::pow( d , 2.0 ) );
    }
    else if ( cell_ptr->get_vert_num() == 4 )
    f = -( 1.0 - d );
    vec += ( vert_ptr->get_point().get_pos() - vert->get_point().get_pos() ).get_unit() * f ;
    ++n;
  }	
  }

  del_vec = vec * avg / static_cast<double>( n );

  div = div_old = 4;
  do {
    div_old = div;
    del_vec /= static_cast<double>( div );
    vec =  del_vec + vert->get_point().get_pos();
    if ( vert->get_type() == 1 )
      vec = vert->move_vertex_bnd( vec );
    for ( unsigned i = 0 ; i < vert->get_cell_num() ; ++i ) {
      cell_ptr = static_cast<c_cell*>( vert->get_cell( i ) );
      if ( cell_ptr->is_valid_move( vert, vec ) == false ) {
        div *= 2; // the move results in a degenerate mesh
        break;
      }
    }
  } while ( div != div_old );

  del_vec = vec - vert->get_point().get_pos();
  vert->set_point( c_point_2d( vec ) ); 

  for ( unsigned i = 0 ; i < vert->get_cell_num() ; ++i ) {
    cell_ptr = static_cast<c_cell*>( vert->get_cell( i ) );
    cell_ptr->update_geometry();
    for ( unsigned j = 0 ; j < cell_ptr->get_face_num() ; ++j ) 
      static_cast<c_face*>( cell_ptr->get_face( j ) )->update_geometry();  
  }

  
  return ( del_vec ).get_magnit();
}


double 
NIFS::c_mesh_adap::smooth()
{
  double max_move;
  double d = 0.0;
  int n = 0;
  do {
    max_move = 0.0;
    for ( m_vert_iter = m_verts.begin(); 
          m_vert_iter != m_verts.end(); 
          ++m_vert_iter ) {
      d = smooth_vertex( *m_vert_iter );
//      delaunay( *m_vert_iter );
      if ( d > max_move ) {
        max_move = d;
      }
    }
    delaunay();
    n++;
    std::cout << n << ' ' << max_move << std::endl;
  } while ( ( n < m_info.get_params()->get_smt_iter() ) &&
            ( max_move > m_info.get_params()->get_smt_trsh() ) );
  return max_move;
}


void 
NIFS::c_mesh_adap::modify()
{
  unsigned n = 0;
  unsigned min_index;
  unsigned max_index;
  double max_val;
  double min_val;
  bool cond = false;
  double r;
  double avg;
  std::vector< c_vertex* > verts( get_vert_num() );
  std::vector< c_vertex* > vert_temp;
  c_face* face_ptr;
  c_vertex* vert_ptr;
  c_vertex* vert_ptr_0;
  c_vertex* vert_ptr_1;
  
  n = 0;
  for ( m_vert_iter = m_verts.begin(); 
        m_vert_iter != m_verts.end(); 
        ++m_vert_iter ) {
    (*m_vert_iter)->set_active( 0 );
    verts[n] = (*m_vert_iter);
    ++n;
  }
  
  std::random_shuffle( verts.begin(), verts.end() );
  
  n = 0;
  for ( int j = 0 ; j < get_vert_num() ; ++j ) {
    vert_ptr = verts[j];
    cond = false;
    for ( unsigned i = 0 ; i < vert_ptr->get_vert_num() ; ++i ) {
      if ( static_cast<c_vertex*>( vert_ptr->get_vert( i ) )->get_active() != 0 )
        cond = true;
    }
    if ( cond == false ) {
      verts[j]->set_active( 1 );
      vert_temp.push_back( verts[j] );
      ++n;
    }
  }

  
  
  for ( unsigned i = 0 ; i < n ; ++i ) {
    min_index = max_index = 0;
    min_val = max_val = 1.0;
    avg = 0.0;
    for ( unsigned j = 0 ; j < vert_temp[i]->get_face_num() ; ++j ) {
      face_ptr = static_cast<c_face*>( vert_temp[i]->get_face( j ) );
      vert_ptr_0 = static_cast<c_vertex*>( face_ptr->get_vert( 0 ) );
      vert_ptr_1 = static_cast<c_vertex*>( face_ptr->get_vert( 1 ) );

      r = m_bg_mesh->get_length_ratio( vert_ptr_0->get_point().get_pos(),
                                 vert_ptr_1->get_point().get_pos() ); 
      avg += r;

      if ( r > max_val ) {
        max_index = j;
        max_val = r;
      } 

      if ( r < min_val ) {
        min_index = j;
        min_val = r;
      } 
      
    }
    avg /= static_cast<double>( vert_temp[i]->get_face_num() );
    
    if (  ( max_val > 1.4142 ) /*&& ( avg > 1.05)*/ ) {
      split_face( static_cast<c_face*>( vert_temp[i]->get_face( max_index ) ) );    
    }    
    else if (  ( min_val < 0.7071 ) /*&& (avg < 0.95)*/ ) {
      remove_vert( vert_temp[i] );
    }
  }
}


void 
NIFS::c_mesh_adap::refine()
{
	double r;
	unsigned m = m_faces.size();
	m_face_iter = m_faces.begin();
	c_face* face_ptr;
  c_vertex* vert_ptr_0;
  c_vertex* vert_ptr_1;
	
	for ( unsigned n = 0 ; n <= m ; n++ ) 
	{
      face_ptr = static_cast<c_face*>( *m_face_iter );
      vert_ptr_0 = static_cast<c_vertex*>( face_ptr->get_vert( 0 ) );
      vert_ptr_1 = static_cast<c_vertex*>( face_ptr->get_vert( 1 ) );

      r = m_bg_mesh->get_length_ratio( vert_ptr_0->get_point().get_pos(),
                                       vert_ptr_1->get_point().get_pos() 
                                     ); 
  		if ( r > 1.4143 ) split_face( face_ptr );
			++m_face_iter;
	}
}


bool
NIFS::c_mesh_adap::delaunay()
{
  bool cond = false;
	unsigned swap_case = 0;
	unsigned swap_case_old;
	do {
	  swap_case_old = swap_case;
		swap_case = 0;
		for ( m_face_iter = m_faces.begin(); m_face_iter != m_faces.end(); m_face_iter++) 
			if ( delaunay( *m_face_iter ) == 1 ) {
				swap_case++;
				cond = true;
		}
//		std::cout << swap_case << std::endl;
	} while ( swap_case != swap_case_old );
	
	return cond;
}


void 
NIFS::c_mesh_adap::switch_quad()
{
	unsigned m = m_faces.size();
	c_cell* cell_ptr_0;
	c_cell* cell_ptr_1;
	c_face* face_ptr;
	c_vertex* vert_ptr;
	c_vertex* vert_ptr_0;
	c_vertex* vert_ptr_1;
	c_vector_2d vec_0;
	c_vector_2d vec_1;
	c_vector_2d vec;
	c_matrix_2d eig;
	c_matrix_2d trans;
	double theta,eig_1,eig_2,cndno;
	double pi=3.14159;
  
	m_face_iter = m_faces.begin();
	for ( unsigned n = 0 ; n < m ; n++ ) 
	{
    face_ptr = static_cast<c_face*>( *m_face_iter );
    ++m_face_iter;
    if ( face_ptr->get_cell_num() != 2 ) 
      continue; // the face is at the boundary
    cell_ptr_0 = static_cast<c_cell*>( face_ptr->get_cell( 0 ) );
    cell_ptr_1 = static_cast<c_cell*>( face_ptr->get_cell( 1 ) );
    if ( ( cell_ptr_0->get_face_num() != 3 ) ||
         ( cell_ptr_1->get_face_num() != 3 ) )
      continue; // one of the neighbouring cells is not triangle

//     std::cout << "test" <<std::endl;
    trans =  (m_bg_mesh->get_size_matrix( face_ptr->get_fmp().get_pos() ));
    eig = trans.eig_val();
    eig_1 = std::abs( eig(0)(0) );
    eig_2 = std::abs( eig(1)(1) );
    cndno = ( eig_1 < eig_2 ) ? (eig_2/eig_1) : ( eig_1/eig_2);
    if ( cndno < 4.0 ) continue;
    
    eig = trans.eig_vec();
    vec = trans * face_ptr->get_utv();
    vec_0 = eig(0).get_unit();
    vec_1 = eig(1).get_unit();
    theta = atan2( fabs(vec*vec_1),fabs(vec*vec_0) );           
    if ( fabs( theta - pi/4.0) > 0.1 )
      continue; // one of the neighbouring cells has low aspect ratio

         
//    std::cout << n << std::endl;
    merge_tri_to_quad( face_ptr );
	}
	return;
}


void 
NIFS::c_mesh_adap::cleanup_tri( void )
{
	c_cell* cell_ptr_0;
	c_cell* cell_ptr_1;
	c_face* face_ptr;
	c_vertex* vert_ptr = NULL;
	c_vertex* vert_ptr_0 = NULL;
	c_vertex* vert_ptr_1 = NULL;
	c_vector_2d vec_0;
	c_vector_2d vec_1;
	bool cond = true;

	m_face_iter = m_faces.begin();
	for ( unsigned n = 0 ; n < m_faces.size() ; n++ ) 
	{
    if ( cond == false )
      ++m_face_iter;
    face_ptr = dynamic_cast<c_face*>( *m_face_iter );
    cond = false;

    if ( face_ptr->get_cell_num() == 2 ) {
      cell_ptr_0 = static_cast<c_cell*>( face_ptr->get_cell( 0 ) );
      cell_ptr_1 = static_cast<c_cell*>( face_ptr->get_cell( 1 ) );
      if ( ( cell_ptr_0->get_face_num() != 3 ) ||
           ( cell_ptr_1->get_face_num() != 3 ) )
        continue; // one of the neighbouring cells is not triangle

      if ( ( cell_ptr_0->get_aspect_ratio() < 1.5 ) ||
           ( cell_ptr_1->get_aspect_ratio() < 1.5 ) )
        continue; // one of the neighbouring cells has low aspect ratio


      for ( unsigned i = 0 ; i < 3 ; ++i ) {
          if ( face_ptr->is_vert_nbr( cell_ptr_0->get_vert( i ) ) == false )
            vert_ptr_0 = static_cast<c_vertex*>( cell_ptr_0->get_vert( i ) );
          if ( face_ptr->is_vert_nbr( cell_ptr_1->get_vert( i ) ) == false )
            vert_ptr_1 = static_cast<c_vertex*>( cell_ptr_1->get_vert( i ) );
      }
    
      vert_ptr = static_cast<c_vertex*>( face_ptr->get_vert( 0 ) );
      vec_0 = vert_ptr_0->get_point().get_pos() - vert_ptr->get_point().get_pos();
      vec_1 = vert_ptr_1->get_point().get_pos() - vert_ptr->get_point().get_pos();
      if ( vec_0.get_unit() * vec_1.get_unit() > -0.5 )
        continue;

      vert_ptr = static_cast<c_vertex*>( face_ptr->get_vert( 1 ) );
      vec_0 = vert_ptr_0->get_point().get_pos() - vert_ptr->get_point().get_pos();
      vec_1 = vert_ptr_1->get_point().get_pos() - vert_ptr->get_point().get_pos();
      if (  vec_0.get_unit() * vec_1.get_unit() > -0.5 )
        continue;

      m_face_iter = collapse_face( m_face_iter );
      cond = true;
    } 
    else {
      cell_ptr_0 = static_cast<c_cell*>( face_ptr->get_cell( 0 ) );
      if ( cell_ptr_0->get_face_num() != 3 )
        continue; // one of the neighbouring cells is not triangle
      if ( cell_ptr_0->get_aspect_ratio() < 1.5 )
        continue; // one of the neighbouring cells has low aspect ratio

      for ( unsigned i = 0 ; i < 3 ; ++i ) {
          if ( face_ptr->is_vert_nbr( cell_ptr_0->get_vert( i ) ) == false )
            vert_ptr = static_cast<c_vertex*>( cell_ptr_0->get_vert( i ) );
      }
      vert_ptr_0 = static_cast<c_vertex*>( face_ptr->get_vert( 0 ) );
      vert_ptr_1 = static_cast<c_vertex*>( face_ptr->get_vert( 1 ) );
      vec_0 = vert_ptr_0->get_point().get_pos() - vert_ptr->get_point().get_pos();
      vec_1 = vert_ptr_1->get_point().get_pos() - vert_ptr->get_point().get_pos();
      if ( vec_0.get_unit() * vec_1.get_unit() < 0.5 )
        continue;

      m_face_iter = collapse_face( m_face_iter );
      cond = true;
    }
  }
	return;
}

void 
NIFS::c_mesh_adap::pushout_tri( void )
{
//	unsigned m = m_faces.size();
	c_cell* cell_ptr_0;
	c_cell* cell_ptr_1;
	c_cell* cell_ptr;
	c_face* face_ptr;
	c_vertex* vert_ptr = NULL;
	c_vertex* vert_ptr_0 = NULL;
	c_vertex* vert_ptr_1 = NULL;
	c_vector_2d vec_0;
	c_vector_2d vec_1;
	bool cond = true;

	m_face_iter = m_faces.begin();
	for ( unsigned n = 0 ; n < m_faces.size() ; n++ ) 
	{
    if ( cond == false )
      ++m_face_iter;
    face_ptr = dynamic_cast<c_face*>( *m_face_iter );
    cond = false;
    
    if ( face_ptr->get_cell_num() == 1 ) 
      continue;
    else {
      cell_ptr_0 = static_cast<c_cell*>( face_ptr->get_cell( 0 ) );
      cell_ptr_1 = static_cast<c_cell*>( face_ptr->get_cell( 1 ) );
      if ( ( cell_ptr_0->get_face_num() +
             cell_ptr_1->get_face_num() ) != 7 ) 
        continue; // 3+4=7 tri and quad

      if ( ( cell_ptr_0->get_aspect_ratio() < 1.5 ) ||
           ( cell_ptr_1->get_aspect_ratio() < 1.5 ) )
        continue; // one of the neighbouring cells has low aspect ratio

      if ( cell_ptr_0->get_face_num() != 3 ) {
        cell_ptr   = cell_ptr_0;
        cell_ptr_0 = cell_ptr_1;
        cell_ptr_1 = cell_ptr;
      }
      

      for ( unsigned i = 0 ; i < 3 ; ++i ) {
          if ( face_ptr->is_vert_nbr( cell_ptr_0->get_vert( i ) ) == false )
            vert_ptr = static_cast<c_vertex*>( cell_ptr_0->get_vert( i ) );
      }
      vert_ptr_0 = static_cast<c_vertex*>( face_ptr->get_vert( 0 ) );
      vert_ptr_1 = static_cast<c_vertex*>( face_ptr->get_vert( 1 ) );
      vec_0 = vert_ptr_0->get_point().get_pos() - vert_ptr->get_point().get_pos();
      vec_1 = vert_ptr_1->get_point().get_pos() - vert_ptr->get_point().get_pos();
      if ( vec_0.get_unit() * vec_1.get_unit() < 0.5 )
        continue;

      split_quad_to_tri( cell_ptr_1 );
      m_face_iter = collapse_face( m_face_iter );
      cond = true;
    } 
  }
	return;
}


int 
NIFS::c_mesh_adap::enhance_topolog( void )
{
  c_face* face_ptr;
//  c_vertex* vert_ptr;
  int n = 1;

  for ( int i = 0 ; ( i < 100 ) && ( n != 0 ) ; ++ i ) {
  n = 0;
  for ( m_face_iter = m_faces.begin(); 
        m_face_iter != m_faces.end(); 
        m_face_iter++ ) {
    face_ptr = dynamic_cast<c_face*>( *m_face_iter );
    if ( enhance_topolog( face_ptr ) > 0 )
      ++n;
  }
/*
  for ( m_vert_iter = m_verts.begin(); 
        m_vert_iter != m_verts.end(); 
        m_vert_iter++ ) {
    vert_ptr = dynamic_cast<c_vertex*>( *m_vert_iter );
    if ( enhance_topolog( vert_ptr ) > 0 )
      ++n;
  }
*/  
  }

  return 0;
}


void 
NIFS::c_mesh_adap::adapt()
{
  int n,k;
  double max_move = 1.0;
  
  update_solution();
  k=0;	
  while ( ( k < m_info.get_params()->get_adap_iter() ) && 
          ( max_move > m_info.get_params()->get_adap_trsh()  ) ) {
    std::cout << k << std::endl;
  	do {
      n = get_vert_num();
      modify();
      delaunay();
    } while ( n != get_vert_num() );
//    std::cout << i << ' ' << k << ' ' << n << ' ' << max_move << std::endl;
    update_geometry();
    max_move = smooth();
    delaunay();
    update_geometry();
    k++;
  } 
//  switch_quad();
//  cleanup_tri();  
}
