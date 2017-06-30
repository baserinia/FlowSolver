// Incompressible Flow Solver

// Copyright (C) 2005  Amir R. Baserinia

//----------------------------------------------------------------------------
// Created on:	17 Sep 2004
// Last update:	30 Jan 2006
//-----------------------------------------------------------------------------

#include <vector>
#include <algorithm>
#include <iostream>
#include "mesh_move.h"
#include "face_int.h"
#include "vector_2d.h"

#define PI 3.1415926535

//--- IMPLEMENTATION ---


// This a simple smoother based on the face error
void 
NIFS::c_mesh_move::move_vertex()
{
	int i;
	int face_num;	// number of neighbor vertices
	c_vector_2d vec;
	double	k;
	double	sum_k;
	c_face*	face_ptr;
	bool 	cond = false;

	for (m_vert_iter = m_verts.begin(); m_vert_iter != m_verts.end(); ++m_vert_iter ) {
		// interior points
//		if ( (*m_vert_iter)->get_type() == 0) {
			face_num = (*m_vert_iter)->get_face_num();
			sum_k = 0.0;
			vec.set_vector( 0.0, 0.0 );
			cond = false;
			for ( i = 0 ; i < face_num ; ++i ) {
				face_ptr = static_cast<c_face*>( (*m_vert_iter)->get_face( i ) );
				if ( face_ptr->get_cell_num() == 1 ) {
					cond = true;
				}
				
				if ( face_ptr->get_cell_num() == 1 )
					k = static_cast< c_cell* >(face_ptr->get_cell( 0 ))->get_residual()(1);
				else
					k = ( static_cast< c_cell* >(face_ptr->get_cell( 0 ))->get_residual()(1) +
						  static_cast< c_cell* >(face_ptr->get_cell( 1 ))->get_residual()(1) 
						) / 2.0;
				
				k /= sqr( face_ptr->get_area() );
				sum_k += k;
				if ( face_ptr->get_vert( 0 ) ==  (*m_vert_iter) ) {
					vec = vec + ( static_cast<c_vertex*>( face_ptr->get_vert( 1 ) )->get_point().get_pos() -
								  static_cast<c_vertex*>( face_ptr->get_vert( 0 ) )->get_point().get_pos()		
								) * k;
				}
				else {
					vec = vec + ( static_cast<c_vertex*>( face_ptr->get_vert( 0 ) )->get_point().get_pos() -
								  static_cast<c_vertex*>( face_ptr->get_vert( 1 ) )->get_point().get_pos()		
								) * k;
				}
			}	
			vec /= sum_k;
			
			vec = vec * 0.05;// + (*m_vert_iter)->get_point().get_pos() * 0.9;
			vec = vec + (*m_vert_iter)->get_point().get_pos();
			
			unsigned vert_type = (*m_vert_iter)->get_type();
			if ( vert_type == 0 )
				(*m_vert_iter)->set_point( c_point_2d( vec ) );
			else if ( vert_type == 1 )
			{
				(*m_vert_iter)->move_vertex_bnd( vec );
//				std::cout << vec.get_x() << ' ' << vec.get_y() << std::endl;
			}
	}
}


// 	swaps a face
//	return valuse:
//		-1:	boundary face
//		-2: one or both neighbouring cells are not triangle
//		-3:	the associated quadrilateral is concave
//		 n: number of swaped faces ( n >= 0 )
//
int 
NIFS::c_mesh_move::swap_face( c_face* face )
{
	int i;

	// ignore boundary faces
	if ( face->get_cell_num() < 2 ) 
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
	for ( i=0 ; i < 3 ; ++i) {
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
	for ( i=0; i<3; i++) {
		if  ( cell_temp[0]->get_face(i) != face ) {
			face_ptr = static_cast<c_face*>(cell_temp[0]->get_face(i));
			if ( static_cast<c_face*>(cell_temp[0]->get_face(i))->is_vert_nbr( bnd_vert_temp[0] ) )
				bnd_face_temp[0] =  static_cast<c_face*>( cell_temp[0]->get_face(i) );
			else
				bnd_face_temp[1] =  static_cast<c_face*>( cell_temp[0]->get_face(i) );
		}

		if  ( cell_temp[1]->get_face(i) != face ) {
			if ( static_cast<c_face*>(cell_temp[1]->get_face(i))->is_vert_nbr( bnd_vert_temp[0] ) )
				bnd_face_temp[2] =  static_cast<c_face*>( cell_temp[1]->get_face(i) );
			else
				bnd_face_temp[3] =  static_cast<c_face*>( cell_temp[1]->get_face(i) );
		}
	}

	// ignore concave quadrilaterals
	c_vector_2d vec_0 = bnd_vert_temp[1]->get_point().get_pos() - bnd_vert_temp[0]->get_point().get_pos();
	c_vector_2d vec_1 = bnd_vert_temp[2]->get_point().get_pos() - bnd_vert_temp[0]->get_point().get_pos();
	c_vector_2d vec_2 = bnd_vert_temp[3]->get_point().get_pos() - bnd_vert_temp[0]->get_point().get_pos();
	double theta_1 = vec_0.get_angle( vec_1 ) + vec_0.get_angle( vec_2 );
	if ( theta_1 > PI )
		return -3;
	// Delaunay
//	theta_1 =  std::acos( ( ( face->GetMetric() * vec_1 ) * vec_2 ) /
//							( sqrt( ( face->GetMetric() * vec_1 ) * vec_1 ) *
//							  sqrt( ( face->GetMetric() * vec_2 ) * vec_2 ) ) );

	theta_1 =  std::acos( (  vec_1  * vec_2 ) /
						  ( std::sqrt( vec_1 * vec_1 ) * 
						  	std::sqrt( vec_2 * vec_2 ) ) 
						);


	vec_0 = bnd_vert_temp[0]->get_point().get_pos() - bnd_vert_temp[1]->get_point().get_pos();
	vec_1 = bnd_vert_temp[2]->get_point().get_pos() - bnd_vert_temp[1]->get_point().get_pos();
	vec_2 = bnd_vert_temp[3]->get_point().get_pos() - bnd_vert_temp[1]->get_point().get_pos();
	double theta_2 = vec_0.get_angle( vec_1 ) + vec_0.get_angle( vec_2 );
	if ( theta_2 > PI )
		return -3;
	// Delaunay
	theta_2 =    std::acos( (  vec_1  * vec_2 ) /
						  ( std::sqrt( vec_1 * vec_1 ) * 
						  	std::sqrt( vec_2 * vec_2 ) ) 
						);

	if ( (theta_1 + theta_2) > PI )
		return 0;
		
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

//	m_iSwaface++;

	return 1;
}


void 
NIFS::c_mesh_move::swap()
{
	int swap_case = 1;
	for ( ; swap_case != 0; ) {
		swap_case = 0;
		for ( m_face_iter = m_faces.begin(); m_face_iter != m_faces.end(); m_face_iter++) 
			if ( swap_face( *m_face_iter ) == 1 )
				swap_case++;
	}
}


void 
NIFS::c_mesh_move::adapt()
{
	update_solution();
	for (int j=0; j<20; j++ ) { 
		move_vertex();
		update_geometry();
		swap();					
		for (int i=0; i<20; i++ )
			update_solution();
	}
}

