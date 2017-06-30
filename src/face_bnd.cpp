// Incompressible Flow Solver

// Copyright (C) 2005  Amir R. Baserinia

//----------------------------------------------------------------------------
// Date:   10 Apr 2004
// Update: 01 Apr 2005
//-----------------------------------------------------------------------------

#include <cstdlib>
#include <iostream>
#include "face_bnd.h"
#include "cell.h"
#include "vertex.h"


//--- IMPLEMENTATION ---


void 
NIFS::c_face_bnd::update_geometry( void )
{
	c_face::update_geometry();
	calc_bc();
}

void 
NIFS::c_face_bnd::calc_uv12()
{
	m_uv12 = ( m_fmp.get_pos() - 
			   static_cast<c_cell*>( get_cell(0) )->get_node().get_pos()
			 ).get_unit();    
}

void 
NIFS::c_face_bnd::calc_s12()
{
    m_s12 = ( m_fmp.get_pos() - 
               static_cast<c_cell*>(get_cell(0))->get_node().get_pos()
             ).get_magnit();    
}

void 
NIFS::c_face_bnd::calc_df()
{
	c_cell* cell_ptr_1 = static_cast<c_cell*>( get_cell(0) );

	m_df =  cell_ptr_1->get_volume() / cell_ptr_1->get_df();
}


double 
NIFS::c_face_bnd::project( const c_vector_2d& pnt, c_vector_2d& prj )
{
	// distance along the face line
	double xi = ( pnt - static_cast<c_vertex*>( get_vert( 0 ) )->get_point().get_pos() ) * 
				get_utv(); 	
	double eta;	// distance from the face 
	double len = get_area(); // face length ( area in 2-D problems )

	c_vector_2d tan_v1 = static_cast<c_vertex*>( get_vert( 0 ) )->get_boundary_tangent();
	c_vector_2d tan_v2 = static_cast<c_vertex*>( get_vert( 1 ) )->get_boundary_tangent();
	double vt1 = 1.0; //tan_v1 * get_utv();
	double vn1 = 0.0; //tan_v1 * get_unv();
	double vt2 = 1.0; //tan_v2 * get_utv();
	double vn2 = 0.0; //tan_v2 * get_unv();
	
	double A1 = 0.0;
	double A2 = 0.0;
	double A3 = 0.0;
	
	bool cnd1 = ( std::abs( vn1 ) < std::abs( vt1 ) ) && ( 0.0 < std::abs( vt1 ) );
	bool cnd2 = ( std::abs( vn2 ) < std::abs( vt2 ) ) && ( 0.0 < std::abs( vt2 ) );
	
	if ( cnd1 && cnd2 )
	{
		A1 =  vn1 / vt1;
		A2 =  -( 2 * vn1 / vt1 + vn2 / vt2 ) / len;
		A3 =  ( vn1 / vt1 + vn2 / vt2 ) / ( len * len );
	}
	else if ( cnd1 )
	{
		A1 =  -0.5 * vn2 / vt2;
		A2 =  0.0;
		A3 =  0.5 * ( vn2 / vt2 ) / ( len * len );
	}
	else if ( cnd2 )
	{
		A1 =  vn1 / vt1;
		A2 =  -1.5 * ( vn1 / vt1 ) / len;
		A3 =  0.5 * ( vn1 / vt1 ) / ( len * len );
	}
	else
	{
		A1 = A2 = A3 = 0.0;
	}
	
	eta = A1 * ( xi ) + A2 * ( xi * xi ) + A3 * ( xi * xi * xi );
	
	prj = static_cast<c_vertex*>( get_vert( 0 ) )->get_point().get_pos() +
		  get_utv() * xi +
		  get_unv() * eta;
		  
	return xi/len;
}

