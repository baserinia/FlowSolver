// Incompressible Flow Solver

// Copyright (C) 2005  Amir R. Baserinia

//-----------------------------------------------------------------------------
// Date:   10 Apr 2004
// Update: 25 Feb 2005
//-----------------------------------------------------------------------------

#include "face.h"
#include "vertex.h"
#include "cell.h"
#include <iostream>

//--- IMPLEMENTATION ---

NIFS::c_face::c_face()
{
	m_state		= 0.0;
	m_grad		= c_vector_2d( 0.0 );
	m_coef_1	= 0.0;
	m_coef_2	= 0.0;
	m_coef_c	= 0.0;
}


void 
NIFS::c_face::update_geometry( void )
{
	calc_area();
	calc_utv();
	calc_fmp();
	calc_r1();
	calc_unv();
	calc_r2();		// pure virtual
	calc_rc();		// pure virtual
	calc_uv12();	// pure virtual
	calc_s12();		// pure virtual
	calc_alpha();	// pure virtual
//	std::cout << m_area << std::endl;
}

void 
NIFS::c_face::update_solution( void )
{
  calc_df();
  calc_state();
  calc_mfr();
  calc_coef_1();
  calc_coef_2();
  calc_coef_c();
  calc_error(); // new
}

void
NIFS::c_face::calc_area()
{
	m_area = ( dynamic_cast<c_vertex*>( get_vert(1) )->get_point().get_pos() - 
			   dynamic_cast<c_vertex*>( get_vert(0) )->get_point().get_pos()
			 ).get_magnit();
}

void 
NIFS::c_face::calc_utv()
{
  m_utv = ( dynamic_cast<c_vertex*>( get_vert(1) )->get_point().get_pos() - 
            dynamic_cast<c_vertex*>( get_vert(0) )->get_point().get_pos()
          ).get_unit();
}

void 
NIFS::c_face::calc_unv()
{
    c_vector_2d vec = m_utv.rotate_ccw();
	m_unv = ( vec * m_r1 ) > 0.0 ? vec : -vec;
}

void 
NIFS::c_face::calc_fmp()
{
	m_fmp.set_pos( ( dynamic_cast<c_vertex*>( get_vert(1) )->get_point().get_pos() + 
					 dynamic_cast<c_vertex*>( get_vert(0) )->get_point().get_pos() 
				  ) / 2.0 );
}

void 
NIFS::c_face::calc_r1()
{
	m_r1 = m_fmp.get_pos() - 
		   dynamic_cast<c_cell*>( get_cell(0) )->get_node().get_pos();
}

void
NIFS::c_face::calc_alpha()
{
	m_alpha = m_uv12 * m_unv; //  alpha = s_hat . n_hat
}

NIFS::c_vector_2d 
NIFS::c_face::get_collapse_pos( void )
{
	c_vector_2d vec;
	unsigned type0 = dynamic_cast< c_vertex* >( get_vert( 0 ) )->get_type();
	unsigned type1 = dynamic_cast< c_vertex* >( get_vert( 1 ) )->get_type();
	
	if ( ( type0 == type1 ) && ( type0 > 0 ) )
		project( get_fmp().get_pos() , vec );
	else if ( type0 == type1 )
		vec = get_fmp().get_pos();
	else if ( type0 > type1 )
		vec = dynamic_cast< c_vertex* >( get_vert( 0 ) )->get_point().get_pos();
	else
		vec = dynamic_cast< c_vertex* >( get_vert( 1 ) )->get_point().get_pos();
	
	return vec;
}

double 
NIFS::c_face::get_size( void )
{
	double d = (get_fmp().get_pos() - c_vector_2d(0.5,0.5)).get_magnit();

	return 0.05-0.04*exp(-200.0*(d-0.25)*(d-0.25));

}
