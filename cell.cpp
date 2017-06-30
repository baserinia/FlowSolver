//	Copyright (C) 2006  Amir R. Baserinia
//
// 	This file is part of SMAIF (Simple Mesh Adaptor for Incompressible Flow)
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
// Date:   12 Apr 2004
// Update: 27 Jun 2006
//----------------------------------------------------------------------------

#include <cmath>
#include <cstdlib>
#include "cell.h"
#include "face.h"
#include "vertex.h"
#include "list_nbr.h"
#include "vector_block.h"
#include "vector_2d.h"
#include <iostream>


//--- IMPLEMENTATION ---

//---( c_cell constructor )-----------------------------------------------------
NIFS::c_cell::c_cell()
{
	m_state		= 0.0;
	m_state_old = 0.0;
	m_grad		= c_vector_2d( 0.0 );
	m_hess		= c_matrix_2d( 0.0, 0.0, 0.0, 0.0 );
	m_coef_slf	= 0.0;
	m_coef_slf[1][1] = m_coef_slf[2][2] = 1.0;
	m_coef_nbr.clear();
}

//---( calculate volume and centroid )------------------------------------------
void 
NIFS::c_cell::calc_volume_and_cent()
{
	double area = 0.0;
	c_vector_2d centroid( 0.0 , 0.0 );
	double a;
	c_face* face_ptr;

	for ( unsigned i = 0 ; i < get_face_num() ; ++i ) {
		face_ptr = dynamic_cast<c_face*>( get_face( i ) );
		a = 0.5 * std::abs( 
            ( dynamic_cast<c_vertex*>( face_ptr->get_vert( 0 ) )->get_point( ).get_pos( ) - 
              dynamic_cast<c_vertex*>( get_vert( 0 ) )->get_point( ).get_pos( ) ) ^ 
            ( dynamic_cast<c_vertex*>( face_ptr->get_vert( 1 ) )->get_point( ).get_pos( ) - 
              dynamic_cast<c_vertex*>( get_vert( 0 ) )->get_point( ).get_pos( ) ) );
		area += a;
		centroid += ( 
              dynamic_cast<c_vertex*>( get_vert( 0 ) )->get_point( ).get_pos( ) + 
              dynamic_cast<c_vertex*>( face_ptr->get_vert( 0 ) )->get_point( ).get_pos( ) + 
              dynamic_cast<c_vertex*>( face_ptr->get_vert( 1 ) )->get_point( ).get_pos( ) 
                 ) * ( a / 3.0 );
	}
	m_node_point.set_pos( centroid / area );
	m_volume = area;	// volume per unit depth 
}

//---( update geometry )--------------------------------------------------------
void 
NIFS::c_cell::update_geometry( void )
{
	calc_volume_and_cent( );
}


//---( gradient reconstruction )------------------------------------------------
void 
NIFS::c_cell::calc_gradient( void )
{
	vector_block< double > A[2][2];
	vector_block< double > b[2];
	vector_block< double > result[2];
	double	x;
	double	y;
	vector_block< double > state;
	c_cell* cell_ptr;
	c_face* face_ptr;
	c_vertex* vert_ptr;
	vector_block< c_vector_2d > grad_temp;
	vector_block< c_vector_2d >	grad_old;
	unsigned int i;
	unsigned int j;
	unsigned int k;
	const vector_block< double > dummy( 1.0 );
	double weight;  

	grad_old = m_grad;
	if ( static_cast<int>( get_info().get_params()->get_grad_swc() ) != 1 )
		return;	//--- first order accurate solution
	
	int hess_acc = static_cast<int>( get_info().get_params()->get_hess_acc() );
	
	for ( i = 0 ; i < 2 ; ++i ) {
		b[i] = 0.0;
		for ( j = 0 ; j < 2 ; ++j )
			A[i][j] = 0.0;
	}

	double x_p	 = get_node().get_pos().get_x();
	double y_p	 = get_node().get_pos().get_y();
	vector_block< double > state_p = get_state();
	double power = get_info().get_params()->get_grad_pow();
	
	// contribution of immediate neighbours
	for ( i = 0 ; i < get_cell_num(); ++i ) {
		cell_ptr = dynamic_cast< c_cell* >( get_cell( i ) );
		x	= cell_ptr->get_node().get_pos().get_x() - x_p;
		y	= cell_ptr->get_node().get_pos().get_y() - y_p;
		state	= cell_ptr->get_state() - state_p;
		weight = std::pow( x * x + y * y, power / 2.0 );
		

		for ( j = 0 ; j < 4 ; ++j ) {
			if ( hess_acc == 1 ) {
				state[j] -= x * x * m_hess(j)(0)(0) / 2.0 + 
							x * y * m_hess(j)(0)(1) + 
							y * y * m_hess(j)(1)(1) / 2.0 ;
			}
			A[0][0][j] += x * x * weight;
			A[0][1][j] += x * y * weight;
			A[1][0][j] += y * x * weight;
			A[1][1][j] += y * y * weight;
			b[0][j] += state(j) * x * weight;
			b[1][j] += state(j) * y * weight; 
		}
	}

	// contribution of vertex neighbours
	if ( static_cast<int>( get_info().get_params()->get_grad_vrt() ) == 1 ) {
		for ( i = 0 ; i < get_vert_num(); ++i ) {
			vert_ptr = dynamic_cast< c_vertex* >( get_vert( i ) );
			for ( j = 0 ; j < vert_ptr->get_cell_num(); ++j ) {
				cell_ptr = dynamic_cast< c_cell* >( vert_ptr->get_cell( j ) );
				if ( is_cell_nbr( cell_ptr ) )
					continue;
				x	= cell_ptr->get_node().get_pos().get_x() - x_p;
				y	= cell_ptr->get_node().get_pos().get_y() - y_p;
				state	= cell_ptr->get_state() - state_p;
				weight = std::pow( x * x + y * y, power / 2.0 );
	
				for ( k = 0 ; k < 4 ; ++k ) {
					if ( hess_acc == 1 ) {
						state[k] -= x * x * m_hess(k)(0)(0) / 2.0 + 
									x * y * m_hess(k)(0)(1) + 
									y * y * m_hess(k)(1)(1) / 2.0 ;
					}
					A[0][0][k] += x * x * weight;
					A[0][1][k] += x * y * weight;
					A[1][0][k] += y * x * weight;
					A[1][1][k] += y * y * weight;
					b[0][k] += state(k) * x * weight;
					b[1][k] += state(k) * y * weight; 
				}
			}
		}
	}

	// contribution of boundary conditions
	if ( get_face_num( ) != get_cell_num( ) ) {
		for ( i = 0 ; i < get_face_num(); ++i ) {
			face_ptr = dynamic_cast< c_face* >( get_face( i ) );
			if ( face_ptr->get_cell_num() == 1 ) { // detects boundary face
				x = face_ptr->get_fmp().get_pos().get_x() - x_p;
				y = face_ptr->get_fmp().get_pos().get_y() - y_p;
				state = face_ptr->get_state( ) - state_p;
				weight = std::pow( x * x + y * y, power / 2.0 );

				for ( j = 0 ; j < 4 ; ++j ) {
					if ( hess_acc == 1 ) {
						state[j] -= x * x * m_hess(j)(0)(0) / 2.0 + 
									x * y * m_hess(j)(0)(1) + 
									y * y * m_hess(j)(1)(1) / 2.0 ;
					}
					A[0][0][j] += face_ptr->get_weight( )(j) * ( x * x ) * weight;
					A[0][1][j] += face_ptr->get_weight( )(j) * ( x * y ) * weight;
					A[1][0][j] += face_ptr->get_weight( )(j) * ( y * x ) * weight;
					A[1][1][j] += face_ptr->get_weight( )(j) * ( y * y ) * weight;
					b[0][j] += face_ptr->get_weight( )(j) * state(j) * x * weight;
					b[1][j] += face_ptr->get_weight( )(j) * state(j) * y * weight;
				}
			}
		}    
	}

	double det; 
	for ( i = 0 ; i < 4 ; ++i ) {
		det = A[0][0](i) * A[1][1](i) - A[0][1](i) * A[1][0](i);
		if ( std::fabs( det ) > 1E-20 ) {
			result[0][i] = ( b[0](i) * A[1][1](i) - b[1](i) * A[0][1](i) ) / det;
			result[1][i] = ( b[1](i) * A[0][0](i) - b[0](i) * A[1][0](i) ) / det;
		} 
		else {
			vector_block< double > v = b[0] / 
						std::sqrt( A[0][0](i)*A[0][0](i) + A[0][1](i)*A[0][1](i) );
			result[0][i] = v(i) * A[0][0](i);
			result[1][i] = v(i) * A[0][1](i);
		}
	}
	
	double grad_rlx = get_info().get_params()->get_grad_rlx();
	for ( i = 0; i < 4; ++i ) {
		grad_temp[i] = c_vector_2d( result[0](i) , result[1](i) );
		m_grad[i] = grad_temp(i) * grad_rlx + grad_old(i) * ( 1.0 - grad_rlx );
	}
}


//---( Hessian reconstruction )-----------------------------------------------
void 
NIFS::c_cell::calc_hessian()
{
	vector_block< double > A[3][3];
	vector_block< double > b[3];
	vector_block< double > result[3];
	double	x;
	double	y;
	vector_block< c_vector_2d > grad;
	c_cell* cell_ptr;
//	c_face* face_ptr;
	c_vertex* vert_ptr;
	vector_block< c_matrix_2d > hess_temp;
	vector_block< c_matrix_2d > hess_old;
	unsigned int i;
	unsigned int j;
	unsigned int k;
	const vector_block< double > dummy( 1.0 );
	double weight;  

	hess_old = m_hess;

	if ( static_cast<int>( get_info().get_params()->get_hess_swc() ) != 1 )
		return;	//--- first order accurate solution

	for (i=0; i<3; i++) {
		b[i] = 0.0;
		for (j=0; j<3; j++)
			A[i][j] = 0.0;
	}

	double x_p	 = get_node().get_pos().get_x();
	double y_p	 = get_node().get_pos().get_y();
	vector_block< c_vector_2d > grad_p = get_grad();
	double power = get_info().get_params()->get_hess_pow();

	// contribution of immediate neighbours
	for ( i = 0 ; i < get_cell_num(); ++i ) {
		cell_ptr = dynamic_cast< c_cell* >( get_cell( i ) );
		x =  cell_ptr->get_node().get_pos().get_x() - x_p;
		y =  cell_ptr->get_node().get_pos().get_y() - y_p;
		grad = cell_ptr->get_grad() - grad_p;
		weight = std::pow( x * x + y * y, power / 2.0 );

		for ( j = 0 ; j < 4 ; ++j ) {
			
			A[0][0][j] += x * x * weight;
			A[0][1][j] += x * y * weight;
			A[1][0][j] += x * y * weight;
			A[1][1][j] += ( x * x + y * y ) * weight;
			A[1][2][j] += x * y * weight;
			A[2][1][j] += x * y * weight;
			A[2][2][j] += y * y * weight;
			b[0][j] += grad(j).get_x() * x * weight;
			b[1][j] += ( grad(j).get_x() * y + grad(j).get_y() * x ) * weight;
			b[2][j] += grad(j).get_y() * y * weight;
		}
	}

	// contribution of vertex neighbours
	if ( static_cast<int>( get_info().get_params()->get_hess_vrt() ) == 1 ) {
		for ( i = 0 ; i < get_vert_num(); ++i ) {
			vert_ptr = dynamic_cast< c_vertex* >( get_vert( i ) );
			for ( j = 0 ; j < vert_ptr->get_cell_num(); ++j ) {
				cell_ptr = dynamic_cast< c_cell* >( vert_ptr->get_cell( j ) );
				if ( is_cell_nbr( cell_ptr ) )
					continue;

				x =  cell_ptr->get_node().get_pos().get_x() - x_p;
				y =  cell_ptr->get_node().get_pos().get_y() - y_p;
				grad = cell_ptr->get_grad() - grad_p;
				weight = std::pow( x * x + y * y, power / 2.0 );
	
				for ( k = 0 ; k < 4 ; ++k ) {
					A[0][0][k] += x * x * weight;
					A[0][1][k] += x * y * weight;
					A[1][0][k] += x * y * weight;
					A[1][1][k] += ( x * x + y * y ) * weight;
					A[1][2][k] += x * y * weight;
					A[2][1][k] += x * y * weight;
					A[2][2][k] += y * y * weight;
					b[0][k] += grad(k).get_x() * x * weight;
					b[1][k] += ( grad(k).get_x() * y + grad(k).get_y() * x ) * weight;
					b[2][k] += grad(k).get_y() * y * weight;
				}
			}
		}
	}

	double det; 
	for ( i = 0 ; i < 4 ; ++i ) {
		det = A[0][0](i) * A[1][1](i) * A[2][2](i) +
			  A[0][1](i) * A[1][2](i) * A[2][0](i) +
			  A[0][2](i) * A[1][0](i) * A[2][1](i) -
			  A[0][0](i) * A[1][2](i) * A[2][1](i) -
			  A[0][1](i) * A[1][0](i) * A[2][2](i) -
			  A[0][2](i) * A[1][1](i) * A[2][0](i) ;

		if ( std::fabs(det) > 1.0E-100 ) {

			// Gauss Elimination, Forward Sweep
			double d = A[1][0](i) / A[0][0](i);
			for (j=0; j<3; ++j) {
				A[1][j][i] -= A[0][j](i) * d;
			}
			b[1][i] -= b[0](i) * d;

			d = A[2][1](i) / A[1][1](i);
			for (j=1; j<3; ++j) {
				A[2][j][i] -=  A[1][j](i) * d;
			}
			b[2][i] -= b[1](i) * d;

			// Back Substitution
			result[2][i] = b[2](i) / A[2][2](i);
			result[1][i] = ( b[1](i) - A[1][2](i) * result[2](i) ) / A[1][1](i);
			result[0][i] = ( b[0](i) - A[0][2](i) * result[2](i) - A[0][1](i) * result[1][i]) / A[0][0](i);
		} else {
			result[2][i] = 0.0;
			result[1][i] = 0.0;
			result[0][i] = 0.0;
		}
	}

	double hess_rlx = get_info().get_params()->get_hess_rlx();
	for ( i = 0; i < 4; ++i ) {
		hess_temp[i].set_matrix( result[0](i) , result[1](i), result[1](i), result[2](i) );
		m_hess[i] = hess_temp(i) * hess_rlx + hess_old(i) * ( 1.0 - hess_rlx );
	}

}


//---( calculate coefficients of the control volume equation )------------------
void 
NIFS::c_cell::calc_coefs( void )
{
	c_face* face_ptr;
	c_cell* cell_ptr;

	m_coef_slf = 0.0;
	m_coef_rhs = 0.0;
	m_resid = 0.0;	// new
	m_coef_nbr.clear();
	m_coef_nbr.resize( get_cell_num(),0.0 );

	for ( unsigned i = 0; i < get_face_num(); ++i ) {
		face_ptr = dynamic_cast< c_face* >( get_face( i ) );
		if ( this == face_ptr->get_cell( 0 ) ) {
			m_coef_slf += face_ptr->get_coef_1();
			m_coef_rhs -= face_ptr->get_coef_c();
			m_resid    += face_ptr->get_error();
		}
		else {
			m_coef_slf -= face_ptr->get_coef_2();
			m_coef_rhs += face_ptr->get_coef_c();
			m_resid    -= face_ptr->get_error();
		}

		if ( face_ptr->get_cell_num() == 2 ) {
			for ( unsigned j = 0 ; j < get_cell_num(); ++j ) {
				cell_ptr = dynamic_cast< c_cell* >( get_cell( j ) );
				if ( cell_ptr == face_ptr->get_cell( 0 ) ) 
					m_coef_nbr[j] = -face_ptr->get_coef_1();
				else if ( cell_ptr == face_ptr->get_cell( 1 ) ) 
					m_coef_nbr[j] = face_ptr->get_coef_2();
			} //--- for
		} //--- if
	} //--- for
	
	// buoyancy terms
	m_coef_slf[1][3] += get_info().get_params()->get_rho() * 
						get_info().get_params()->get_beta() * 
						get_info().get_params()->get_g_x() * 
						get_volume();
	m_coef_rhs[1]	 += get_info().get_params()->get_rho() * 
						get_info().get_params()->get_beta() * 
						get_info().get_params()->get_g_x() * 
						get_info().get_params()->get_ref_temp() * 
						get_volume();
	m_coef_slf[2][3] += get_info().get_params()->get_rho() * 
						get_info().get_params()->get_beta() * 
						get_info().get_params()->get_g_y() * 
						get_volume();
	m_coef_rhs[2]	 += get_info().get_params()->get_rho() * 
						get_info().get_params()->get_beta() * 
						get_info().get_params()->get_g_y() * 
						get_info().get_params()->get_ref_temp() * 
						get_volume();
}

//---( calculate and retrieve residual )----------------------------------------
NIFS::vector_block< double > 
NIFS::c_cell::get_residual()
{
	vector_block< double > residual( 0.0 );
	for ( int i = 0 ; i < 4 ; ++i )
		residual[i] = fabs( m_state(i) - m_state_old(i) );
	return residual;
}

//---( prepare results for printing )-------------------------------------------
void 
NIFS::c_cell::print_result( void )
{
	c_face* face_ptr;
	c_vertex* vert_ptr;
	vector_block< double > dummy = 0.0;

  if ( get_info().get_params()->get_soln_out() != 0.0 )
    get_info().get_res_file()->print_line_soln( get_node().get_pos() , m_state );

  if ( get_info().get_params()->get_grad_out() != 0.0 )
    get_info().get_res_file()->print_line_grad( get_node().get_pos() , m_grad );

  if ( get_info().get_params()->get_hess_out() != 0.0 )
    get_info().get_res_file()->print_line_hess( get_node().get_pos() , m_hess );

  if ( get_info().get_params()->get_solb_out() != 0.0 )
    get_info().get_res_file()->print_line_solb( get_node().get_pos() , m_state );

  if ( get_info().get_params()->get_geom_out() != 0.0 )
    get_info().get_res_file()->print_line_geom( get_node().get_pos() , get_volume() );

  if ( get_info().get_params()->get_rsdl_out() != 0.0 )
    get_info().get_res_file()->print_line_rsdl( get_node().get_pos() , m_resid );


  // if the cell has no boundary face or the flag is not set
  if ( ( get_cell_num() == get_face_num() ) ||
       ( get_info().get_params()->get_bnd_out() == 0.0 ) )
    return;
		
  for ( unsigned i = 0 ; i < get_face_num() ; ++i ) {
    face_ptr = static_cast< c_face* >( get_face( i ) );
    if ( face_ptr->get_cell_num() == 1 ) {
      if ( get_info().get_params()->get_soln_out() != 0.0 )
        get_info().get_res_file()->print_line_soln( face_ptr->get_fmp().get_pos() , 
                                                    face_ptr->get_state() );

      if ( get_info().get_params()->get_grad_out() != 0.0 )
		   get_info().get_res_file()->print_line_grad( face_ptr->get_fmp().get_pos() , 
                                                       face_ptr->get_grad() );

      if ( get_info().get_params()->get_hess_out() != 0.0 )
           get_info().get_res_file()->print_line_hess( face_ptr->get_fmp().get_pos() , 
                                                       face_ptr->get_hess() );

      if ( get_info().get_params()->get_solb_out() != 0.0 )
           get_info().get_res_file()->print_line_solb( face_ptr->get_fmp().get_pos() , 
                                                       face_ptr->get_state() );
    }
  }
  
  // if cell has corner vertices
  for ( unsigned i = 0 ; i < get_vert_num() ; ++i ) {
    vert_ptr = static_cast< c_vertex* >( get_vert( i ) );
    if ( vert_ptr->get_type() == 2 ) {
      vector_block< double > state;
      c_vector_2d r = vert_ptr->get_point().get_pos() - get_node().get_pos();
	    for ( int j = 0 ; j < 4 ; ++j )
		    state[j] = get_state()(j) + get_grad()(j) * r;
      
      if ( get_info().get_params()->get_soln_out() != 0.0 )
        get_info().get_res_file()->print_line_soln( vert_ptr->get_point().get_pos() , state );

      if ( get_info().get_params()->get_grad_out() != 0.0 )
		   get_info().get_res_file()->print_line_grad( vert_ptr->get_point().get_pos() , get_grad() );

      if ( get_info().get_params()->get_hess_out() != 0.0 )
           get_info().get_res_file()->print_line_hess( vert_ptr->get_point().get_pos() , get_hess() );

      if ( get_info().get_params()->get_solb_out() != 0.0 )
           get_info().get_res_file()->print_line_solb( vert_ptr->get_point().get_pos() , state );
    }
  }
  
}

//---( set state of all variables )---------------------------------------------
void 
NIFS::c_cell::set_state( const vector_block< double >& state ) 
{ 
  double state_rlx = get_info().get_params()->get_soln_rlx();
	m_state_old = m_state;
	m_state = state * state_rlx + m_state_old * ( 1.0 - state_rlx);
} 

void 
NIFS::c_cell::set_force_state( const vector_block< double >& state ) 
{ 
	m_state = state;
} 

//---( set state for some variables )-------------------------------------------
void 
NIFS::c_cell::set_state( const vector_block< double >& state, 
						 unsigned index_beg, unsigned index_end)
{ 
  double state_rlx = get_info().get_params()->get_soln_rlx();
	m_state_old = m_state;
	for ( unsigned i = index_beg; i < index_end; ++i )
		m_state[i] = state(i) * state_rlx + m_state_old(i) * ( 1.0 - state_rlx);
} 

void 
NIFS::c_cell::set_force_state( const vector_block< double >& state, 
						 unsigned index_beg, unsigned index_end)
{ 
  m_state_old = m_state;
	for ( unsigned i = index_beg; i < index_end; ++i )
		m_state[i] = state(i);
} 


bool 
NIFS::c_cell::is_valid_move( c_mesh_obj* vert , c_vector_2d vec )
{
  bool cond = true;
	c_vertex* vert_ptr = static_cast< c_vertex* >( vert );
	c_face* face_ptr = NULL;
	for ( unsigned i = 0 ; i < get_face_num() ; ++i ) {
		face_ptr = dynamic_cast< c_face* >( get_face( i ) );
		if ( face_ptr->is_vert_nbr( vert_ptr ) == false ) {
		 if ( ( ( face_ptr->get_fmp().get_pos() - vert_ptr->get_point().get_pos() ) * 
		        ( face_ptr->get_unv() ) ) *
          ( ( face_ptr->get_fmp().get_pos() - vec ) *
            ( face_ptr->get_unv() ) ) < 0.0 )
      		cond = false;
		}
	}
	
	return cond;
}

double 
NIFS::c_cell::get_charac_size( void )
{
	return 1.5196714 * std::sqrt( get_volume() ); // sqrt( 4*volume / sqrt(3) )
}

NIFS::vector_block< double > 
NIFS::c_cell::get_err_iso_coef( void )
{
  vector_block< double > vec = 1.0;  
	
  vec[0] = std::sqrt( 
    sqr(-m_hess(1)(0)(0) + 
         m_hess(1)(1)(1) +
         4 * m_hess(1)(1)(0) +
         m_hess(2)(0)(0) -
         m_hess(2)(1)(1)
       ) +
    sqr( m_hess(2)(0)(0) - 
         m_hess(2)(1)(1) +
         4 * m_hess(2)(1)(0) -
         m_hess(1)(0)(0) +
         m_hess(1)(1)(1)
       ) ) / 32.0;

  vec[1] = 0.75 * std::sqrt( 
    sqr( get_info().get_params()->get_rho() * 
         ( m_grad(1)(1) * m_grad(2)(1) -
           m_grad(1)(0) * m_grad(2)(0) +
           2.0 * m_grad(2)(0) * m_grad(2)(1) 
         ) +
         2.0 * m_hess(0)(0)(1) 
       ) +
    sqr( get_info().get_params()->get_rho() * 
         (-m_grad(2)(1) * m_grad(2)(1) +
           m_grad(1)(1) * m_grad(2)(0) +
           m_grad(1)(0) * m_grad(2)(1) +
           m_grad(2)(0) * m_grad(2)(0) 
         ) +
         m_hess(0)(0)(0) - 
         m_hess(0)(1)(1)
       ) 
	);

  vec[2] = 0.75 * std::sqrt( 
    sqr( get_info().get_params()->get_rho() * 
         (-m_grad(1)(0) * m_grad(1)(0) +
           m_grad(1)(1) * m_grad(2)(0) +
           m_grad(1)(0) * m_grad(2)(1) +
           m_grad(1)(1) * m_grad(1)(1)  
         ) +
         m_hess(0)(1)(1) - m_hess(0)(0)(0)
       ) +
    sqr( get_info().get_params()->get_rho() * 
         ( m_grad(1)(0) * m_grad(2)(0) -
           m_grad(1)(1) * m_grad(2)(1) +
           2.0 * m_grad(1)(0) * m_grad(1)(1) 
         ) +
         m_hess(0)(0)(1)  
       ) 
	);
	   
	return vec;
}


double 
NIFS::c_cell::get_aspect_ratio()
{
  double prm = 0.0;

	for ( unsigned i = 0 ; i < get_face_num() ; ++i ) {
    prm += dynamic_cast<c_face*>( get_face( i ) )->get_area();
    
  }
  
  if ( get_face_num() == 3 )
    return sqr( prm ) / ( 12.0 * sqrt( 3.0 ) * get_volume() );
  else if ( get_face_num() == 4 )
    return sqr( prm ) / ( 16.0 * get_volume() );
  else
    return 1.0;
}


void 
NIFS::c_cell::calc_anis_spacing( void )
{
   m_ans_spc[0] = c_matrix_2d( 0.0, 0.0, 0.0, 0.0 );
   m_ans_spc[1] = get_anis_spacing_memu( );
   m_ans_spc[2] = c_matrix_2d( 0.0, 0.0, 0.0, 0.0 );      
   m_ans_spc[3] = c_matrix_2d( 0.0, 0.0, 0.0, 0.0 );
}


NIFS::c_matrix_2d
NIFS::c_cell::get_anis_spacing_memu( void )
{
  double e;
  double k; // aspect ratio
  double t; // orientation
  double rho = get_info().get_params()->get_rho();
  
  e = 4.0 * sqr( 
              rho * sqr(k) * ( 
                m_state(1) * m_hess(1)(0)(0) * pow( cos(t) , 3.0 ) + 
					      ( 2.0 * m_state(1) * m_hess(1)(0)(1) + 
					        m_state(2) * m_hess(1)(0)(0) ) * pow( cos(t) , 2.0 
					      ) * sin(t) + 
                ( m_state(1) * m_hess(1)(1)(1) + 
	 				        2.0 * m_state(2)*m_hess(1)(0)(1) 
                ) * cos(t) * pow( sin(t) , 2.0 ) + 
                  m_state(2) * m_hess(1)(1)(1) * pow( sin(t) , 3.0 ) ) +
              rho * ( 1.0 / 3.0 ) * ( 
                ( m_grad(1)(1) * m_grad(1)(1) + 
                  m_state(1) * m_hess(1)(1)(1) 
                ) * pow( cos(t) , 3.0 ) + 
                ( m_grad(1)(1) * m_grad(2)(1) - 
                  2.0 * m_grad(1)(0) * m_grad(1)(1) + 
                  m_state(2) * m_hess(1)(1)(1) - 
                  2.0 * m_state(1) * m_hess(1)(0)(1) 
                ) * pow( cos(t) , 2.0 ) * sin(t) +
                ( -m_grad(1)(1) * m_grad(2)(0) - 
                  m_grad(1)(0) * m_grad(2)(1) + 
                  m_grad(1)(0) * m_grad(1)(0) + 
                  m_state(1) * m_hess(1)(0)(0) -
                  2.0 * m_state(2) * m_hess(1)(0)(1)
                ) * cos(t) * pow( sin(t) , 2.0 ) +
                ( m_grad(1)(0) * m_grad(2)(0) + 
                  m_state(2) * m_hess(1)(0)(0)
                ) * pow( sin(t) , 3.0 )
              )+
              ( 1.0 / 6.0 ) * cos(t) *
                ( m_hess(0)(0)(0) * pow( sin(t) , 2.0 ) - 
                  2.0 * m_hess(0)(0)(1) * cos(t) * sin(t) + 
                  m_hess(0)(1)(1) * pow( cos(t) , 2.0 )
                ) -
              ( 1.0 / 2.0 ) * sqr(k) * cos(t) *
                ( m_hess(0)(0)(0) * pow( cos(t) , 2.0 ) + 
                  2.0 * m_hess(0)(0)(1) * cos(t) * sin(t) + 
                  m_hess(0)(1)(1) * pow( sin(t) , 2.0 ) 
                )
            ) +
      4.0 * sqr(k) * sqr( 
              rho * ( 
                m_state(2) * m_hess(1)(1)(1) * pow( cos(t) , 3.0 ) +
                ( -m_state(1) * m_hess(1)(1)(1) - 
                  2.0 * m_state(2) * m_hess(1)(0)(1) 
                ) * pow( cos(t) , 2.0 ) * sin(t) +
                ( m_state(2) * m_hess(1)(0)(0) + 
                  2.0 * m_state(1) * m_hess(1)(0)(1) 
                ) * cos(t) * pow( sin(t) , 2.0 ) +
                ( -m_state(1) * m_hess(1)(0)(0) ) * pow( sin(t) , 3.0 )
              ) + 
              rho * ( 1.0 / 3.0 ) * sqr(k) * ( 
                ( m_grad(1)(0) * m_grad(2)(0) + 
                  m_state(2)*m_hess(1)(0)(0)
                ) * pow( cos(t) , 3.0 ) +
                ( m_grad(1)(0) * m_grad(2)(1) - 
                  m_grad(1)(0) * m_grad(1)(0) + 
                  m_grad(1)(1) * m_grad(2)(0) - 
                  m_state(1) * m_hess(1)(0)(0) + 
                  2.0 * m_state(2) * m_hess(1)(0)(1)
                ) * pow( cos(t) , 2.0 ) * sin(t) +
                ( -2.0 * m_grad(1)(0) * m_grad(1)(1) + 
                  m_grad(1)(1) * m_grad(2)(1) + 
                  m_state(2) * m_hess(1)(1)(1) - 
                  2.0 * m_state(1) * m_hess(1)(0)(1)
                ) * cos(t) * pow( sin(t) , 2.0 ) +
                ( -m_grad(1)(1) * m_grad(1)(1) - 
                  m_state(1) * m_hess(1)(1)(1)
                ) * pow( sin(t) , 3.0 )
              ) - 
            ( 1.0 / 6.0 ) * sqr(k) * sin(t) *
              ( m_hess(0)(0)(0) * pow( cos(t) , 2.0 ) + 
                2.0 * m_hess(0)(0)(1) * cos(t) * sin(t) + 
                m_hess(0)(1)(1) * pow( sin(t) , 2.0 ) 
              ) + 
            ( 1.0 / 2.0 ) * sin(t) *
              ( m_hess(0)(0)(0) * pow( sin(t) , 2.0 ) - 
                2.0 * m_hess(0)(0)(1) * cos(t) * sin(t) + 
                m_hess(0)(1)(1) * pow( cos(t) , 2.0 ) 
              )
            );
      return c_matrix_2d( 0.0, 0.0, 0.0, 0.0 );
}
