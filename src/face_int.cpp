// Incompressible Flow Solver

// Copyright (C) 2005  Amir R. Baserinia

//------------------------------------------------------------------------------
// Date:   10 Apr 2004
// Update: 31 Mar 2005
//------------------------------------------------------------------------------

#include "face_int.h"
#include "cell.h"
#include "params.h"
#include "vector_block.h"
#include "vector_2d.h"

//--- IMPLEMENTATION ---


//---( Calculate unit vector from node 1 to node 2 )----------------------------
void 
NIFS::c_face_int::calc_uv12( void )
{
    m_uv12 = ( static_cast<c_cell*>( get_cell(1) )->get_node().get_pos() - 
			   static_cast<c_cell*>( get_cell(0) )->get_node().get_pos() 
			 ).get_unit();
}

//---( Calculate distance between node 1 to node 2 )----------------------------
void 
NIFS::c_face_int::calc_s12( void )
{
	m_s12 = ( static_cast<c_cell*>( get_cell(1) )->get_node().get_pos() - 
			  static_cast<c_cell*>( get_cell(0) )->get_node().get_pos() 
			).get_magnit();
}

void 
NIFS::c_face_int::calc_r2()
{
	m_r2 = m_fmp.get_pos() - 
		   static_cast<c_cell*>( get_cell(1) )->get_node().get_pos();
}

void 
NIFS::c_face_int::calc_rc()
{
    m_rc = m_fmp.get_pos() - 
		  ( static_cast<c_cell*>( get_cell(0) )->get_node().get_pos() + 
			static_cast<c_cell*>( get_cell(1) )->get_node().get_pos() 
		  ) / 2.0; 
}

void 
NIFS::c_face_int::calc_df()
{
	c_cell* cell_ptr_1 = static_cast<c_cell*>( get_cell(0) );
	c_cell* cell_ptr_2 = static_cast<c_cell*>( get_cell(1) );

	m_df = 0.5 * ( cell_ptr_1->get_volume() / cell_ptr_1->get_df() + 
				   cell_ptr_2->get_volume() / cell_ptr_2->get_df() );
}

void 
NIFS::c_face_int::calc_coef_1()
{
	m_coef_1[0][0] = m_alpha * m_df * m_area / m_s12;
	m_coef_1[0][1] = 0.5 * m_unv.get_x() * m_area;
	m_coef_1[0][2] = 0.5 * m_unv.get_y() * m_area;
	m_coef_1[0][3] = 0.0;
	
	m_coef_1[1][0] = 0.5 * m_unv.get_x() * m_area;
	m_coef_1[1][1] = get_info().get_params()->get_mu() * 
					 m_alpha * m_area / m_s12 +
					 ( m_mfr > 0.0 ? m_mfr : 0.0 );
	m_coef_1[1][2] = 0.0;
	m_coef_1[1][3] = 0.0;
	
	m_coef_1[2][0] = 0.5 * m_unv.get_y() * m_area;
	m_coef_1[2][1] = 0.0;
	m_coef_1[2][2] = get_info().get_params()->get_mu() * 
					 m_alpha * m_area / m_s12 + 
					 ( m_mfr > 0.0 ? m_mfr : 0.0 );
	m_coef_1[2][3] = 0.0;

	m_coef_1[3][0] = 0.0;
	m_coef_1[3][1] = 0.0;
	m_coef_1[3][2] = 0.0;
	m_coef_1[3][3] = get_info().get_params()->get_gamma() *
					 m_alpha * m_area / m_s12 + 
					 ( m_mfr > 0.0 ? m_mfr * get_info().get_params()->get_c_m() : 0.0 );
}

void 
NIFS::c_face_int::calc_coef_2()
{
	m_coef_2[0][0] = - m_alpha * m_df * m_area / m_s12;
	m_coef_2[0][1] = 0.5 * m_unv.get_x() * m_area;
	m_coef_2[0][2] = 0.5 * m_unv.get_y() * m_area;
	m_coef_2[0][3] = 0.0;

	m_coef_2[1][0] = 0.5 * m_unv.get_x() * m_area;
	
	m_coef_2[1][1] = - get_info().get_params()->get_mu() * 
					 m_alpha * m_area / m_s12 +
					 ( m_mfr < 0.0 ? m_mfr : 0.0 );
	m_coef_2[1][2] = 0.0;
	m_coef_2[1][3] = 0.0;

	m_coef_2[2][0] = 0.5 * m_unv.get_y() * m_area;
	m_coef_2[2][1] = 0.0;
	m_coef_2[2][2] = - get_info().get_params()->get_mu() * 
					 m_alpha * m_area / m_s12 +
					 ( m_mfr < 0.0 ? m_mfr : 0.0 );
	m_coef_2[2][3] = 0.0;

	m_coef_2[3][0] = 0.0;
	m_coef_2[3][1] = 0.0;
	m_coef_2[3][2] = 0.0;
	m_coef_2[3][3] = - get_info().get_params()->get_gamma() *
					 m_alpha * m_area / m_s12 + 
					 ( m_mfr < 0.0 ? m_mfr * get_info().get_params()->get_c_m() : 0.0 );
}

void 
NIFS::c_face_int::calc_coef_c()
{
	c_cell* cell_ptr_1 = static_cast<c_cell*>( get_cell( 0 ) );
	c_cell* cell_ptr_2 = static_cast<c_cell*>( get_cell( 1 ) );

	m_coef_c[0] = 
		0.5 * m_area * ( 
		( cell_ptr_1->get_grad().get_velu() +
		  cell_ptr_2->get_grad().get_velu() ) * m_rc * m_unv.get_x() +
		( cell_ptr_1->get_grad().get_velv() +
		  cell_ptr_2->get_grad().get_velv() ) * m_rc * m_unv.get_y() 
					   ) + 
		0.5 * m_alpha * m_df * m_area * (
			( cell_ptr_1->get_grad().get_pres() +
			  cell_ptr_2->get_grad().get_pres() ) * m_uv12 );

	m_coef_c[1] = 
		0.5 * m_unv.get_x() * m_area * 
			( m_rc * cell_ptr_1->get_grad().get_pres() +
			  m_rc * cell_ptr_2->get_grad().get_pres() 
			) -
		0.5 * get_info().get_params()->get_mu() * m_area * 
			( ( cell_ptr_1->get_grad().get_velu() +
				cell_ptr_2->get_grad().get_velu() ) *
			  ( m_unv - m_uv12 * m_alpha ) 
			) -
		0.5 * get_info().get_params()->get_mu() * m_area * 
			( ( cell_ptr_1->get_grad().get_velu().get_x() +
				cell_ptr_2->get_grad().get_velu().get_x() ) *
			  m_unv.get_x() +
			  ( cell_ptr_1->get_grad().get_velv().get_x() +
				cell_ptr_2->get_grad().get_velv().get_x() ) *
			  m_unv.get_y() 
			) + 
		( m_mfr >= 0.0 ?
			m_mfr * ( m_r1 * cell_ptr_1->get_grad().get_velu() ) :
			m_mfr * ( m_r2 * cell_ptr_2->get_grad().get_velu() )
		);

	m_coef_c[2] = 
		0.5 * m_unv.get_y() * m_area * 
			( m_rc * cell_ptr_1->get_grad().get_pres() +
			  m_rc * cell_ptr_2->get_grad().get_pres() 
			) -
		0.5 * get_info().get_params()->get_mu() * m_area * 
			( ( cell_ptr_1->get_grad().get_velv() +
				cell_ptr_2->get_grad().get_velv() ) *
			  ( m_unv - m_uv12 * m_alpha ) 
			) - 
		0.5 * get_info().get_params()->get_mu() * m_area * 
			( ( cell_ptr_1->get_grad().get_velu().get_y() +
				cell_ptr_2->get_grad().get_velu().get_y() 
			  ) * m_unv.get_x() +
			  ( cell_ptr_1->get_grad().get_velv().get_y() +
				cell_ptr_2->get_grad().get_velv().get_y() 
			  ) * m_unv.get_y()
			) + 
		( m_mfr >= 0.0 ?
			m_mfr * ( m_r1 * cell_ptr_1->get_grad().get_velv() ) :
			m_mfr * ( m_r2 * cell_ptr_2->get_grad().get_velv() )
		);

	m_coef_c[3] = -0.5 * get_info().get_params()->get_gamma() * m_area *
		( ( cell_ptr_1->get_grad().get_temp() +
			cell_ptr_2->get_grad().get_temp() ) * ( m_unv - m_uv12 * m_alpha )
		) + 
		m_mfr * get_info().get_params()->get_c_m() * 
		( m_mfr >= 0.0 ?
			( m_r1 * cell_ptr_1->get_grad().get_temp() ) :
			( m_r2 * cell_ptr_2->get_grad().get_temp() )
		);
}

void 
NIFS::c_face_int::calc_mfr()
{
  m_flow = ( m_coef_1 * static_cast<c_cell*>( get_cell(0) )->get_state() +
             m_coef_2 * static_cast<c_cell*>( get_cell(1) )->get_state() +
             m_coef_c ); 
	m_mfr = m_flow( 0 ) * get_info().get_params()->get_rho();
}

void 
NIFS::c_face_int::calc_state()
{
	c_cell* cell_ptr_1 = static_cast<c_cell*>( get_cell(0) );
	c_cell* cell_ptr_2 = static_cast<c_cell*>( get_cell(1) );

	for ( int i = 0 ; i < 4 ; ++i ) {
		m_state[i] = (  cell_ptr_1->get_state()(i) + 
						cell_ptr_2->get_state()(i) + 
						( cell_ptr_1->get_grad()(i) +
						  cell_ptr_2->get_grad()(i) ) * m_rc  
					 ) / 2.0;

		m_grad[i] = ( cell_ptr_1->get_grad()(i) +
					  cell_ptr_2->get_grad()(i) 
					) / 2.0;

		m_hess[i] = ( cell_ptr_1->get_hess()(i) +
					  cell_ptr_2->get_hess()(i) 
					) / 2.0;
	} 
}

void 
NIFS::c_face_int::calc_error()
{
	c_cell* cell_ptr_1 = static_cast<c_cell*>( get_cell( 0 ) );
	c_cell* cell_ptr_2 = static_cast<c_cell*>( get_cell( 1 ) );
	
	m_error[0] = 
		0.125 * m_area * (
			( cell_ptr_1->get_hess().get_velu() +
			  cell_ptr_2->get_hess().get_velu() 
			).gen_dot_prod( m_r1 , m_r2 ) * m_unv.get_x() +
			( cell_ptr_1->get_hess().get_velv() +
			  cell_ptr_2->get_hess().get_velv() 
			).gen_dot_prod( m_r1 , m_r2 ) * m_unv.get_y() +
			( cell_ptr_1->get_hess().get_velu() +
			  cell_ptr_2->get_hess().get_velu() 
			).gen_dot_prod( m_r2 , m_r1 ) * m_unv.get_x() +
			( cell_ptr_1->get_hess().get_velv() +
			  cell_ptr_2->get_hess().get_velv() 
			).gen_dot_prod( m_r2 , m_r1 ) * m_unv.get_y() +
			( cell_ptr_1->get_hess().get_velu() +
			  cell_ptr_2->get_hess().get_velu() 
			).gen_dot_prod( m_utv * m_area , m_utv * m_area ) * m_unv.get_x() / 6.0 +
			( cell_ptr_1->get_hess().get_velv() +
			  cell_ptr_2->get_hess().get_velv() 
			).gen_dot_prod( m_utv * m_area , m_utv * m_area ) * m_unv.get_y() / 6.0 ) -
		0.25 * m_alpha * m_df * m_area / m_s12 * (
			( cell_ptr_1->get_hess().get_pres() +
			  cell_ptr_2->get_hess().get_pres() 
			).gen_dot_prod( m_r1 , m_r2 ) - 
			( cell_ptr_1->get_hess().get_pres() +
			  cell_ptr_2->get_hess().get_pres() 
			).gen_dot_prod( m_r2 , m_r1 ) ); 

	double mfr = m_mfr + get_info().get_params()->get_rho() * m_error(0);

	m_error[1] = 
		( m_mfr >= 0.0 ?
			0.5 * mfr * (cell_ptr_1->get_hess().get_velu()).gen_dot_prod( m_r1 , m_r1 ) + 
			get_info().get_params()->get_rho() * m_area / 24.0 * (
				( cell_ptr_1->get_grad().get_velu() +
				  cell_ptr_2->get_grad().get_velu() ) * m_utv * m_area * m_unv.get_x() +
				( cell_ptr_1->get_grad().get_velv() +
				  cell_ptr_2->get_grad().get_velv() ) * m_utv * m_area * m_unv.get_y() ) *
				( cell_ptr_1->get_grad().get_velu() * m_utv * m_area ) +
			mfr / 24.0 *
				( cell_ptr_1->get_hess().get_velu()).gen_dot_prod( m_utv * m_area , m_utv * m_area )
			:
			0.5 * mfr * (cell_ptr_2->get_hess().get_velu()).gen_dot_prod( m_r2 , m_r2 ) + 
			get_info().get_params()->get_rho() * m_area / 24.0 * (
				( cell_ptr_1->get_grad().get_velu() +
				  cell_ptr_2->get_grad().get_velu() ) * m_utv * m_area * m_unv.get_x() +
				( cell_ptr_1->get_grad().get_velv() +
				  cell_ptr_2->get_grad().get_velv() ) * m_utv * m_area * m_unv.get_y() ) *
				( cell_ptr_2->get_grad().get_velu() * m_utv * m_area ) +
			mfr / 24.0 *
				( cell_ptr_2->get_hess().get_velu()).gen_dot_prod( m_utv * m_area , m_utv * m_area )
		) + 
		0.125 * m_unv.get_x() * m_area * (
			( cell_ptr_1->get_hess().get_pres() +
			  cell_ptr_2->get_hess().get_pres() ).gen_dot_prod( m_utv * m_area , m_utv * m_area ) / 6.0 +
			( cell_ptr_1->get_hess().get_pres() +
			  cell_ptr_2->get_hess().get_pres() ).gen_dot_prod( m_r1 , m_r2 ) +
			( cell_ptr_1->get_hess().get_pres() +
			  cell_ptr_2->get_hess().get_pres() ).gen_dot_prod( m_r2 , m_r1 ) 
		) -   
		0.25 * get_info().get_params()->get_mu() * m_alpha * m_area / m_s12 * (
			( cell_ptr_1->get_hess().get_velu() +
			  cell_ptr_2->get_hess().get_velu() ).gen_dot_prod( m_r1 , m_r2 ) -
			( cell_ptr_1->get_hess().get_velu() +
			  cell_ptr_2->get_hess().get_velu() ).gen_dot_prod( m_r2 , m_r1 ) 
		) -   
		0.5 * get_info().get_params()->get_mu() * m_area * (
			( cell_ptr_1->get_hess().get_velu() +
			  cell_ptr_2->get_hess().get_velu() ).gen_dot_prod( m_rc , m_unv )
		) -
		0.5 * get_info().get_params()->get_mu() * m_area * (
			( ( cell_ptr_1->get_hess().get_velu()(0)(0) + 
				cell_ptr_2->get_hess().get_velu()(0)(0) ) * m_rc.get_x() +
			  ( cell_ptr_1->get_hess().get_velu()(0)(1) + 
				cell_ptr_2->get_hess().get_velu()(0)(1) ) * m_rc.get_y() 
			) * m_unv.get_x() + 
			( ( cell_ptr_1->get_hess().get_velv()(0)(0) + 
				cell_ptr_2->get_hess().get_velv()(0)(0) ) * m_rc.get_x() +
			  ( cell_ptr_1->get_hess().get_velv()(0)(1) + 
				cell_ptr_2->get_hess().get_velv()(0)(1) ) * m_rc.get_y() 
			) * m_unv.get_y()
		);

	m_error[2] = 
		( m_mfr >= 0.0 ?
			0.5 * mfr * (cell_ptr_1->get_hess().get_velv()).gen_dot_prod( m_r1 , m_r1 ) + 
			get_info().get_params()->get_rho() * m_area / 24.0 * (
				( cell_ptr_1->get_grad().get_velu() +
				  cell_ptr_2->get_grad().get_velu() ) * m_utv * m_area * m_unv.get_x() +
				( cell_ptr_1->get_grad().get_velv() +
				  cell_ptr_2->get_grad().get_velv() ) * m_utv * m_area * m_unv.get_y() ) *
				( cell_ptr_1->get_grad().get_velv() * m_utv * m_area ) +
			mfr / 24.0 *
				( cell_ptr_1->get_hess().get_velv()).gen_dot_prod( m_utv * m_area , m_utv * m_area )
			:
			0.5 * mfr * (cell_ptr_2->get_hess().get_velv()).gen_dot_prod( m_r2 , m_r2 ) + 
			get_info().get_params()->get_rho() * m_area / 24.0 * (
				( cell_ptr_1->get_grad().get_velu() +
				  cell_ptr_2->get_grad().get_velu() ) * m_utv * m_area * m_unv.get_x() +
				( cell_ptr_1->get_grad().get_velv() +
				  cell_ptr_2->get_grad().get_velv() ) * m_utv * m_area * m_unv.get_y() ) *
				( cell_ptr_2->get_grad().get_velv() * m_utv * m_area ) +
			mfr / 24.0 *
				( cell_ptr_2->get_hess().get_velv()).gen_dot_prod( m_utv * m_area , m_utv * m_area )
		) + 
		0.125 * m_unv.get_y() * m_area * (
			( cell_ptr_1->get_hess().get_pres() +
			  cell_ptr_2->get_hess().get_pres() ).gen_dot_prod( m_utv * m_area , m_utv * m_area ) / 6.0 +
			( cell_ptr_1->get_hess().get_pres() +
			  cell_ptr_2->get_hess().get_pres() ).gen_dot_prod( m_r1 , m_r2 ) +
			( cell_ptr_1->get_hess().get_pres() +
			  cell_ptr_2->get_hess().get_pres() ).gen_dot_prod( m_r2 , m_r1 ) 
		) -   
		0.25 * get_info().get_params()->get_mu() * m_alpha * m_area / m_s12 * (
			( cell_ptr_1->get_hess().get_velv() +
			  cell_ptr_2->get_hess().get_velv() ).gen_dot_prod( m_r1 , m_r2 ) -
			( cell_ptr_1->get_hess().get_velv() +
			  cell_ptr_2->get_hess().get_velv() ).gen_dot_prod( m_r2 , m_r1 ) 
		) -   
		0.5 * get_info().get_params()->get_mu() * m_area * (
			( cell_ptr_1->get_hess().get_velv() +
			  cell_ptr_2->get_hess().get_velv() ).gen_dot_prod( m_rc , m_unv )
		) -
		0.5 * get_info().get_params()->get_mu() * m_area * (
			( ( cell_ptr_1->get_hess().get_velu()(1)(0) + 
				cell_ptr_2->get_hess().get_velu()(1)(0) ) * m_rc.get_x() +
			  ( cell_ptr_1->get_hess().get_velu()(1)(1) + 
				cell_ptr_2->get_hess().get_velu()(1)(1) ) * m_rc.get_y() 
			) * m_unv.get_x() + 
			( ( cell_ptr_1->get_hess().get_velv()(1)(0) + 
				cell_ptr_2->get_hess().get_velv()(1)(0) ) * m_rc.get_x() +
			  ( cell_ptr_1->get_hess().get_velv()(1)(1) + 
				cell_ptr_2->get_hess().get_velv()(1)(1) ) * m_rc.get_y() 
			) * m_unv.get_y()
		);

	m_error[3] = 0.0;
}


