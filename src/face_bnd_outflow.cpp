// Incompressible Flow Solver

// Copyright (C) 2005  Amir R. Baserinia

//------------------------------------------------------------------------------
// Date:   14 Apr 2004
// Update: 12 Apr 2005
//------------------------------------------------------------------------------

#include "face_bnd_outflow.h"
#include "cell.h"
#include <iostream>

//--- IMPLEMENTATION ---

NIFS::c_face_bnd_outflow::c_face_bnd_outflow() : c_face_bnd() 
{
	m_weight.set_pres( 1.0 );
	m_weight.set_velu( 0.0 );
	m_weight.set_velv( 0.0 );
	m_weight.set_temp( 0.0 );
}

void 
NIFS::c_face_bnd_outflow::calc_coef_1()
{
	m_coef_1[0][0] = m_alpha * m_df * m_area / m_s12;
	m_coef_1[0][1] = m_unv.get_x() * m_area;
	m_coef_1[0][2] = m_unv.get_y() * m_area;
	m_coef_1[0][3] = 0.0;

	m_coef_1[1][0] = 0.0;
	m_coef_1[1][1] = m_mfr;
	m_coef_1[1][2] = 0.0;
	m_coef_1[1][3] = 0.0;

	m_coef_1[2][0] = 0.0;
	m_coef_1[2][1] = 0.0;
	m_coef_1[2][2] = m_mfr;
	m_coef_1[2][3] = 0.0;

	m_coef_1[3][0] = 0.0;
	m_coef_1[3][1] = 0.0;
	m_coef_1[3][2] = 0.0;
	m_coef_1[3][3] = m_mfr * get_info().get_params()->get_c_m();
}


void 
NIFS::c_face_bnd_outflow::calc_coef_c()
{
	c_cell* cell_ptr_1 = static_cast<c_cell*>( get_cell(0) );

	m_coef_c[0] = m_alpha * m_df * m_area * 
				  ( ( m_uv12 * cell_ptr_1->get_grad().get_pres() ) -
					m_state.get_pres() / m_s12 ) + 
				  m_area * 
				  ( m_unv.get_x() *
				  	( m_r1 * cell_ptr_1->get_grad().get_velu() ) +
				  	m_unv.get_y() *
				  	( m_r1 * cell_ptr_1->get_grad().get_velv() ) );
/*
	// with outflow normal stress
	m_coef_c[1] = m_mfr * ( m_r1 * cell_ptr_1->get_grad().get_velu() ) +
				  m_unv.get_x() * m_area * m_state.get_pres() - 
				  get_info().get_params()->get_mu() * m_area *
					( m_unv * cell_ptr_1->get_grad().get_velu() ) - 
				  get_info().get_params()->get_mu() * m_area *
					( m_unv.get_x() * cell_ptr_1->get_grad().get_velu().get_x() +
					  m_unv.get_y() * cell_ptr_1->get_grad().get_velv().get_x() 
					);

	m_coef_c[2] = m_mfr * ( m_r1 * cell_ptr_1->get_grad().get_velv() ) +
				  m_unv.get_y() * m_area * m_state.get_pres() - 
				  get_info().get_params()->get_mu() * m_area *
					( m_unv * cell_ptr_1->get_grad().get_velv() ) - 
				  get_info().get_params()->get_mu() * m_area *
					( m_unv.get_x() * cell_ptr_1->get_grad().get_velu().get_y() +
					  m_unv.get_y() * cell_ptr_1->get_grad().get_velv().get_y() 
					);
*/

	// without outflow normal stress
	m_coef_c[1] = m_mfr * ( m_r1 * cell_ptr_1->get_grad().get_velu() ) +
				  m_unv.get_x() * m_area * m_state.get_pres() - 
				  get_info().get_params()->get_mu() * m_area * m_utv.get_x() *
				  (  cell_ptr_1->get_grad().get_velu() *
						( m_utv * m_unv.get_x() + m_unv * m_utv.get_x() ) +
					 cell_ptr_1->get_grad().get_velv() *
						( m_utv * m_unv.get_y() + m_unv * m_utv.get_y() )
				  );  

	m_coef_c[2] = m_mfr * ( m_r1 * cell_ptr_1->get_grad().get_velv() ) +
				  m_unv.get_y() * m_area * m_state.get_pres() - 
				  get_info().get_params()->get_mu() * m_area * m_utv.get_y() *
				  (  cell_ptr_1->get_grad().get_velu() *
						( m_utv * m_unv.get_x() + m_unv * m_utv.get_x() ) +
					 cell_ptr_1->get_grad().get_velv() *
						( m_utv * m_unv.get_y() + m_unv * m_utv.get_y() )
				  );  

	m_coef_c[3] = m_mfr * get_info().get_params()->get_c_m() *
		( cell_ptr_1->get_grad().get_temp() * m_r1 )
		- get_info().get_params()->get_gamma() * m_area * 
		( cell_ptr_1->get_grad().get_temp() * m_unv );

}

void 
NIFS::c_face_bnd_outflow::calc_bc()
{
	//--- calculate pressure at inlet
	m_state[0] = get_info().get_bnd_cond()->get_pres( 
			get_code(), m_fmp.get_pos().get_x(), m_fmp.get_pos().get_y() );
}

void 
NIFS::c_face_bnd_outflow::calc_state()
{
	c_cell* cell_ptr_1 = static_cast<c_cell*>( get_cell(0) );
	
	m_state[1] = cell_ptr_1->get_state().get_velu() + 
				 cell_ptr_1->get_grad().get_velu() * m_r1 + 
				 ( cell_ptr_1->get_hess().get_velu() * m_r1 ) * m_r1 * 0.5;
	m_state[2] = cell_ptr_1->get_state().get_velv() + 
				 cell_ptr_1->get_grad().get_velv() * m_r1 +
				 ( cell_ptr_1->get_hess().get_velv() * m_r1 ) * m_r1 * 0.5;
	m_state[3] = cell_ptr_1->get_state().get_temp() + 
				 cell_ptr_1->get_grad().get_temp() * m_r1 +
				 ( cell_ptr_1->get_hess().get_temp() * m_r1 ) * m_r1 * 0.5;
				 
	for ( int i = 0 ; i < 4 ; ++i ) 
		m_grad[i] =  cell_ptr_1->get_grad()(i) + 
					 cell_ptr_1->get_hess()(i) * m_r1;
	m_hess =  cell_ptr_1->get_hess();
}

void 
NIFS::c_face_bnd_outflow::calc_mfr()
{
	m_flow = ( m_coef_1 * static_cast<c_cell*>( get_cell(0) )->get_state() +
			   m_coef_c ); 
	m_mfr = m_flow( 0 ) * get_info().get_params()->get_rho();

	if ( m_mfr < 0.0 ) {
		m_mfr = 0.0;
	}
}
