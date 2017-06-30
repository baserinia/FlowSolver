// Incompressible Flow Solver

// Copyright (C) 2005  Amir R. Baserinia

//------------------------------------------------------------------------------
// Date:   28 Apr 2005
// Update: 28 Apr 2005
//------------------------------------------------------------------------------

#include "face_bnd_wall_dir_pres.h"
#include "cell.h"

//--- IMPLEMENTATION ---

NIFS::c_face_bnd_wall_dir_pres::c_face_bnd_wall_dir_pres() : 
									c_face_bnd_wall()
{
	m_weight.set_pres( 1.0 );
	m_weight.set_velu( 1.0 );
	m_weight.set_velv( 1.0 );
	m_weight.set_temp( 1.0 );
}

void 
NIFS::c_face_bnd_wall_dir_pres::calc_coef_1()
{
	m_coef_1[0][0] = m_alpha * m_df * m_area / m_s12;
	m_coef_1[0][1] = m_unv.get_x() * m_area;
	m_coef_1[0][2] = m_unv.get_y() * m_area;
	m_coef_1[0][3] = 0.0;

	m_coef_1[1][0] = 0.0;
	m_coef_1[1][1] = get_info().get_params()->get_mu() * 
					 m_area * m_alpha * sqr( m_utv.get_x() ) / m_s12;
	m_coef_1[1][2] = get_info().get_params()->get_mu() * 
					 m_area * m_alpha * m_utv.get_x() * m_utv.get_y() / m_s12;
	m_coef_1[1][3] = 0.0;

	m_coef_1[2][0] = 0.0;
	m_coef_1[2][1] = get_info().get_params()->get_mu() * 
					 m_area * m_alpha * m_utv.get_x() * m_utv.get_y() / m_s12;
	m_coef_1[2][2] = get_info().get_params()->get_mu() * 
					 m_area * m_alpha * sqr( m_utv.get_y() ) / m_s12;
	m_coef_1[2][3] = 0.0;

	m_coef_1[3][0] = 0.0;
	m_coef_1[3][1] = 0.0;
	m_coef_1[3][2] = 0.0;
	m_coef_1[3][3] = get_info().get_params()->get_gamma() * m_area * m_alpha / m_s12;
}


void 
NIFS::c_face_bnd_wall_dir_pres::calc_coef_c()
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

	m_coef_c[1] = m_mfr * m_state.get_velu() -
				  ( get_info().get_params()->get_mu() * 
					m_area * m_alpha * m_utv.get_x() / m_s12 ) * 
					( m_state.get_velu() * m_utv.get_x() +
					  m_state.get_velv() * m_utv.get_y() ) +
				  m_area * m_unv.get_x() * m_state.get_pres() -
				  get_info().get_params()->get_mu() * m_area * m_utv.get_x() *
					( cell_ptr_1->get_grad().get_velu() *
					  ( m_unv * m_utv.get_x() + 
						m_utv * m_unv.get_x() - 
						m_uv12 * m_alpha * m_utv.get_x() ) +
					  cell_ptr_1->get_grad().get_velv() *
					  ( m_unv * m_utv.get_y() + 
						m_utv * m_unv.get_y() - 
						m_uv12 * m_alpha * m_utv.get_y() ) 
					);

	m_coef_c[2] = m_mfr * m_state.get_velv() -
				  ( get_info().get_params()->get_mu() * 
					m_area * m_alpha * m_utv.get_y() / m_s12 ) * 
					( m_state.get_velu() * m_utv.get_x() +
					  m_state.get_velv() * m_utv.get_y() ) +
				  m_area * m_unv.get_y() * m_state.get_pres() -
				  get_info().get_params()->get_mu() * m_area * m_utv.get_y() *
					( cell_ptr_1->get_grad().get_velu() *
					  ( m_unv * m_utv.get_x() + 
						m_utv * m_unv.get_x() - 
						m_uv12 * m_alpha * m_utv.get_x() ) +
					  cell_ptr_1->get_grad().get_velv() *
					  ( m_unv * m_utv.get_y() + 
						m_utv * m_unv.get_y() - 
						m_uv12 * m_alpha * m_utv.get_y() )
					);
					
	m_coef_c[3] = -get_info().get_params()->get_gamma() * m_area * 
			( cell_ptr_1->get_grad().get_temp() * 
			  ( m_unv - m_uv12 * m_alpha ) +
			  ( m_state.get_temp() * m_alpha / m_s12 ) 
			) +
			m_mfr * m_state.get_temp() * 
			get_info().get_params()->get_c_m();
}


void 
NIFS::c_face_bnd_wall_dir_pres::calc_bc()
{
	//--- calculate pressure at pressure specified wall
	m_state[0] = get_info().get_bnd_cond()->get_pres( 
			get_code(), m_fmp.get_pos().get_x(), m_fmp.get_pos().get_y() );
	m_state[1] = get_info().get_bnd_cond()->get_velu( 
			get_code(), m_fmp.get_pos().get_x(), m_fmp.get_pos().get_y() );
	//--- calculate vel-v component at wall
	m_state[2] = get_info().get_bnd_cond()->get_velv( 
			get_code(), m_fmp.get_pos().get_x(), m_fmp.get_pos().get_y() );
	//--- calculate temperature at wall
	m_state[3] = get_info().get_bnd_cond()->get_temp( 
			get_code(), m_fmp.get_pos().get_x(), m_fmp.get_pos().get_y() );
}

void 
NIFS::c_face_bnd_wall_dir_pres::calc_state()
{
	c_cell* cell_ptr_1 = static_cast<c_cell*>( get_cell(0) );
	
	for ( int i = 0 ; i < 4 ; ++i ) 
		m_grad[i] =  cell_ptr_1->get_grad()(i) + 
					 cell_ptr_1->get_hess()(i) * m_r1;

	m_hess =  cell_ptr_1->get_hess();
}
