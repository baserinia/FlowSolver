// Incompressible Flow Solver

// Copyright (C) 2005  Amir R. Baserinia

//------------------------------------------------------------------------------
// Date:   05 Apr 2005
// Update: 05 Apr 2005
//------------------------------------------------------------------------------

#include "face_bnd_sym.h"
#include "cell.h"

//--- IMPLEMENTATION ---


NIFS::c_face_bnd_sym::c_face_bnd_sym() : c_face_bnd() 
{
	m_weight.set_pres( 0.0 );
	m_weight.set_velu( 0.0 );
	m_weight.set_velv( 0.0 );
	m_weight.set_temp( 0.0 );
}

void 
NIFS::c_face_bnd_sym::calc_coef_1()
{
	m_coef_1[0][0] = 0.0;
	m_coef_1[0][1] = 0.0;
	m_coef_1[0][2] = 0.0;
	m_coef_1[0][3] = 0.0;

	m_coef_1[1][0] = m_unv.get_x() * m_area;
	m_coef_1[1][1] = 2.0 * get_info().get_params()->get_mu() * 
					 m_area * m_alpha * sqr( m_unv.get_x() ) / m_s12;
	m_coef_1[1][2] = 2.0 * get_info().get_params()->get_mu() * 
					 m_area * m_alpha * m_unv.get_x() * m_unv.get_y() / m_s12;
	m_coef_1[1][3] = 0.0;

	m_coef_1[2][0] = m_unv.get_y() * m_area;
	m_coef_1[2][1] = 2.0 * get_info().get_params()->get_mu() * 
					 m_area * m_alpha * m_unv.get_x() * m_unv.get_y() / m_s12;
	m_coef_1[2][2] = 2.0 * get_info().get_params()->get_mu() * 
					 m_area * m_alpha * sqr( m_unv.get_y() ) / m_s12;
	m_coef_1[2][3] = 0.0;

	m_coef_1[3][0] = 0.0;
	m_coef_1[3][1] = 0.0;
	m_coef_1[3][2] = 0.0;
	m_coef_1[3][3] = 0.0;
}


void 
NIFS::c_face_bnd_sym::calc_coef_c()
{
	c_cell* cell_ptr_1 = static_cast<c_cell*>( get_cell(0) );
	
	m_coef_c[0] = 0.0;

	m_coef_c[1] = m_area * m_unv.get_x() * 
					( m_r1 * cell_ptr_1->get_grad().get_pres() ) -
				  2.0 * get_info().get_params()->get_mu() * m_area * m_unv.get_x() *
					( ( cell_ptr_1->get_grad().get_velu() *
						( m_unv - m_uv12 * m_alpha ) ) * m_unv.get_x() +
					  ( cell_ptr_1->get_grad().get_velv() *
						( m_unv - m_uv12 * m_alpha ) ) * m_unv.get_y() 
					);

	m_coef_c[2] = m_area * m_unv.get_y() * 
					( m_r1 * cell_ptr_1->get_grad().get_pres() ) -
				  2.0 * get_info().get_params()->get_mu() * m_area * m_unv.get_y() *
					( ( cell_ptr_1->get_grad().get_velu() *
						( m_unv - m_uv12 * m_alpha ) ) * m_unv.get_x() +
					  ( cell_ptr_1->get_grad().get_velv() *
						( m_unv - m_uv12 * m_alpha ) ) * m_unv.get_y() 
					);

	m_coef_c[3] = 0.0;
}

void 
NIFS::c_face_bnd_sym::calc_mfr()
{
	m_mfr = 0.0;
}

void 
NIFS::c_face_bnd_sym::calc_state()
{
	c_cell* cell_ptr_1 = static_cast<c_cell*>( get_cell(0) );
	
	for ( int i = 0 ; i < 4 ; ++i ) {
		m_state[i] = cell_ptr_1->get_state()(i) + 
				 cell_ptr_1->get_grad()(i) * m_r1; // +
//				 ( cell_ptr_1->get_hess()(i) * m_r1 ) * m_r1 * 0.5;
		m_grad[i] =  cell_ptr_1->get_grad()(i) + 
					 cell_ptr_1->get_hess()(i) * m_r1;
	}
	m_hess =  cell_ptr_1->get_hess();
}
