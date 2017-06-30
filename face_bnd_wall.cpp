// Incompressible Flow Solver

// Copyright (C) 2005  Amir R. Baserinia

//------------------------------------------------------------------------------
// Date:   12 Apr 2005
// Update: 12 Apr 2005
//------------------------------------------------------------------------------

#include "face_bnd_wall.h"
#include "cell.h"


//--- IMPLEMENTATION ---

void 
NIFS::c_face_bnd_wall::calc_mfr()
{
	m_flow = ( m_coef_1 * static_cast<c_cell*>( get_cell(0) )->get_state() +
			   m_coef_c ); 
	m_mfr = m_flow( 0 ) * get_info().get_params()->get_rho();
}
