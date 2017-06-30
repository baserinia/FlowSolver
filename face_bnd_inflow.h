// Incompressible Flow Solver

// Copyright (C) 2005  Amir R. Baserinia

//------------------------------------------------------------------------------
// Date:   17 Feb 2005
// Update: 17 Feb 2005
//------------------------------------------------------------------------------

#ifndef CLASS_NIFS_C_FACE_BND_INFLOW_H
#define CLASS_NIFS_C_FACE_BND_INFLOW_H

#include "face_bnd.h"

namespace NIFS{    
	
//--- INTERFACE ---

class c_face_bnd_inflow : public c_face_bnd {
	public:
		c_face_bnd_inflow(); 
		virtual ~c_face_bnd_inflow() {}

	protected:
		virtual void calc_coef_1();
		virtual void calc_coef_c();
		virtual void calc_bc(); 
		virtual void calc_state();
		virtual void calc_mfr();
};

} //--- namespace NIFS

#endif //--- CLASS_NIFS_C_FACE_BND_INFLOW_H

