// Incompressible Flow Solver

// Copyright (C) 2005  Amir R. Baserinia

//----------------------------------------------------------------------------
// Date:   17 Feb 2005
// Update: 17 Feb 2005
//-----------------------------------------------------------------------------

#ifndef CLASS_NIFS_C_FACE_BND_WALL_H
#define CLASS_NIFS_C_FACE_BND_WALL_H

#include "face_bnd.h"

namespace NIFS{    
	
//--- INTERFACE ---

class c_face_bnd_wall : public c_face_bnd {
	public:
		c_face_bnd_wall() : c_face_bnd() {}
		virtual ~c_face_bnd_wall() {}
		virtual void calc_error()  { m_error = 0.0; }
	
	protected:
		virtual void calc_coef_1() {}
		virtual void calc_coef_c() {}
		virtual void calc_bc() {} 
		virtual void calc_state() {}
		virtual void calc_mfr();
		double sqr( double r ) { return r * r; }
};

} //--- namespace NIFS ---

#endif //--- CLASS_NIFS_C_FACE_BND_WALL_H

