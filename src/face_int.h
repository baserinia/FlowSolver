// Incompressible Flow Solver

// Copyright (C) 2005  Amir R. Baserinia

//----------------------------------------------------------------------------
// Date:   10 Apr 2004
// Update: 31 Mar 2005
//-----------------------------------------------------------------------------

#ifndef CLASS_NIFS_C_FACE_INT_H
#define CLASS_NIFS_C_FACE_INT_H

#include "face.h"

namespace NIFS{    

//--- INTERFACE ---

class c_face_int : public c_face {
	public:
		c_face_int() : c_face() {}
		virtual ~c_face_int() {}
		virtual unsigned  get_type() { return 0; }
		// project is not used for internal faces
		virtual double project( const c_vector_2d& pnt, c_vector_2d& prj ) { return 0.0; } 
		
	protected:
		virtual void calc_r2();
		virtual void calc_rc();
		virtual void calc_uv12();
		virtual void calc_s12();
		virtual void calc_coef_1();
		virtual void calc_coef_2();
		virtual void calc_coef_c();        
		virtual void calc_mfr();
		virtual void calc_df();
		virtual void calc_state();
		virtual void calc_error();
};

} //--- namespace NIFS

#endif //--- CLASS_NIFS_C_FACE_INT_H

