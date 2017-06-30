// Incompressible Flow Solver

// Copyright (C) 2005  Amir R. Baserinia

//------------------------------------------------------------------------------
// Date:   17 Feb 2005
// Update: 05 Apr 2005
//------------------------------------------------------------------------------

#ifndef CLASS_NIFS_C_FACE_BND_WALL_DIR_H
#define CLASS_NIFS_C_FACE_BND_WALL_DIR_H

#include "face_bnd_wall.h"

namespace NIFS{    

//--- INTERFACE ---

class c_face_bnd_wall_dir : public c_face_bnd_wall {
	public:
		c_face_bnd_wall_dir();
		virtual ~c_face_bnd_wall_dir() {}
		virtual void calc_error();		
	protected:
		virtual void calc_coef_1();
		virtual void calc_coef_c();
		virtual void calc_bc();
		virtual void calc_state();
};

} //--- namespace NIFS

#endif //--- CLASS_NIFS_C_FACE_BND_WALL_DIR_H

