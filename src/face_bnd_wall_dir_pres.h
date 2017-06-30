// Incompressible Flow Solver

// Copyright (C) 2005  Amir R. Baserinia

//------------------------------------------------------------------------------
// Date:   28 Apr 2005
// Update: 28 Apr 2005
//------------------------------------------------------------------------------

#ifndef CLASS_NIFS_C_FACE_BND_WALL_DIR_PRES_H
#define CLASS_NIFS_C_FACE_BND_WALL_DIR_PRES_H

#include "face_bnd_wall.h"

namespace NIFS{    

//--- INTERFACE ---

class c_face_bnd_wall_dir_pres : public c_face_bnd_wall {
	public:
		c_face_bnd_wall_dir_pres();
		virtual ~c_face_bnd_wall_dir_pres() {}

	protected:
		virtual void calc_coef_1();
		virtual void calc_coef_c();
		virtual void calc_bc();
		virtual void calc_state();
};

} //--- namespace NIFS

#endif //--- CLASS_NIFS_C_FACE_BND_WALL_DIR_PRES_H

