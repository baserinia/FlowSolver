//	Copyright (C) 2006  Amir R. Baserinia
//
// 	This file is part of SMAIF (Simple Mesh Adaptor for Incompressible Flows)
//
//	SMAIF is free software; you can redistribute it and/or modify
//	it under the terms of the GNU General Public License as published by
//	the Free Software Foundation; either version 2 of the License, or
//	(at your option) any later version.
//
//	SMAIF is distributed in the hope that it will be useful,
//	but WITHOUT ANY WARRANTY; without even the implied warranty of
//	MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//	GNU General Public License for more details.
//
//	You should have received a copy of the GNU General Public License
//	along with SMAIF; if not, write to the Free Software
//	Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA

//------------------------------------------------------------------------------
// Created on:	07 Feb 2006
// Last update:	07 Feb 2006
//------------------------------------------------------------------------------

#ifndef CLASS_NIFS_C_VERTEX_INT_H
#define CLASS_NIFS_C_VERTEX_INT_H

#include "vertex.h"

namespace NIFS{    

//--- INTERFACE ---

class c_vertex_int : public c_vertex
{
	public:
		c_vertex_int() { }
		virtual ~c_vertex_int() { }
		virtual int get_type() { return 0; } // internal vertex
};

} //--- namespace NIFS ---

#endif //--- CLASS_NIFS_C_VERTEX_INT_H ---

