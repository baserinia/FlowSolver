// Incompressible Flow Solver

// Copyright (C) 2005  Amir R. Baserinia

//------------------------------------------------------------------------------
// Created on:	17 Sep 2004
// Last update:	01 Feb 2006
//------------------------------------------------------------------------------

#ifndef CLASS_NIFS_C_MESH_MOVE_H
#define CLASS_NIFS_C_MESH_MOVE_H

#include "mesh.h"
#include "face.h"

namespace NIFS{    

//--- INTERFACE ---

class c_mesh_move : public c_mesh  
{
public:
	c_mesh_move() {}
	virtual ~c_mesh_move() {}
	void move_vertex();
	int swap_face( c_face* face );
	void swap();
	void adapt();
private:
	double pi() { return 3.1415926535; }
	double sqr( double x ) { return x*x; }
};

} // namespace NBCFD

#endif // CLASS_NIFS_C_MESH_MOVE_H
