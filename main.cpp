//	Copyright (C) 2006  Amir R. Baserinia
//
// 	This file is part of SMAIF (Simple Mesh Adaptor for Incompressible Flow)
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
//
//----------------------------------------------------------------------------
// Created on:	15 Mar 2005
// Last update:	20 Jun 2005
//-----------------------------------------------------------------------------

#include <iostream>
#include <exception>
#include <stdexcept>
#include <cstdlib>
#include "my_except.h"
#include "control.h"

//---( main )---
// the main function of the code
int 
main( int argc , char *argv[] )
{
	NIFS::c_control control;

	try {	
		control.initialize( argc, argv ); // initilize the solver
		control.solve();                  // solve the problem
//		control.shutdown();               // shutdown the process
	}
	catch ( const NIFS::my_except& e ) {
		std::cerr << "Error: " << e.what() << std::endl;
		return EXIT_FAILURE;
	}
	catch ( const std::exception& e ) {
		std::cerr << "System Error: " << e.what() << std::endl;
		return EXIT_FAILURE;
	}
	catch ( ... ) {
		std::cerr << "Unknown Error." << std::endl;
		return EXIT_FAILURE;
	}

	std::cout << "Execution terminated normally." << std::endl;

  return EXIT_SUCCESS;
}
