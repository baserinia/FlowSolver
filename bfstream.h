// Incompressible Flow Solver

// Copyright (C) 2005  Amir R. Baserinia

//----------------------------------------------------------------------------
// Date:   25 Nov 2005
// Update: 25 Nov 2005
//----------------------------------------------------------------------------

#ifndef CLASS_NIFS_BFSTREAM_H
#define CLASS_NIFS_BFSTREAM_H

#include <iostream>
#include <fstream>

namespace NIFS{    

//--- INTERFACE ---
	
class bofstream : public std::ofstream {
 public:
  bofstream( ) : std::ofstream( ) { }
  bofstream( const char* file_name )
	: std::ofstream( file_name , std::ios::out | std::ios::binary ) { }
  void write_bytes( const void* , int );
  //  friend bofstream& operator<<( bofstream& , double );
  bofstream& operator<<( double );
};

	
class bifstream: public std::ifstream {
 public:
  bifstream( const char* file_name )
	: std::ifstream( file_name , std::ios::in | std::ios::binary ) { }
  void read_bytes( void* , int );
  friend bifstream& operator>>( bifstream& , double& );
};

} //--- namespace NIFS ---

#endif // CLASS_NIFS_BFSTREAM_H

