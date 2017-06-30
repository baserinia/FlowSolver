// Incompressible Flow Solver

// Copyright (C) 2005  Amir R. Baserinia

//----------------------------------------------------------------------------
// Date:   25 Nov 2005
// Update: 25 Nov 2005
//----------------------------------------------------------------------------

#include "bfstream.h"

//--- IMPLEMENTATION ---

//---( generic output operator )-----------------------------------------------
NIFS::bofstream&
NIFS::bofstream::operator<<( double obj )
{
  write_bytes( &obj , sizeof( double ) );
  return *this;
}

//---( generic input operator )------------------------------------------------
NIFS::bifstream&
operator>>( NIFS::bifstream& bifs , double& obj )
{
  bifs.read_bytes( &obj , sizeof( obj) );
  return bifs;
}


//---( binary write to file )------------------------------------------------
void
NIFS::bofstream::write_bytes( const void* p , int length )
{
  if ( ( !p ) || ( length <= 0 ) ) 
	return;
  write( static_cast< const char* >( p ) , length );
}


//---( binary read from file )------------------------------------------------
void
NIFS::bifstream::read_bytes( void* p , int length )
{
  if ( ( !p ) || ( length <= 0 ) ) 
	return;
  read( static_cast< char* >( p ) , length );
}



