// Incompressible Flow lin_sys

// Copyright (C) 2005  Amir R. Baserinia

//------------------------------------------------------------------------------
// Created on:	16 Apr 2004
// Last update:	12 Apr 2005
//------------------------------------------------------------------------------

#ifndef CLASS_C_elem_H
#define CLASS_C_elem_H

namespace NIFS{    

// Application Class
class c_elem {
	public:
		c_elem() {}
		~c_elem() {} 
		bool operator<(const c_elem& elem) const;
		bool operator>(const c_elem& elem) const;
		void set_val(double dVal) { m_val = dVal; }
		void set_row( unsigned row ) { m_row = row; }
		void set_col( unsigned col ) { m_col = col; }
		double get_val() { return m_val; }
		unsigned get_row() { return m_row; }
		unsigned get_col() { return m_col; }
	private:
		double	m_val;	// value
		unsigned int m_row;	// row index
		unsigned int m_col;	// column index
};

inline bool 
c_elem::operator<( const c_elem& elem ) const
{
	if ( m_col < elem.m_col )
		return true;
	else if ( m_col == elem.m_col )
		return ( m_row < elem.m_row );
	else
		return false;
} 

inline bool 
c_elem::operator>(const c_elem& elem) const
{
	if ( m_col > elem.m_col )
		return true;
	else if ( m_col == elem.m_col )
		return ( m_row > elem.m_row );
	else
		return false;
} 

} // namespace NBCFD

#endif // CLASS_C_elem_H	
