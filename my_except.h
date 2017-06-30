// Incompressible Flow Solver

// Copyright (C) 2005  Amir R. Baserinia

//------------------------------------------------------------------------------
// Created on:	18 Mar 2005
// Last update:	18 Mar 2005
//------------------------------------------------------------------------------

#ifndef CLASS_NIFS_MY_EXCEPT_H
#define CLASS_NIFS_MY_EXCEPT_H

#include <string>

namespace NIFS{
	
//--- INTERFACE ---	
	
class my_except {
	public:
		my_except( const std::string& msg ) { m_msg = msg; }
		virtual ~my_except( ) { }
		virtual const char* what() const throw() { return m_msg.c_str(); }
		virtual void set_msg( const std::string& msg ) throw() { m_msg = msg; }
		std::string ui_to_str( unsigned n );
	protected:
		std::string m_msg;
};	


class param_except : public my_except {
	public:
		param_except( const std::string& msg );
		virtual ~param_except( ) { }
		virtual const char* what() const throw();
		void set_line_num( unsigned line_num ) { m_line_num = line_num; }
		void set_file_name( std::string file_name ) { m_file_name = file_name; }
	private:
		std::string m_file_name;
		unsigned 	m_line_num;
};


}

#endif // CLASS_NIFS_MY_EXCEPT_H
