// Incompressible Flow Solver

// Copyright (C) 2005  Amir R. Baserinia

//----------------------------------------------------------------------------
// Created on:	07 Mar 2005
// Last update: 07 Mar 2005
//----------------------------------------------------------------------------

#ifndef CLASS_NIFS_PARSER_H
#define CLASS_NIFS_PARSER_H

#include <string>
#include <list>
#include <map>
#include <cmath>

namespace NIFS{
		  

class c_parser {
	
	public:
		c_parser();
		~c_parser() {}
		
		void set_param( std::string name, double value );
		void eval_postfix( const std::string& pfstr);
		void clear_param( void ) { m_param.clear(); }
		
		double get_param( std::string name );
		double get_answer( void ) { return m_answer; };
		std::string get_token( const std::string& pfstr, unsigned n );
		
    private:
        // binary operators
		double fadd( double v1, double v2 ) { return v1 + v2; }
		double fsub( double v1, double v2 ) { return v1 - v2; }
		double fmul( double v1, double v2 ) { return v1 * v2; }
		double fdiv( double v1, double v2 ) { return v1 / v2; }
		double fpow( double v1, double v2 ) { return std::pow( v1, v2); }

		// unary operators
		double fpos( double v )   { return v; }
		double fneg( double v )   { return -v; }
		double fabs( double v )   { return std::abs( v ); }
		double ffloor( double v ) { return std::floor( v ); }
		double fceil( double v )  { return std::ceil( v ); }
		double fsqrt( double v )  { return std::sqrt( v ); }
		double fln( double v )    { return std::log( v ); }
		double flog( double v )   { return std::log10( v ); }
		double fexp( double v )   { return std::exp( v ); }
		double fcos( double v )   { return std::cos( v ); }
		double fsin( double v )   { return std::sin( v ); }
		double ftan( double v )   { return std::tan( v ); }
		double facos( double v )  { return std::acos( v ); }
		double fasin( double v )  { return std::asin( v ); }
		double fatan( double v )  { return std::atan( v ); }
		double ferf( double v );
		double fhmz( double v );
		double fhmzd( double v );

		void clear( void );
		void set_postfix( std::string pfstr );
		void tokenize_postfix();
		void eval_postfix_list();		
		void convert_lowcase( std::string& str );
		void set_error( unsigned int error_code, std::string error_token );
		
		const std::string::size_type m_max_size;
		const std::string m_delims;
		const std::string m_numers;
		
		std::string m_pfstr; // postfix string
		std::list< std::string > m_pflist; // postfix list
		std::map< std::string, double > m_param; // parameters
		
		typedef double (NIFS::c_parser::*Tfunc_ptr2)( double, double );
		std::map< std::string, Tfunc_ptr2 > func_ptr2;
		
		typedef double (NIFS::c_parser::*Tfunc_ptr1)( double );
		std::map< std::string, Tfunc_ptr1 > func_ptr1;
		
		double m_answer;
};

} // namespace NIFS


#endif // CLASS_NIFS_PARSER_H
