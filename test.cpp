
#include <iostream>
#include "vector_block.h"
#include "vector_2d.h"

int main( int argc, char *argv[] )
{
	NIFS::vector_block< NIFS::c_vector_2d > vec;
	NIFS::c_vector_2d  vec1(0.0,1.0);
	NIFS::c_vector_2d  vec2(1.0,0.0);
	NIFS::c_vector_2d  vec3(1.0,0.0);
	
	double d; 
	d =  (vec3+vec1) ^ vec2;//( vec1.get_unit() ^ vec2 );
//	d = vec2^vec1; 
	
	std::cout << d << std::endl;	
}
