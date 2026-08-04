%module ecorrel
%include "std_vector.i"
%template(DoubleVector) std::vector<double>;
%template(IntVector) std::vector<int>;
%{
	#define SWIG_FILE_WITH_INIT
 	#include <fastjet/PseudoJet.hh>
	#define SWIG
	#include "ecorrel.hh"
%}

%include "ecorrel.hh" 
