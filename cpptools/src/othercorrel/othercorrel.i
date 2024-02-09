%module othercorrel
%include "std_vector.i"
%template(DoubleVector) std::vector<double>;
%{
	#define SWIG_FILE_WITH_INIT
 	#include <fastjet/PseudoJet.hh>
	#define SWIG
	#include "othercorrel.hh"
%}

%include "othercorrel.hh" 
