
%module hadamardKernel
%{
  #define SWIG_FILE_WITH_INIT
  #include "hadamardKernel.h"
%}

%include "numpy.i"
%init %{
    import_array();
%}

%apply (int DIM1, double* INPLACE_ARRAY1) {(int n, double *arr)};
%include "hadamardKernel.h"


