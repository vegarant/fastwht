
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
%apply (int DIM1, int DIM2, double* INPLACE_ARRAY2) {(int n, int m, double *arr)};
%apply (int DIM1, int DIM2, long* INPLACE_ARRAY2) {(int n, int m, long *arr)};
%include "hadamardKernel.h"

