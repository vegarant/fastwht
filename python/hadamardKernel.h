#include<complex>
/*

This header file contains an interface to the code which is used from the python
applications.  

*/

int PAL(unsigned int N, unsigned int n, unsigned int t);
int WAL(unsigned int N, unsigned int n, unsigned int t); 

void fwhtKernelSequency(int n, double *arr);
void fwhtKernelOrdinary(int n, double *arr);
void fwhtKernelPaley(int n, double *arr);




