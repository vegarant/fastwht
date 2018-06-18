#include "../hadamard.h"

inline bool isPowOf2(ulong x);
inline void zeroOut(double * x, unsigned long N);
inline unsigned long upper_power_of_two(unsigned long v);


/*
The computational core of wal and had_max_idx.

For an N × N Hadamard matrix, the values of n_long and k_long represent the 
row and column indices respectively. This function insert the right values 
matrix values into the pointer xd.


:param n_long: Array of frequencies
:param n_size: Size of n_long 
:param k_size: Array of indices 
:param k_size: Size of k_long 
:param xd: Pointer to Matlab output array
:param N: Max frequency  

:return: Nothing

*/
void wal_comp_core(unsigned long *n_long, unsigned long n_size, 
                     unsigned long *k_long, unsigned long k_size,
                     double * xd, unsigned long N, HadamardOrder order=SEQUENCY) 
{
    if (k_size == 1 && n_size == 1) {
        
        int val = WAL( N, n_long[0], k_long[0] );
        xd[0] = (double) val;
        
        //std::cout << "n_size: " << n_size << ", k_size: " << k_size << std::endl;
        //std::cout << "n_long[0]: " << n_long[0] << std::endl;
        //std::cout << "k_long[0]: " << k_long[0] << std::endl;

        //std::cout << "n_double[0]: " << n_double[0] << std::endl;
        //std::cout << "k_double[0]: " << k_double[0] << std::endl;

    } else if (k_size == 1) {

        double * x  = new double[N];
        zeroOut(x,N);

        x[k_long[0]] = 1;
        hadamardTransform<double>(x, N, SEQUENCY);

        for (unsigned long i = 0; i < n_size; i++) {
            xd[i] = x[n_long[i]];
        }

        delete [] x;    
        
    } else if (n_size == 1) {
        
        double * x  = new double[N];
        zeroOut(x,N);
        
        x[n_long[0]] = 1;
        hadamardTransform<double>(x, N, SEQUENCY);

        for (unsigned long i = 0; i < k_size; i++) {
            xd[i] = x[k_long[i]];
        }

        delete [] x;    
        
    } else { // n_size > 1 and k_size > 1

        double * x  = new double[N];
        zeroOut(x,N);
        
        if (n_size > k_size) {
            
            for (long j = 0; j < k_size; j++ ) {
                x[k_long[j]] = 1;
                hadamardTransform<double>(x, N, SEQUENCY);
                for (unsigned long i = 0; i < n_size; i++) {
                    xd[j*n_size + i] = x[n_long[i]];
                }
                zeroOut(x,N);
            }
            
        } else { // k_size > n_size
            
            for (long j = 0; j < n_size; j++ ) {
                x[n_long[j]] = 1;
                hadamardTransform<double>(x, N, SEQUENCY);
                for (unsigned long i = 0; i < k_size; i++) {
                    xd[i*n_size + j] = x[k_long[i]];
                }
                zeroOut(x,N);
            }
        
        }

        delete [] x;    

    }

}

/*
   Is `x` on the form 2^nu;
*/
inline bool isPowOf2(ulong x)
{
    return  !(x & (x-1));
}

inline void zeroOut(double * x, unsigned long N) 
{
    for (unsigned long i = 0; i < N; i++) {
        x[i] = 0;
    }
}

inline unsigned long upper_power_of_two(unsigned long v) 
{
    v--;
    v |= v >> 1;
    v |= v >> 2;
    v |= v >> 4;
    v |= v >> 8;
    v |= v >> 16;
    v |= v >> 32;
    v++;
    
    v += (v == 0); // if v == 0 return 1
    
    return v;
}
