// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, version 3 of the License.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program.  If not, see <http://www.gnu.org/licenses/>.
//
// Copyright 2016 Vegard Antun
//

/*

This program in the computational core of the Walsh-Hadamard transform.

*/


#include "hadamard.h"
#include "python/hadamardKernel.h"


/*

Computes the element in the sequency ordered Walsh-Hadamard matrix.

N - Dimension of the matrix N × N
n - row number i.e., Walsh-Hadamard function ψ_n
t - column i.e., the input ψ_n(t/N)

*/
int WAL(uint32_t N, uint32_t n, uint32_t t) {
    const uint32_t dyadicPower = findMostSignificantBit(N) - 1;
    uint32_t s = 0;
    uint32_t n_pr = 0x00000001;
    uint32_t t_r  = 0x00000001;
    n_pr = n_pr << (dyadicPower-1);
    uint32_t n_r = 0;
    uint32_t t_s = 0;

    for (uint32_t r = 0; r < dyadicPower; r++) {

        n_r = (((n_pr >> r) & n) >> (dyadicPower - r-1));
        t_s = (((t_r << r) & t) >> r) - (((t_r << (r+1)) & t) >> (r+1));
        s += n_r*t_s;

    }

    if ( s%2 ) { // s is odd
        return -1;
    } else {
        return 1;
    }
}

/*

Computes the element in the Paley ordered Walsh-Hadamard matrix.

N - Dimension of the matrix N × N
n - row number i.e., Walsh-Hadamard function ψ_n
t - column i.e., the input ψ_n(t/N)

*/
int PAL(uint32_t N, uint32_t n, uint32_t x) {
    const int dyadicPower = findMostSignificantBit(N) - 1;
    uint32_t s = 0;
    uint32_t ONE = 0x00000001;
    uint32_t n_j = 0;
    uint32_t x_jp1 = 0;

    for (uint32_t j = 0; j < dyadicPower; j++) {

        n_j = (n & (ONE << j)) >> j ;
        x_jp1 = (x & (ONE << (dyadicPower - j-1))) >> (dyadicPower - j-1);

        s += n_j*x_jp1;

    }

    if ( s%2 ) { // s is odd
        return -1;
    } else {
        return 1;
    }
}

/*

Finds the position of the most significant bit in the uint32_t 'a'. If
no bit is set to 1, it will return 0.

*/
unsigned int findMostSignificantBit(uint32_t a) {
    uint32_t k =  0x80000000;
    int i = sizeof(uint32_t)*8-1;
    while(i >= 0) {
        if (a & k) {
            return ffs(k);
        }
        i--;
        k = k >> 1;
    }
    return 0;
}

/*

Reverses the bit sequence of x, assuming one is using R-bits with N = 2^R.

*/
uint32_t reverseBitSequence(const uint32_t N, uint32_t x) {
    const unsigned int n = findMostSignificantBit(N) - 1;

    uint32_t out = 0;
    uint32_t ONE = 0x00000001;

    uint32_t ONE_rev = ONE << (n-1);
    for (int i = 0; i < n; i++) {
        if (x & ONE) {
            out = out | ONE_rev;
        }

        ONE = ONE << 1;
        ONE_rev = ONE_rev >> 1;
    }

    return out;
}

/*

This function is a direct copy of Wikipedia's implementation

It converts a binary number to gray code.

*/
uint32_t binaryToGrayCode( uint32_t x ) {
    return x ^ (x >> 1);
}

/*

This function is a direct copy of Wikipedia's implementation

It converts a gray coded number into a binary ordered number.

*/
uint32_t grayCodeToBinary(uint32_t x) {
    x = x ^ (x >> 16);
    x = x ^ (x >> 8);
    x = x ^ (x >> 4);
    x = x ^ (x >> 2);
    x = x ^ (x >> 1);
    return x;
}

/*

Find the new position of index 'a' in a hadamard matrix in a sequency ordered
haramard matrix. 2^nu is the total number of entries in the hadamard matrix.

*/
uint32_t idxFromOrdinaryToSequency(uint32_t a, uint32_t nu) {
    uint32_t value = 0;
    uint32_t k = 1;
    uint32_t l = 1 << (nu-1);

    bool keepOne = true;

    for ( uint32_t i = 0; i < nu; i++) {
        if ( a & k ) { // bit is one, keep it
            if (keepOne) {
                keepOne = false;
                value = value | l;
            } else {
                keepOne = true;
                value = value & ~l;
            }
        } else { // bit is zero
            if (!keepOne) {// complement it
                value = value | l; // insert a one
            } else {
                // keep the zero
            }
        }
        //printf("k: 0x%x\n", k);
        k = k << 1;
        l = l >> 1;
    }
    return value;
}

/*

Perform the fast Walsh-Hadamard transform in sequency order.

*/
template <typename T>
void hadamardSequency(T * x, const uint32_t N) {
    if (N < 2) 
    {
        return;
    }
    
    hadamardOrdinary<T>(x, N);
    T *y = new T[N];
    uint32_t pos;

    const int dyadicPower = findMostSignificantBit(N) - 1;

    for ( uint32_t i = 0; i < N; i++ ) 
    {
        pos = idxFromOrdinaryToSequency(i,dyadicPower);
        y[pos] = x[i];
    }

    std::memcpy(x,y,N*sizeof(T));
    delete [] y;
}

template void hadamardSequency<>(short * x, const uint32_t N);
template void hadamardSequency<>(int * x, const uint32_t N);
template void hadamardSequency<>(long * x, const uint32_t N);
template void hadamardSequency<>(float * x, const uint32_t N);
template void hadamardSequency<>(double * x, const uint32_t N);
template void hadamardSequency<>(long double * x, const uint32_t N);
template void hadamardSequency<>(std::complex<double>* x, const uint32_t N);

/*

Perform the fast Walsh-Hadamard transform in Paley order.

The bit-reversal permutation step of this algorithm is a copy of Øyvind Ryan's 
code found at http://folk.uio.no/oyvindry/matinf2360/code/python/fft.py

*/
template <typename T>
void hadamardPaley(T * x, const uint32_t N) {
    if (N < 2) 
    {
        return;
    }

    T tmpElem = 0;
    uint32_t j = 0;

    for (int i = 0; i < N/2; i += 2) {
        if (j > i) 
        {
            tmpElem = x[j];
            x[j] = x[i];
            x[i] = tmpElem;

            tmpElem = x[N - j - 1];
            x[N - j - 1] = x[N - i - 1];
            x[N - i - 1] = tmpElem;
        }

        tmpElem = x[i+1];
        x[i+1] = x[j + N/2];
        x[j + N/2] = tmpElem;

        uint32_t m = N/4;

        while (m >= 1 and j >= m)
        {
            j -= m;
            m /= 2;
        }

        j += m;

    }    

    hadamardOrdinary<T>(x, N);
}

template void hadamardPaley<>(short * x, const uint32_t N);
template void hadamardPaley<>(int * x, const uint32_t N);
template void hadamardPaley<>(long * x, const uint32_t N);
template void hadamardPaley<>(float * x, const uint32_t N);
template void hadamardPaley<>(double * x, const uint32_t N);
template void hadamardPaley<>(long double * x, const uint32_t N);
template void hadamardPaley<>(std::complex<double>* x, const uint32_t N);

/*

Performs the Walsh-Hadamard transform using the ordinary order.
The calculations happens in-place.

INPUT:
x - The vector one would like to transform
N - Length of the vector. N must be a power of 2.

*/
template <typename T>
void hadamardOrdinary(T *x, const uint32_t N) {

    if (N > 65537) // N > 2^16 + 1
    {
        hadamardParallel<T>(x, N);
    } 
    else 
    {
        hadamardArndt<T>(x, N);
    }
}

// Specialization
template void hadamardOrdinary<>(int* x, const uint32_t N);
template void hadamardOrdinary<>(short* x, const uint32_t N);
template void hadamardOrdinary<>(long* x, const uint32_t N);
template void hadamardOrdinary<>(float* x, const uint32_t N);
template void hadamardOrdinary<>(double* x, const uint32_t N);
template void hadamardOrdinary<>(long double* x, const uint32_t N);
template void hadamardOrdinary<>(std::complex<double>* x, const uint32_t N);

/*

Returns  2^k

*/
uint32_t powDyadic(const unsigned int k) {
    uint32_t x = 1;
    return (x<<k);
}

/*

Various kernels used for the python binding.

*/
void fwhtKernelSequency(int n, double *arr) {
     hadamardSequency<double>(arr, n);
}

void fwhtKernelOrdinary(int n, double *arr) {
     hadamardOrdinary<double>(arr, n);
}

void fwhtKernelPaley(int n, double *arr) {
     hadamardPaley<double>(arr, n);
}

int PAL_kernel(unsigned int N, unsigned int n, unsigned int x) {
    return PAL(N,n,x);
}

int WAL_kernel(unsigned int N, unsigned int n, unsigned int x) {
    return WAL(N,n,x);
}


///////////////////////////////////////////////////////////////////////////////
///          Parallel version of the ordinary Hadamard transform            ///
///////////////////////////////////////////////////////////////////////////////

/*

Parallel version of the ordinary Hadamard transform.

INPUT:
x - The vector one would like to transform
N - Length of the vector

It is possible to use more 

*/
template <typename T>
void hadamardParallel(T* x, const uint32_t N) {
    
    // Find the number of possible threads and choose the lower dyadic bound
    const unsigned int numThreadSupported = std::thread::hardware_concurrency();

    const unsigned int msb = findMostSignificantBit(numThreadSupported) - 1;
    uint32_t numThreadUsed = powDyadic(msb);

    // Compute the main part of the transform using all threads 
    std::thread threadArray[numThreadUsed-1];

    for (int i = 1; i < numThreadUsed; i++) 
    {
        threadArray[i-1] = std::thread(hadamardArndt<T>, 
                                       &x[i*(N/numThreadUsed)], 
                                       N/numThreadUsed);
    }
    
    // Let the main thread do some work as well
    hadamardArndt<T>(x, N/numThreadUsed);

    for (int i = 0; i < numThreadUsed-1; i++) 
    {
        threadArray[i].join();
    }
    
    // Gradually reduce the number of active treads
    const unsigned int Nmsb = findMostSignificantBit(N)-1;

    //uint32_t ldm = Nmsb-msb+1;
    //while (ldm != Nmsb) 
    //{

    //    numThreadUsed /= 2;

    //    for (int i = 0; i < numThreadUsed; i++) 
    //    {
    //        threadArray[i] = std::thread(hadamardArndtOneStep<T>, 
    //                                     &x[i*(N/numThreadUsed)], 
    //                                     N/numThreadUsed, ldm);
    //    }

    //    for (int i = 0; i < numThreadUsed; i++) 
    //    {
    //        threadArray[i].join();
    //    }
    //    
    //    ldm++;
    //}
    
    // For the final part of the transform we do all the work with one thread
    for (uint32_t ldm=Nmsb-msb+1; ldm<=Nmsb; ++ldm) 
    {
        const uint32_t m = (1 << ldm);         // 2^ldm
        const uint32_t mh = (m >> 1);          // m/2
        for (uint32_t r = 0; r < N; r += m) 
        {
            uint32_t t1 = r;
            uint32_t t2 = r + mh;
            for (uint32_t j = 0; j < mh; ++j, ++t1, ++t2) 
            {
                T u = x[t1];
                T v = x[t2];
                x[t1] = u + v;
                x[t2] = u - v;
            }
        }
    }
}

template void hadamardParallel<>(short * x, const uint32_t N);
template void hadamardParallel<>(int * x, const uint32_t N);
template void hadamardParallel<>(long * x, const uint32_t N);
template void hadamardParallel<>(float * x, const uint32_t N);
template void hadamardParallel<>(double * x, const uint32_t N);
template void hadamardParallel<>(long double * x, const uint32_t N);
template void hadamardParallel<>(std::complex<double>* x, const uint32_t N);



/*

Performs the Walsh-Hadamard transform using the ordinary order.
The calculations happens in-place.

This formulation of the code is found in the book "Arndt Computational" by 
Jörg Arndt, Springer 2011. 

*/
template <typename T>
void hadamardArndt(T *x, const uint32_t N) {

    const uint32_t ldn = findMostSignificantBit(N) - 1;

    for (uint32_t ldm=1; ldm<=ldn; ++ldm) 
    {
        const uint32_t m = (1 << ldm);         // 2^ldm
        const uint32_t mh = (m >> 1);          // m/2
        for (uint32_t r = 0; r < N; r += m) 
        {
            uint32_t t1 = r;
            uint32_t t2 = r + mh;
            for (uint32_t j = 0; j < mh; ++j, ++t1, ++t2) 
            {
                T u = x[t1];
                T v = x[t2];
                x[t1] = u + v;
                x[t2] = u - v;
            }
        }
    }
}

template void hadamardArndt<>(short * x, const uint32_t N);
template void hadamardArndt<>(int * x, const uint32_t N);
template void hadamardArndt<>(long * x, const uint32_t N);
template void hadamardArndt<>(float * x, const uint32_t N);
template void hadamardArndt<>(double * x, const uint32_t N);
template void hadamardArndt<>(long double * x, const uint32_t N);
template void hadamardArndt<>(std::complex<double>* x, const uint32_t N);


/*

Preform one step in the outer loop of the HadamardArndt<T> function. It is 
only used in the parallel version of the program. 

*/
template <typename T>
void hadamardArndtOneStep(T *x, const uint32_t N, const uint32_t ldm ) {

    const uint32_t ldn = findMostSignificantBit(N) - 1;

    const uint32_t m = (1 << ldm);         // 2^ldm
    const uint32_t mh = (m >> 1);          // m/2
    for (uint32_t r = 0; r < N; r += m) 
    {
        uint32_t t1 = r;
        uint32_t t2 = r + mh;
        for (uint32_t j = 0; j < mh; ++j, ++t1, ++t2) 
        {
            T u = x[t1];
            T v = x[t2];
            x[t1] = u + v;
            x[t2] = u - v;
        }
    }
}

template void hadamardArndtOneStep<>(short * x, const uint32_t N, const uint32_t ldm );
template void hadamardArndtOneStep<>(int * x, const uint32_t N, const uint32_t ldm );
template void hadamardArndtOneStep<>(long * x, const uint32_t N, const uint32_t ldm );
template void hadamardArndtOneStep<>(float * x, const uint32_t N, const uint32_t ldm );
template void hadamardArndtOneStep<>(double * x, const uint32_t N, const uint32_t ldm );
template void hadamardArndtOneStep<>(long double * x, const uint32_t N, const uint32_t ldm );
template void hadamardArndtOneStep<>(std::complex<double>* x, const uint32_t N, const uint32_t ldm );


///////////////////////////////////////////////////////////////////////////////
///                        Functions not in use                             ///
///////////////////////////////////////////////////////////////////////////////
/*

Recursive formulation provided Anders Matheson

*/
template <typename T>
void hadamardRecursive(T *x, const unsigned int N) {
    T elem1, elem2;
    if(N == 2)
    {
        elem1 = x[0] + x[1];
        elem2 = x[0] - x[1];

        x[0] = elem1;
        x[1] = elem2;
        return;
    }

    hadamardRecursive(x, N/2);
    hadamardRecursive(x + N/2, N/2);

    for (unsigned i=0; i < N/2; i++) 
    {
        elem1 = x[i] + x[i + N/2];
        elem2 = x[i] - x[i + N/2];

        x[i] = elem1;
        x[i + N/2] = elem2;
    }
}

// Specialization
template void hadamardRecursive<>(int* x, const unsigned int N);
template void hadamardRecursive<>(short* x, const unsigned int N);
template void hadamardRecursive<>(long* x, const unsigned int N);
template void hadamardRecursive<>(float* x, const unsigned int N);
template void hadamardRecursive<>(double* x, const unsigned int N);
template void hadamardRecursive<>(long double* x, const unsigned int N);
template void hadamardRecursive<>(std::complex<double>* x, const unsigned int N);


/*

Iterative formulation provided by Anders Matheson

*/
template <typename T>
void hadamardDepthFirst(T *x, const unsigned int N) {
    T elem1, elem2;

    unsigned stack = 0;
    unsigned lN = 1;
    int offset = 0;

    do
    {
        // Recurse down the stack
        stack *= lN/2;
        lN = 2;

        // Do inner case
        elem1 = x[offset] + x[offset + 1];
        elem2 = x[offset] - x[offset + 1];

        x[offset] = elem1;
        x[offset + 1] = elem2;

        // Reduce back up again
        while(stack & 1) {
            offset -= lN;

            for (unsigned i=offset; i != offset + lN; i++) {
                elem1 = x[i] + x[i + lN];
                elem2 = x[i] - x[i + lN];

                x[i] = elem1;
                x[i + lN] = elem2;
            }

            stack /= 2;
            lN *= 2;
        }

        stack |= 1;
        offset += lN;
    }
    while(lN != N);
}

template void hadamardDepthFirst<>(int* x, const unsigned int N);
template void hadamardDepthFirst<>(short* x, const unsigned int N);
template void hadamardDepthFirst<>(long* x, const unsigned int N);
template void hadamardDepthFirst<>(float* x, const unsigned int N);
template void hadamardDepthFirst<>(double* x, const unsigned int N);
template void hadamardDepthFirst<>(long double* x, const unsigned int N);
template void hadamardDepthFirst<>(std::complex<double>* x, const unsigned int N);






