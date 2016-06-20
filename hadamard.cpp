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

uint32_t idx_from_ordinary_to_sequency(uint32_t a, uint32_t N);

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
uint32_t idx_from_ordinary_to_sequency(uint32_t a, uint32_t nu) {
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
    if (N < 2) {
        return;
    }

    hadamardOrdinary<T>(x, N);
    T *y = new T[N];
    uint32_t pos;

    const int dyadicPower = findMostSignificantBit(N) - 1;

    for ( uint32_t i = 0; i < N; i++ ) {
        pos = idx_from_ordinary_to_sequency(i,dyadicPower);
        y[pos] = x[i];
    }

    std::memcpy(x,y,N*sizeof(T));
    delete [] y;
}

template void hadamardSequency<>(short * x, const uint32_t N);
template void hadamardSequency<>(unsigned int * x, const uint32_t N);
template void hadamardSequency<>(int * x, const uint32_t N);
template void hadamardSequency<>(long * x, const uint32_t N);
template void hadamardSequency<>(float * x, const uint32_t N);
template void hadamardSequency<>(double * x, const uint32_t N);
template void hadamardSequency<>(long double * x, const uint32_t N);
template void hadamardSequency<>(std::complex<double>* x, const uint32_t N);

/*

Perform the fast Walsh-Hadamard transform in Paley order.

*/
template <typename T>
void hadamardPaley(T * x, const uint32_t N) {
    if (N < 2) {
        return;
    }

    T tmpElement;
    uint32_t pos;

    // Permute the vector such that all indices are changed with their
    // bit-reversed version.
    for (uint32_t i = 1; i < N; i++) {
        pos = reverseBitSequence(N, i);
        if (i < pos) { // swap elements
            tmpElement = x[pos];
            x[pos] = x[i];
            x[i] = tmpElement;
        }
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

Does the Walsh-Hadamard transform using the ordinary order.
The calculations happens in-place.

*/
template <typename T>
void hadamardOrdinary(T *x, const uint32_t N) {
    int R = findMostSignificantBit(N) - 1;
    uint32_t scale;
    uint32_t startIdx;
    uint32_t scaleDiv2;
    uint32_t startIdxPlussScaleDiv2;
    T elem1;
    T elem2;
    for (uint32_t k = 1; k <= R; k++) {
        scale = powDyadic(k);
        for (uint32_t step = 0; step < N/scale; step++) {
            startIdx = step*scale;
            scaleDiv2 = scale/2;
            startIdxPlussScaleDiv2 = startIdx + scaleDiv2;

            for (uint32_t i=0; i < scaleDiv2; i++) {
                elem1 = x[startIdx + i] + x[startIdxPlussScaleDiv2 + i];
                elem2 = x[startIdx + i] - x[startIdxPlussScaleDiv2 + i];

                x[startIdx + i] = elem1;
                x[startIdxPlussScaleDiv2 + i] = elem2;
            }
        }
    }
}

// Specialization
template void hadamardOrdinary<>(unsigned int* x, const uint32_t N);
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
uint32_t powDyadic(const uint32_t k) {
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

