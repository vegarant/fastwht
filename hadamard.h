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

This header file is used to interface with any C++ implementation

For the file interfacing with python see python/hadamard.h

*/

#ifndef FASTWHT_HADAMARD_H_
#define FASTWHT_HADAMARD_H_

#include <complex>

// FXT
#include "fxt/bits/bit2pow.h" // ld
#include "fxt/walsh/walshseq.h" // walsh_seq2
#include "fxt/perm/revbinpermute.h" // revbin_permute
#include "fxt/walsh/walshwak.h" // walsh_wak
#include "fxt/walsh/walshwal.h" // walsh_wal

enum HadamardOrder {ORDINARY, PALEY, SEQUENCY};

template <typename T>
void hadamardTransform(T* x, const unsigned long N, const HadamardOrder order);

template<typename T> 
void hadamardTransform2d(T* x, const unsigned long M , 
                               const unsigned long N, 
                               const HadamardOrder order); 

template<typename T> 
void hadamardTransform2dColumn(T* x, const unsigned long M , 
                                     const unsigned long N, 
                                     const HadamardOrder order);
int WAL(unsigned long N, unsigned long n, unsigned long t);
int PAL(unsigned long N, unsigned long n, unsigned long t);

template <typename T>
void hadamardOrdinary(T *x, const unsigned long N);

template <typename T>
void hadamardSequency(T * x, const unsigned long N);

template <typename T>
void hadamardPaley(T * x, const unsigned long N);

/*

Returns  2^k

*/
inline unsigned long powDyadic(const unsigned long k) {
    return (1UL<<k);
}

#endif

