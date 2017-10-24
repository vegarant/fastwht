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

#ifndef FASTWHT_PYTHON_HADAMARDKERNEL_H_
#define FASTWHT_PYTHON_HADAMARDKERNEL_H_

/*

This header file contains an interface to the code which is used from the python
applications.

*/

int PAL_kernel(unsigned int N, unsigned int n, unsigned int t);
int WAL_kernel(unsigned int N, unsigned int n, unsigned int t);

void fwhtKernelSequency(int n, double *arr);
void fwhtKernelOrdinary(int n, double *arr);
void fwhtKernelPaley(int n, double *arr);

#endif
