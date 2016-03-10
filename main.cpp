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

#include "python/hadamardKernel.h" // Python interface
#include "hadamard.h"


int main(int argc, char *argv[]) {
     
    unsigned int N = 4;
    int x[N];
    int *y;
    std::memset(x,0, N*sizeof(int));
    x[2] = 1;
    
    hadamardSequency<int>(x, N);
    
    for (int i = 0; i < N; i++) {
        std::cout << x[i] << std::endl;
    }
    
    
    
}









