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
#include "timer.h"

int main(int argc, char *argv[]) {

    const unsigned long nu = 26;
    const unsigned long N = powDyadic(nu);
    const unsigned long M = powDyadic(nu);
    const unsigned long MN = M*N;
    int reps = 5;
    
    timeit("Sequency      ", reps, N, hadamardSequency<double>);
    timeit("Paley         ", reps, N, hadamardPaley<double>);
    timeit("Ordinary      ", reps, N, hadamardOrdinary<double>);

    //double *x = new double[MN];
    //double *y = new double[MN];
    //for (int i = 0; i < MN; i++) {
    //    x[i] = i;
    //    y[i] = i;
    //}
    //
    //x[4] = -30;
    //y[4] = -30;
    //x[3] = 50;
    //y[3] = 50;


    //hadamard2d<double>(x,N,M);
    //hadamardTransform2d<double>(y,M,N, SEQUENCY);
    //
    //for (int i = 0; i < MN; i++) {
    //    std::cout << x[i] - y[i] 
    //    << ", x: " << std::setw(10) << std::left << x[i] 
    //    << ", y: " << std::setw(10) << std::left << y[i] << std::endl;
    //}
    
    

}

    
    
    
    
    
    
    
    
    
    
    
    //timeit("Sequency    ", reps, N, hadamardSequency<double>);
    //timeit("Paley       ", reps, N, hadamardPaley<double>);
    //timeit("Ordinary    ", reps, N, hadamardOrdinary<double>);
    //timeit("Recursive   ", reps, N, hadamardRecursive<double>);
    //timeit("Depth First ", reps, N, hadamardDepthFirst<double>);
    //timeit("Arndt       ", reps, N, hadamardArndt<double>);
    //timeit("Parallel    ", reps, N, hadamardParallel<double>);
