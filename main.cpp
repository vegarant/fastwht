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

using duration = std::chrono::duration<double>;

int main(int argc, char *argv[]) {

    int nu = 23;
    std::cout << "N = 2^" << nu << std::endl;
    
    const unsigned int N = powDyadic(nu);

    timeit("Sequency    ", 7, N, hadamardSequency<double>);
    timeit("Paley       ", 7, N, hadamardPaley<double>);
    timeit("Ordinary    ", 7, N, hadamardOrdinary<double>);
    timeit("Recursive   ", 7, N, hadamardRecursive<double>);
    timeit("Depth First ", 7, N, hadamardDepthFirst<double>);
    timeit("Matters     ", 7, N, hadamardMatters<double>);

}

