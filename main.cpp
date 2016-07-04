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
#include "cycles.h"

int main(int argc, char *argv[]) {

    int nu = 25;

    //for (int nu = 4; nu < nu_max; nu++) {
    //    std::cout << "-------------------------------------------------------\n";

        const unsigned int N = powDyadic(nu);
        Timer t1;
        t1.set_name("deque");
        t1.start();
        std::deque<uint32_t> leaders = detectCycleLeaders(N);
        t1.stop();
        std::cout << "nu: " << nu << ", #leaders: " << leaders.size() << std::endl;
        //printCycleLeadersAndSize(leaders, N);
    //}
}

    
    
    
    
    
    
    
    
    
    
    
    //timeit("Sequency    ", reps, N, hadamardSequency<double>);
    //timeit("Paley       ", reps, N, hadamardPaley<double>);
    //timeit("Ordinary    ", reps, N, hadamardOrdinary<double>);
    //timeit("Recursive   ", reps, N, hadamardRecursive<double>);
    //timeit("Depth First ", reps, N, hadamardDepthFirst<double>);
    //timeit("Arndt       ", reps, N, hadamardArndt<double>);
    //timeit("Parallel    ", reps, N, hadamardParallel<double>);
