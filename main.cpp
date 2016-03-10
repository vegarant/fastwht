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









