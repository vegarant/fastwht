#include "cycles.h"


std::deque<uint32_t> detectCycleLeaders(const uint32_t N) {
    std::vector<uint32_t> elements(N,0); 
    std::deque<uint32_t>  leaders;
    const uint32_t nu = findMostSignificantBit(N)-1;

    uint32_t idx = 1;
    elements[0] = 1;
    while(idx = findIdxOfNextZero(elements)) {
        const uint32_t starter = idx;
        elements[idx] = 1;
        idx = idxFromOrdinaryToSequency(idx, nu);
        while (idx != starter) {
            elements[idx] = 1;
            idx = idxFromOrdinaryToSequency(idx, nu);
        }

        leaders.push_back(starter);
    }

    return leaders;
}

uint32_t findIdxOfNextZero(const std::vector<uint32_t> & x) {
    
    uint32_t i = 1;
    const size_t N = x.size();
    
    while(x[i] and i < N) {
        i++;
    }
    
    if (i == N) {
        return 0;
    } else {
        return i;
    }
}


void printCycles(const std::deque<uint32_t> & x, const uint32_t N) {
    const uint32_t nu = findMostSignificantBit(N)-1;

    for( auto it = x.begin(); it != x.end(); it++) {    
        const uint32_t cycleLeader = *it;
        uint32_t idx = cycleLeader;
        
        do {
            std::cout << idx << " -> ";
            idx = idxFromOrdinaryToSequency(idx , nu);
        } while (idx != cycleLeader);
        std::cout << idx << std::endl;
    }
}


void printCycleLeadersAndSize(const std::deque<uint32_t> & x, const uint32_t N) {
    
    const uint32_t nu = findMostSignificantBit(N)-1;

    for( auto it = x.begin(); it != x.end(); it++) {    
        const uint32_t cycleLeader = *it;
        uint32_t idx = cycleLeader;
        unsigned int i = 0;
        do {
            i++;
            idx = idxFromOrdinaryToSequency(idx , nu);
        } while (idx != cycleLeader);
        std::cout << cycleLeader << ": " << i << std::endl;
    }
    
}



