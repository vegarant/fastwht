

#ifndef FASTWHT_CYCLES_H_
#define FASTWHT_CYCLES_H_

#include <iostream>
#include <iomanip>
#include <complex>
#include <vector>
#include <deque>
#include <cstdint>
#include "hadamard.h"

uint32_t findIdxOfNextZero(const std::vector<uint32_t> & x);
std::deque<uint32_t> detectCycleLeaders(const uint32_t N);
void printCycles(const std::deque<uint32_t> & x, const uint32_t N);
void printCycleLeadersAndSize(const std::deque<uint32_t> & x, const uint32_t N);




#endif



