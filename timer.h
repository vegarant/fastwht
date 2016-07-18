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

This file contain a complete timer class, in addition to a few timer functions.  

*/

#ifndef FASTWHT_TIMER_H_
#define FASTWHT_TIMER_H_

#include <cstdint>
#include <cstring>
#include <cmath>

#include <algorithm>
#include <iostream>
#include <chrono>
#include <random>
#include <vector>
#include <tuple>
#include <iomanip>

/*

This class tries to simplifies the time measurements of the program by 
implementing start, pause and stop functions.

*/
class Timer {

    std::chrono::time_point<std::chrono::high_resolution_clock> start_time, end_time;
    std::chrono::duration<double> elapsed_time;
    char * name;
    bool is_active;
    protected:

    public:

    Timer() {
        name = new char[100];
        name[0] = '\0';
        is_active = false;
    }

    ~Timer() {
        delete [] name;    
    }

    void set_name(const char * new_name);
    void start();
    void pause();
    void stop();

};

////////////////////////////////////////////////////////////////////////////////

std::tuple<double, double, double, double> statistics(const std::vector<double> &val);

void timeit( const char * name, const int runs, const unsigned long N, 
             void (*transform)(double*, const unsigned long));

#endif

