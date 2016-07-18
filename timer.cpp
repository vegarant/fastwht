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

#include "timer.h"

/*

timeit is a function specially designed to test the computational time of 
various transforms.

timeit test the computational time of the function `transform` running it
`runs` number of times. The timeit function will output the computational
statistics from these runs to stdout. The `name` will identify the name of the
function being called. The length of the input vector to the transform
function will be `N`.

INPUT:
name - String identifying which functions is being called. 
runs - Number of times the `transform` functions should be called.
N    - Length of the vector the `transform` function should transform.
transform - The transform function. The first argument to this function will be 
            a random vector of doubles it should transform, the second argument 
            is the length of the transform.

*/
void timeit( const char * name, const int runs, const unsigned long N, 
             void (*transform)(double*, const unsigned long))
{
    // Create vector with random numbers
    std::vector<double> s(N);

    std::mt19937 gen(0);
    std::uniform_real_distribution<double> dis(0, 1);

    std::generate(s.begin(), s.end(), [ &dis, &gen ] () { return dis(gen); });

    // Prepare for time measurements
    std::vector<double> time_measurements; 
    std::chrono::time_point<std::chrono::high_resolution_clock> 
                                               start_time_local, end_time_local;

    for (int i = 0; i < runs; i++) {
        start_time_local = std::chrono::high_resolution_clock::now();    
        transform( &s[0], N );
        end_time_local = std::chrono::high_resolution_clock::now();    
        time_measurements.push_back(std::chrono::duration_cast<
                                    std::chrono::duration<double> >
                                    (end_time_local - start_time_local).count());
    }
    double min=0, mean=0, median=0, sd=0;
    std::tie(min,mean, median, sd) = statistics(time_measurements);
    
    int spaceing = 8;
    std::cout << std::left << std::setw(10) << name 
              << " Min: "    << std::left << std::setw(spaceing) << min 
              << " Mean: "   << std::left << std::setw(spaceing) << mean 
              << " Median: " << std::left << std::setw(spaceing) << median 
              << " sd: "     << std::left << std::setw(spaceing) << sd << std::endl; 

}

/*

statistics computes the computational statistics found in the vector `val`. The 
following entries are being returned as 4-tuple:
min value, mean, median and standard deviation. 

*/
std::tuple<double, double, double, double> statistics(const std::vector<double> &val) 
{
    size_t N = val.size();
    double mean = 0, sd = 0, median = 0, min = 0;

    for(double v : val) mean += v / N;

    for(double v : val) sd += (v - mean)*(v - mean);

    sd = std::sqrt(sd/N);

    std::vector<double> val_sorted(val);
    std::sort(val_sorted.begin(), val_sorted.end());
    min = val_sorted[0];
    median = (N%2) ? val_sorted[N/2+1] : (val_sorted[N/2-1] + val_sorted[N/2-1])/2;

    return std::make_tuple(min, mean, median, sd);

}

/*

Start the time measurement.

*/
void Timer::start() {
    if (is_active) {
        std::cerr << "Using Timer::start() wrong" << std::endl;
    } else {
        start_time = std::chrono::high_resolution_clock::now();    
        is_active = true;
    }
}

/*

Pauses the time measurement.

*/
void Timer::pause() {
    if (is_active) {
        end_time = std::chrono::high_resolution_clock::now();    
        elapsed_time += end_time - start_time;
        is_active = false;
    } else {

        std::cerr << "Using Timer::end() wrong" << std::endl;

    }
}

/*

Sets the name of the timer. It is advisable to set the name of the timer before 
one call the stop() function, to get a more intuitive output. 

*/
void Timer::set_name(const char * new_name) {
    std::strcpy(name, new_name); 
}

/*

Stop the time measurement. The function will write the elapsed time to stdout. 

*/
void Timer::stop() {
    if (is_active) {
        pause();    
    }

    if (name[0] == '\0') {
        std::cout << "Elapsed time: " << elapsed_time.count() << " s" 
                  << std::endl;    
    } else {

        std::cout << "Time used by " << name << ": " << elapsed_time.count() 
                  << " s" << std::endl;    
    } 
}

