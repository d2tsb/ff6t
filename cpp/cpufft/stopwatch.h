#include <iostream>
#include <complex>
#include <vector>
#include <cmath>
#include <algorithm>
#include <array>
#include <string>
#include <climits>
#include <cstdlib>
#include <ctime>
#include <chrono>


class Stopwatch 
{
    private:
        std::chrono::_V2::system_clock::time_point start_, finish_; 
    
    public: 
        Stopwatch()
        {
            reset(); 
        }
        void start () 
        {
            this->start_ = std::chrono::high_resolution_clock::now();
        }

        void finish() 
        {
            this->finish_ = std::chrono::high_resolution_clock::now();
        }
        void print_duration_in_microseconds()
        {
            auto microseconds = std::chrono::duration_cast<std::chrono::microseconds>(finish_ - start_);
            std::cout << microseconds.count() << " microseconds\n"; //µ
        }
        void print_duration_in_milliseconds()
        {
            auto milliseconds = std::chrono::duration_cast<std::chrono::milliseconds>(finish_ - start_);
            std::cout << milliseconds.count() << "ms\n";
        }
 
        void reset () 
        {
            std::chrono::_V2::system_clock::time_point now = std::chrono::high_resolution_clock::now();
            this->start_ = this->finish_ = now; 
        }

};

