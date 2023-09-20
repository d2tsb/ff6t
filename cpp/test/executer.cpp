
#include "../cpufft.h"
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
            std::cout << microseconds.count() << "µs\n";
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

void execute_vec_cooleytuckey_fft()
{
    std::cout << "Test of (recursive) cooleytuckey::recursive::vec(tor)::fft/ifft" << std::endl;
    std::srand(std::time(0));
    //std::vector< double > v({ 1, 1 }); 
    //std::vector< double > v; 
    std::vector< double > v({ 1, 2, 3, 4, 5, 6, 7, 8 }); 

    // for ( unsigned i = 0; i < 1 << ORDER; i++)
    // {
    //     v.push_back(std::rand()%10000);
    // }
    
    
    /*

    v.resize(1<< ORDER);
     for (std::vector<double>::iterator it = v.begin(); it != v.end(); ++it)
    {
       *it = rand();
    }
    
    
    */
    //std::cout << fft(v).size(); o
    std::vector<std::complex<double>> cv = cooleytukey::recursive::vec::fft(v);
        for (unsigned i = 0; i < STATEMENT; i++)
                {
                    std::cout << cv[i] << "," ;
                }
        std::cout << "fft done" << std::endl;
    cooleytukey::recursive::vec::ifft(cv);
        
    std::cout << std::endl; 
    
   //std::cout << "Statement: " << STATEMENT; 
}

void execute_arr_cooleytuckey_fft()
{

    std::cout << "Test of (recursive) cooleytuckey::recursive::arr(std::array)::fft/ifft" << std::endl;
    std::srand(std::time(0));
    //std::vector< double > v({ 1, 1 }); 
    std::array< double, 1 << ORDER> v; 
    //std::vector< double > v({ 1, 2, 3, 4, 5, 6, 7, 8 }); 

    for ( unsigned i = 0; i < 1 << ORDER; i++)
    {
        v[i] = (std::rand()%10000);
    }
    
    
    /*

    v.resize(1<< ORDER);
     for (std::vector<double>::iterator it = v.begin(); it != v.end(); ++it)
    {
       *it = rand();
    }
    
    
    */
    //std::cout << fft(v).size(); o
    std::array<std::complex<double>, 1 << ORDER> cv = cooleytukey::recursive::arr::fft(v);
        for (unsigned i = 0; i < STATEMENT; i++)
                {
                    std::cout << cv[i] << "," ;
                }
        std::cout << "fft done" << std::endl;
    cooleytukey::recursive::arr::ifft(cv);
        
    std::cout << "ifft done\n";
    
   //std::cout << "Statement: " << STATEMENT; 
}


void execute_vec_cooleytuckey_fft_iterative()
{
    std::cout << "Test of (iterative) cooleytuckey::iterative::vec(std::vector)::fft/ifft" << std::endl;
    /*
    std::cout << "sizeof(unsinged):" << sizeof(unsigned) << std::endl; 
    std::cout << "bitinverse of 0:" << inverse(0, 4) << std::endl;
    std::cout << int_as_string(inverse(0,4)) << std::endl; 
    std::cout << int_as_string(7) << std::endl; 
    std::cout << int_as_string(8) << std::endl; 
    std::cout << int_as_string(inverse(8,4)) << std::endl; 
    */
    //std::vector<double> v({1,1});
    std::vector<double> v({1,2,3,4,5,6,7,8});
    auto r = cooleytukey::iterative::vec::fft(v);
    for (unsigned i = 0; i < 8; i++) 
    {
        std::cout << r[i] << ","; 
    }
    std::cout << "fft done\n";
}



int main()
{
    // Stopwatch sw; 
    // sw.start();
    //execute_arr_cooleytuckey_fft();
    // sw.finish();
    // sw.print_duration_in_milliseconds(); 

    // sw.reset(); 
    // sw.start();
    execute_vec_cooleytuckey_fft();
    // sw.finish();
    // sw.print_duration_in_milliseconds(); 
    execute_vec_cooleytuckey_fft_iterative();
    
}

