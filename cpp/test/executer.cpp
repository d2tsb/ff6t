
#include "../cpufft.h"
#include <cstdlib>
#include <ctime>


using namespace cooleytukey;
void execute_fft()
#define ORDER 5

{

    std::srand(std::time(0));
    std::vector< double > v; 
    v.reserve(1<<ORDER);
    for ( unsigned i = 0; i < 1 << ORDER; i++)
    {
        v.push_back(std::rand()%10000);
    }
    
    
    /*

    v.resize(1<< ORDER);
     for (std::vector<double>::iterator it = v.begin(); it != v.end(); ++it)
    {
       *it = rand();
    }
    
    
    */
    //std::cout << fft(v).size(); 
    std::vector<std::complex<double>> cv = fft(v);
       for (unsigned i = 0; i < 10; i++)
            {
                std::cout << cv[i] << "," ;
            }
    std::cout << "fft done" << std::endl;
    cooleytukey::ifft(cv);
}

int main()
{
    execute_fft();
    
}

