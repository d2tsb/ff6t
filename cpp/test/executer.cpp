
#include "../cpufft.h"
#include <cstdlib>
#include <ctime>


using namespace cooleytukey;
void execute_fft()

{

    std::srand(std::time(0));
    //std::vector< double > v({ 1, 1 }); 
    std::vector< double > v; 


    //std::vector< double > v({ 1, 2, 3, 4, 5, 6, 7, 8 }); 

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
    //std::cout << fft(v).size(); o




    std::vector<std::complex<double>> cv = fft(v);
        for (unsigned i = 0; i < STATEMENT; i++)
                {
                    std::cout << cv[i] << "," ;
                }
        std::cout << "fft done" << std::endl;
    cooleytukey::ifft(cv);
        
    
   //std::cout << "Statement: " << STATEMENT; 
}

int main()
{
    execute_fft();
    
}

