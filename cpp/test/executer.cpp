
#include "../cpufft.h"

using namespace cooleytukey;
void execute_fft()
{

    fft(std::vector<double> ({4,3,2,4})); 

}

int main()
{
    execute_fft();
    
}

