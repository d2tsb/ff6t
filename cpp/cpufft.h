#include <iostream>
#include <complex>
#include <vector>
#include <cmath>
#include <algorithm>
#include <string>

#define WRONGSIZE 48

#define M_PIO 3.1415926535897932384626433832795028841971693993751058209749445923078164062862089986280348253421170679

bool IsPowerOfTwo(unsigned long x)
{
    return (x != 0) && ((x & (x - 1)) == 0);
}
//constexpr
std::complex<double> k_nth_root(unsigned k, unsigned n)
{
    return std::exp((std::complex<double>(-2.0*k*M_PIO/n) * std::complex<double>(0, 1)) );

}

std::string complex2string(std::complex<double> n)
{
    return std::to_string( n.real()) + " + " + std::to_string(n.imag()) + "j";
}

namespace cooleytukey
{


    std::complex<double> *
    fftrecursive (std::complex<double> * r, unsigned size)
    {
            if (size == 1)
            {
                return r; 
            }
            else {
                std::complex<double> * even_indices = new std::complex<double>[size/2];
                std::complex<double> * odd_indices = new std::complex<double>[size/2];

                for (unsigned i = 0; i < size/2; i++)
                {
                    even_indices[i] = r[i*2];
                    odd_indices[i] = r[i*2+1];
                }

                std::complex<double> * const processed_even = fftrecursive(even_indices, size/2);
                std::complex<double> * const processed_odd = fftrecursive(odd_indices, size/2);

                delete [] even_indices;
                delete [] odd_indices;

                std::complex<double> * c = new std::complex<double>[size];
                const unsigned n_half = size/2; 
                for ( unsigned i = 0; i < size/2; i++)
                {
                    const std::complex<double> knr = k_nth_root(i, size);
                   c[i] = processed_even[i] + processed_odd[i] * knr;
                   c[i+n_half] = processed_even[i] + processed_odd[i] * knr;

                }
                delete [] processed_even;
                delete [] processed_odd;
                //delete [] r;

                return c; 

            }
        
    }


    std::vector<std::complex<double>>
    fft(std::vector<double> r) {
        bool ispow2 = IsPowerOfTwo(r.size());
        if (ispow2)
        {
            std::vector<std::complex<double>> c; 
            c.reserve(r.size());
            for (unsigned i = 0; i < r.size(); i++)
            {
                c.push_back(std::complex<double>(r[i]));
                //std::cout << complex2string(c[i]) << ","; 
            }
            std::complex<double> * result = fftrecursive(c.data(), c.size());

            for (unsigned i = 0; i < 10; i++)
            {
                std::cout << result[i] << "," ;
            }


            
        }
        else {

            std::cerr << "Could not process r since the size of elements is not pow of 2." << std::endl; 
            throw 48; 


        }

    }





};