#include <iostream>
#include <complex>
#include <vector>
#include <cmath>
#include <algorithm>

#define WRONGSIZE 48

const double M_PI = atan(1)*4; 


bool IsPowerOfTwo(unsigned long x)
{
    return (x != 0) && ((x & (x - 1)) == 0);
}
constexpr std::complex<double> k_nth_root(unsigned k, unsigned n)
{
    return std::exp((std::complex<double>(-2.0*k*M_PI/n) * std::complex<double>(0, 1)) );

}

namespace cooleytukey
{



    std::vector<std::complex<double>>
    fft(std::vector<double> r) {
        bool ispow2 = IsPowerOfTwo(r.size());
        if (ispow2)
        {
            std::vector<std::complex<double>> c; 
            std::copy_n(r.begin(), r.size(), c.begin());


        }
        else {

            std::cerr << "Could not process r since the size of elements is not pow of 2." << std::endl; 
            throw 48; 


        }

    }



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

                return c; 

            }
        
    }



};