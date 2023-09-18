#include <iostream>
#include <complex>
#include <vector>
#include <cmath>
#include <algorithm>
#include <string>

#define WRONGSIZE 48

#define M_PIO 3.1415926535897932384626433832795028841971693993751058209749445923078164062862089986280348253421170679
#define PREFACTOR 2.0

bool IsPowerOfTwo(unsigned long x)
{
    return (x != 0) && ((x & (x - 1)) == 0);
}
//constexpr
std::complex<double> k_nth_root(unsigned k, unsigned n)
{
    return std::exp((std::complex<double>((-1.0*PREFACTOR*k*M_PIO)/n) * std::complex<double>(0, 1)) );

}
std::complex<double> inverse_k_nth_root(unsigned k, unsigned n)
{
    return std::exp((std::complex<double>((PREFACTOR*k*M_PIO)/n) * std::complex<double>(0, 1)) );

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

                // if (size == 32)
                // {
                //     std::cout << std::endl; 
                //     for (unsigned i = 0; i < 10; i++)
                //     {
                //         std::cout << r[i] << "," ;
        
                //         }
                //         std::cout<<std::endl; 
                //     }

                const unsigned n = size; 
                const unsigned n_half = size/2; 
                std::complex<double> * even_indices = new std::complex<double>[n_half];
                std::complex<double> * odd_indices = new std::complex<double>[n_half];

                for (unsigned i = 0; i < n_half; i++)
                {
                    even_indices[i] = r[i*2]; //assign even indices
                    odd_indices[i] = r[i*2+1]; //assign uneven/odd indices
                }

                std::complex<double> * const processed_even = fftrecursive(even_indices, n_half);
                std::complex<double> * const processed_odd = fftrecursive(odd_indices, n_half);

                delete [] even_indices;
                delete [] odd_indices;

                std::complex<double> * c = new std::complex<double>[n];
                for ( unsigned i = 0; i < n_half; i++)
                {

                    std::complex<double> ur = k_nth_root(i, n);
                    c[i] = processed_even[i] + processed_odd[i] * ur;
                    c[i+n_half] = processed_even[i] - processed_odd[i] * ur;

                }
                delete [] processed_even;
                delete [] processed_odd;
                //delete [] r;

                return c; 

            }
        
    }


    std::complex<double> *
    ifftrecursive (std::complex<double> * r, unsigned size)
    {
            if (size == 1)
            {
                return r; 
            }
            else {
                const unsigned n = size; 
                const unsigned n_half = size/2; 
 
                std::complex<double> * even_indices = new std::complex<double>[n_half];
                std::complex<double> * odd_indices = new std::complex<double>[n_half];

                for (unsigned i = 0; i < n_half; i++)
                {
                    even_indices[i] = r[i*2];
                    odd_indices[i] = r[(i*2)+1];
                }

                std::complex<double> * const processed_even = ifftrecursive(even_indices, n_half);
                std::complex<double> * const processed_odd = ifftrecursive(odd_indices, n_half);

                delete [] even_indices;
                delete [] odd_indices;

                std::complex<double> * const c = new std::complex<double>[n];
               for ( unsigned i = 0; i < n_half; i++)
                {
                    std::complex<double> ur = k_nth_root(i, n);
                    c[i] = (processed_even[i] + processed_odd[i] * ur)/std::complex<double>(n, 0);
                    c[i+n_half] = (processed_even[i] - processed_odd[i] * ur)/std::complex<double>(n, 0);

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
            for (unsigned i = 0; i < 10; i++)
            {
                std::cout << r[i] << "," ;
            }
            std::cout << std::endl; 
         
            std::vector<std::complex<double>> c; 
            //c.reserve(r.size());
            for (unsigned i = 0; i < r.size(); i++)
            {
                c.push_back(std::complex<double>(r[i]));
                //std::cout << complex2string(c[i]) << ","; 
            }
            std::complex<double> * result = fftrecursive(c.data(), c.size());

            /*
            for (unsigned i = 0; i < 10; i++)
            {
                std::cout << result[i] << "," ;
            }
            std::cout << std::endl; 
            */
           return std::vector<std::complex<double>> (result, result + c.size());
            
        }
        else {

            std::cerr << "Could not process fft of r since the size of elements is not pow of 2." << std::endl; 
            throw 48; 


        }

    }



    //std::vector<double> 
    void
    ifft( 
        std::vector<std::complex<double>> r) {
        bool ispow2 = IsPowerOfTwo(r.size());
        if (ispow2)
        {
           std::complex<double> * result = ifftrecursive(r.data(), r.size());

            for (unsigned i = 0; i < 10; i++)
            {
                std::cout << result[i] << "," ;
            }


            
        }
        else {

            std::cerr << "Could not process ifft of r since the size of elements is not pow of 2." << std::endl; 
            throw 48; 


        }

    }


};