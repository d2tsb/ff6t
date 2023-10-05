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

#pragma once
#define WRONGSIZE 48
#define SIZENOTMATCH 44
#define ORDER 3
#define STATEMENT ((unsigned)((1 << ORDER) > 10) ? (10) : (1 << ORDER))

//#define M_PIO 3.141592653589793238462643383yyp2795028841971693993751058209749445923078164062862089986280348253421170679
#define M_PIO (double) 3.1415926535897932384626433832795028841971693993751058209749445923078164062862089986280348253421170679
#define PREFACTOR 2.0


std::vector<std::complex<double>> double_to_complex_conversion(std::vector<double> r)
{ //without iterator
        std::vector<std::complex<double>> c; 
                        for ( unsigned i = 0; i < r.size(); i++)
                        {
                            c.push_back(std::complex<double>(r[i]));
                        }
        return c; 

}
bool IsPowerOfTwo(unsigned long x)
{
    return (x != 0) && ((x & (x - 1)) == 0);
}
//constexpr
std::complex<double> k_nth_root(const unsigned k, const unsigned n)
{
    return std::exp((std::complex<double>(((-1.0*PREFACTOR*k*M_PIO)) * std::complex<double>(0, 1)))/std::complex<double>(n) );

}
std::complex<double> inverse_k_nth_root(const unsigned k, const unsigned n)
{
    return std::exp((std::complex<double>(((PREFACTOR*k*M_PIO)) * std::complex<double>(0, 1)))/std::complex<double>(n) );

}

std::string complex2string(std::complex<double> n)
{
    return std::to_string( n.real()) + " + " + std::to_string(n.imag()) + "j";
}


std::string int_as_string(unsigned i)
{
    std::string outcome; 
    long j = ((sizeof(unsigned) << 3) - 1);
    for ( ; j >= 0; j--)
    {
        outcome += std::to_string((bool)(i & (1 << j)));
    }
    //std::cout<< "leaving";
    return outcome; 

}

