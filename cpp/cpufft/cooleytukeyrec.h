#pragma once
#include "defines.h"
namespace cooleytukey 
{
    namespace recursive
    {
    namespace vec 
    {


        std::complex<double> *
        fftrecursive (std::complex<double> * r, unsigned size)
        {
                if (size == 1)
                {
                    return r; 
                }
                else {

                    // if (size == (1 << ORDER))
                    // {
                    //     std::cout << std::endl; 
                    //     for (unsigned i = 0; i < (1<<ORDER); i++)
                    //     {
                    //         std::cout << r[i] << "," ;
            
                    //         }
                    //         std::cout<<std::endl; 
                    // }

                    const unsigned n = size; 
                    const unsigned n_half = size>>1; 

                    // for (unsigned i = 0; i < n; i++)
                    //     {
                    //         std::cout << r[i] << "," ;
            
                    //         }
                        



                    std::complex<double> even_indices[n_half];
                    std::complex<double> odd_indices[n_half];

                    for (unsigned i = 0; i < n_half; i++)
                    {
                        even_indices[i] = r[i*2]; //assign even indices
                        odd_indices[i] = r[i*2+1]; //assign uneven/odd indices
                    }

                    std::complex<double> * processed_even = fftrecursive(even_indices, n_half);
                    std::complex<double> * processed_odd = fftrecursive(odd_indices, n_half);

                    // delete [] even_indices;
                    // delete [] odd_indices;

                    std::complex<double> * c = new std::complex<double>[n];
                    auto knr = k_nth_root(1, n);
                    std::complex<double> ur(1, 0);
                    for ( unsigned i = 0; i < n_half; i++)
                    {

                        c[i] = processed_even[i] + (processed_odd[i] * ur);
                        c[i+n_half] = processed_even[i] - (processed_odd[i] * ur);
                        ur *= knr; 

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
                    const unsigned n_half = size>>1; 
    
                    std::complex<double> even_indices[n_half];
                    std::complex<double> odd_indices[n_half];

                    for (unsigned i = 0; i < n_half; i++)
                    {
                        even_indices[i] = r[i*2];
                        odd_indices[i] = r[(i*2)+1];
                    }

                    std::complex<double> * const processed_even = ifftrecursive(even_indices, n_half);
                    std::complex<double> * const processed_odd = ifftrecursive(odd_indices, n_half);

                    // delete [] even_indices;
                    // delete [] odd_indices;

                    std::complex<double> * const c = new std::complex<double>[n];
                    std::complex<double> ur(1, 0);
                    auto knr = inverse_k_nth_root(1, n);
                    for ( unsigned i = 0; i < n_half; i++)
                    {
                        c[i] = (processed_even[i] + (processed_odd[i] * ur));
                        c[i+n_half] = (processed_even[i] - (processed_odd[i] * ur));
                        ur *= knr; 
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
                for (unsigned i = 0; i < STATEMENT; i++)
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



        std::vector<std::complex<double>>
        ifft( 
            std::vector<std::complex<double>> r) {
            bool ispow2 = IsPowerOfTwo(r.size());
            if (ispow2)
            {
                std::complex<double> * result = ifftrecursive(r.data(), r.size());

                std::complex<double> quotient(r.size(), 0);
                for ( unsigned i = 0; i < r.size(); i++)
                    {
                        result[i] /= quotient; 
                    }

                for (unsigned i = 0; i < STATEMENT; i++)
                    {
                        std::cout << result[i] << "," ;
                    }
        
                std::vector<std::complex<double>> result_v;
                for ( unsigned i = 0; i < r.size(); i++)
                    {
                        result_v.push_back(result[i]);
                    }
                return result_v;

            }
            else {

                std::cerr << "Could not process ifft of r since the size of elements is not pow of 2." << std::endl; 
                throw 48; 


            }

        }


    }

    namespace arr 
    {

    template <size_t S>
    std::array<std::complex<double>, S>
    fftrecursive (std::array<std::complex<double>, S> r)
    {
            if (r.size() == 1)
            {
                return r; 
            }
            else {

                // if (size == (1 << ORDER))
                // {
                //     std::cout << std::endl; 
                //     for (unsigned i = 0; i < (1<<ORDER); i++)
                //     {
                //         std::cout << r[i] << "," ;
        
                //         }
                //         std::cout<<std::endl; 
                // }

                const unsigned n = r.size(); 
                const unsigned n_half = n>>1; 

                // for (unsigned i = 0; i < n; i++)
                //     {
                //         std::cout << r[i] << "," ;
        
                //         }
                    
                std::array<std::complex<double>, n_half> even_indices;
                std::array<std::complex<double>, n_half> odd_indices;

                for (unsigned i = 0; i < n_half; i++)
                {
                    even_indices[i] = r[i*2]; //assign even indices
                    odd_indices[i] = r[i*2+1]; //assign uneven/odd indices
                }

                std::array<std::complex<double>, n_half> processed_even = fftrecursive(even_indices); // with length n_half
                std::array<std::complex<double>, n_half> processed_odd = fftrecursive(odd_indices); // with length n_half

                // delete [] even_indices;
                // delete [] odd_indices;

                std::array<std::complex<double>, n> c;
                auto knr = k_nth_root(1, n);
                std::complex<double> ur(1, 0);
                for ( unsigned i = 0; i < n_half; i++)
                {

                    c[i] = processed_even[i] + (processed_odd[i] * ur);
                    c[i+n_half] = processed_even[i] - (processed_odd[i] * ur);
                    ur *= knr; 

                }
                //delete & processed_even;
                //delete & processed_odd;
                //delete [] r;

                return c; 

            }
        
    }


    template <size_t S>
    std::array<std::complex<double>, S>
    ifftrecursive (std::array<std::complex<double>, S> r)
    {
            if (r.size() == 1)
            {
                return r; 
            }
            else {
                const unsigned n = r.size(); 
                const unsigned n_half = n>>1; 
 
                std::array<std::complex<double>, n_half> even_indices;
                std::array<std::complex<double>, n_half> odd_indices;

                for (unsigned i = 0; i < n_half; i++)
                {
                    even_indices[i] = r[i*2];
                    odd_indices[i] = r[(i*2)+1];
                }

                std::array<std::complex<double>, n_half> processed_even = ifftrecursive(even_indices);
                std::array<std::complex<double>, n_half> processed_odd = ifftrecursive(odd_indices);
                // delete [] odd_indices;

                std::array<std::complex<double>, n> c;
                std::complex<double> ur(1, 0);
                auto knr = inverse_k_nth_root(1, n);
                for ( unsigned i = 0; i < n_half; i++)
                {
                    c[i] = (processed_even[i] + (processed_odd[i] * ur));
                    c[i+n_half] = (processed_even[i] - (processed_odd[i] * ur));
                    ur *= knr; 
                }
                //delete & processed_even;
                //delete & processed_odd;
                //delete [] r;

                return c; 

            }
        
    }




   
    template <size_t S>
    std::array<std::complex<double>, S>
    fft(std::array<double, S> r) {
        bool ispow2 = IsPowerOfTwo(r.size());
        if (ispow2)
        {
            for (unsigned i = 0; i < STATEMENT; i++)
            {
                std::cout << r[i] << "," ;
            }
            std::cout << std::endl; 
         
            std::array<std::complex<double>, S> c; 
            for (unsigned i = 0; i < r.size(); i++)
            {
                c[i] = std::complex<double>(r[i]);
                //std::cout << complex2string(c[i]) << ","; 
            }
            std::array<std::complex<double>, S> result = fftrecursive(c);
            return result;
            
        }
        else {

            std::cerr << "Could not process fft of r since the size of elements is not pow of 2." << std::endl; 
            throw 48; 


        }

    }



    //std::vector<double> 
    template < size_t S >
    void
    ifft( 
        std::array<std::complex<double>, S> r) {
        bool ispow2 = IsPowerOfTwo(r.size());
        if (ispow2)
        {
            std::array<std::complex<double>, r.size()> result = ifftrecursive(r);

            std::complex<double> quotient(r.size(), 0);
            for ( unsigned i = 0; i < r.size(); i++)
                {
                    result[i] /= quotient; 
                }

            for (unsigned i = 0; i < STATEMENT; i++)
                {
                    std::cout << result[i] << "," ;
                }
    

            
        }
        else {

            std::cerr << "Could not process ifft of r since the size of elements is not pow of 2." << std::endl; 
            throw 48; 


        }

    }


    }
       }
    namespace convolution_recursive
    {
        std::vector <std::complex<double>> multiplication_via_fft(std::vector<double> a, 
                std::vector<double> b)
        {

            if ( a.size() != b.size() )
            {
                throw SIZENOTMATCH; 
            }

            if ( IsPowerOfTwo(a.size())
                && IsPowerOfTwo(b.size()) 
            )
            {
                try
                {
                    auto a_coeff = recursive::vec::fft(a);
                    auto b_coeff = recursive::vec::fft(b);
                    for ( unsigned i = 0; i < a_coeff.size(); i++)
                        {
                            a_coeff[i] *= b_coeff[i];
                        }

                    std::vector <std::complex<double>> c_ifft = recursive::vec::ifft(a_coeff);

                    return c_ifft;
                    }

                catch(const std::exception& e)
                {
                    std::cerr << e.what() << '\n';
                }
                
                




            }
            else {
                        throw WRONGSIZE; 
                    }
        }


    }
  
       
       
       }


