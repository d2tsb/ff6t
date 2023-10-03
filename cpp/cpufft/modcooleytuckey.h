
#include "cooleytuckey.h"
#include "defines.h"
namespace cooleytukey 
    {
  
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
    namespace modrecursive337
    {
        //works for only 4 digits padded with 4 zeros   with first primenumber after 4*base(10)==9
        //8 the resulting inputlength is invertable
    namespace vec 
    {
        
        #define INTT int

        #define FINITE 337
        #define FINITEINV 295
        #define EIGHTYFIVEINV 226

        INTT mod337 (INTT b)
        {
            return b < FINITE && b >= 0 ? b : (b < 0 ? b + FINITE : b - FINITE); 
        }
        INTT fast_mod337 (INTT b)
        { //does not check whether b >= 0 because not needed in the (i)fft-task.
            return b < FINITE ? b : b - FINITE; 
        }


        bool checkunitrootproperty(INTT ur, INTT finite)
        {
            INTT sum = 0; 
            for ( unsigned i = 0; i < finite; i++)
            {
                
            }

            return sum ? false : true;

        }

        INTT *
        fftrecursive (INTT * r, unsigned size)
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
                        



                    INTT even_indices[n_half];
                    INTT odd_indices[n_half];

                    for (unsigned i = 0; i < n_half; i++)
                    {
                        even_indices[i] = r[i*2]; //assign even indices
                        odd_indices[i] = r[i*2+1]; //assign uneven/odd indices
                    }

                    INTT * processed_even = fftrecursive(even_indices, n_half);
                    INTT * processed_odd = fftrecursive(odd_indices, n_half);

                    // delete [] even_indices;
                    // delete [] odd_indices;

                    INTT * c = new INTT[n];
                    const INTT knr = 85;
                    INTT ur = 1;
                    for ( unsigned i = 0; i < n_half; i++)
                    {

                        c[i] = processed_even[i] + (processed_odd[i] * ur);
                        c[i+n_half] = processed_even[i] - (processed_odd[i] * ur);
                        ur *= knr; 
                        ur %= FINITE; 

                    }
                    delete [] processed_even;
                    delete [] processed_odd;
                    //delete [] r;

                    return c; 

                }
            
        }


        INTT *
        ifftrecursive (INTT * r, unsigned size)
        {
                if (size == 1)
                {
                    return r; 
                }
                else {
                    const unsigned n = size; 
                    const unsigned n_half = size>>1; 
    
                    INTT even_indices[n_half];
                    INTT odd_indices[n_half];

                    for (unsigned i = 0; i < n_half; i++)
                    {
                        even_indices[i] = r[i*2];
                        odd_indices[i] = r[(i*2)+1];
                    }

                    INTT * const processed_even = ifftrecursive(even_indices, n_half);
                    INTT * const processed_odd = ifftrecursive(odd_indices, n_half);

                    // delete [] even_indices;
                    // delete [] odd_indices;

                    INTT * const c = new INTT[n];
                    INTT ur = 1;
                    INTT knr = 85; //tbd
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



        std::vector<INTT>
        fft(std::vector<INTT> r) {
            bool ispow2 = IsPowerOfTwo(r.size());
            if (ispow2)
            {
                for (unsigned i = 0; i < STATEMENT; i++)
                {
                    std::cout << r[i] << "," ;
                }
                std::cout << std::endl; 
            
                INTT *  result = fftrecursive(r.data(), r.size());

                /*
                for (unsigned i = 0; i < 10; i++)
                {
                    std::cout << result[i] << "," ;
                }
                std::cout << std::endl; 
                */
                return std::vector<INTT> (result, result + r.size());
                
            }
            else {

                std::cerr << "Could not process fft of r since the size of elements is not pow of 2." << std::endl; 
                throw 48; 


            }

        }



        std::vector<INTT>
        ifft( 
            std::vector<INTT> r) {
            bool ispow2 = IsPowerOfTwo(r.size());
            if (ispow2)
            {
                INTT * result = ifftrecursive(r.data(), r.size());

                for ( unsigned i = 0; i < r.size(); i++)
                    {
                        result[i] = ( FINITEINV * result[i]) % FINITE; 
                    }

                for (unsigned i = 0; i < STATEMENT; i++)
                    {
                        std::cout << result[i] << "," ;
                    }
        
                std::vector<INTT> result_v;
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


       }   }


}