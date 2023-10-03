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







    namespace iterative

    {

        namespace vec
        {
            std::vector < std::complex<double> > fft 
            (std::vector<double> c)
            {
                if ( IsPowerOfTwo(c.size()))
                {
                    unsigned o = std::log2(c.size());
                    //conversion to complex<double> array
                    std::vector<std::complex<double>> r; 
                    for ( unsigned i = 0; i < c.size(); i++)
                    {
                        r.push_back(std::complex<double>(c[inverse(i, o)])); //fed with bitinverse index
                        //r.push_back(std::complex<double>(c[i])); //fed with bitinverse index
                    }
                    for ( unsigned i = 0; i < r.size(); i++)
                        std::cout << r[i] << ","; 
                    std::cout << std::endl;

                    //perform iterative fft
                    unsigned N_SIZE = c.size(); 
                    unsigned Partitionsize = 2; //2
                    for (unsigned i = 0; i < o; i++)
                    {
                        for (unsigned Partition = 0; Partition<N_SIZE; Partition+= Partitionsize )
                            {
                                std::complex<double> ur = k_nth_root(1,Partitionsize);
                                std::complex<double> multiplier(1,0); 
                                //extract uneven and even
                                std::complex<double> Even[Partitionsize/2]; //2
                                std::complex<double> Uneven[Partitionsize/2]; //2
                                for (unsigned k=0; k<Partitionsize/2; k++)  //2
                                {
                                    Even[k] = r[2*k+Partition];
                                    Uneven[k] = r[2*k+ Partition + 1];
                                }
                                for(unsigned k=0; k<Partitionsize/2; k++)  //2
                                {
                                    r[k+Partition] = Even[k] + Uneven[k]*multiplier;  
                                    r[k+Partition+Partitionsize/2] = Even[k] - Uneven[k]*multiplier; 
                                    multiplier *= ur; 
                                    
                                }
                            }
                                                            
                        Partitionsize<<=1;
                    }
                    return r; 

                }
                else {
                    throw WRONGSIZE; 
                }
                

            }

        std::vector < std::complex<double> > ifft
                (std::vector<std::complex<double>> c)
                {
                    if ( IsPowerOfTwo(c.size()))
                    {
                        unsigned o = std::log2(c.size());
                        //conversion to complex<double> array
                        std::vector<std::complex<double>> r; 
                        for ( unsigned i = 0; i < c.size(); i++)
                        {
                            //r.push_back(std::complex<double>(c[inverse(i, o)])); //fed with bitinverse index
                            r.push_back(std::complex<double>(c[i])); //fed with bitinverse index
                        }
                                         //perform iterative fft
                        unsigned N_SIZE = c.size(); 
                        unsigned Partitionsize = 2; //2
                        for (unsigned i = 0; i < o; i++)
                        {
                            for (unsigned Partition = 0; Partition<N_SIZE; Partition+= Partitionsize )
                                {
                                    std::complex<double> ur = inverse_k_nth_root(1,Partitionsize);
                                    std::complex<double> multiplier(1,0); 
                                    //extract uneven and even
                                    std::complex<double> Even[Partitionsize/2]; //2
                                    std::complex<double> Uneven[Partitionsize/2]; //2
                                    for (unsigned k=0; k<Partitionsize/2; k++)  //2
                                    {
                                        Even[k] = r[2*k+Partition];
                                        Uneven[k] = r[2*k+ Partition + 1];
                                    }
                                    for(unsigned k=0; k<Partitionsize/2; k++)  //2
                                    {
                                        r[k+Partition] = Even[k] + Uneven[k]*multiplier;  
                                        r[k+Partition+Partitionsize/2] = Even[k] - (Uneven[k]*multiplier); 
                                        multiplier *= ur; 
                                        
                                    }
                                }
                                                                
                            Partitionsize<<=1;
                        }


                        for (unsigned i = 0; i < r.size(); i++)
                        {
                            r[i] /= r.size(); 
                        }

                        return r; 

                    }
                    else {
                        throw WRONGSIZE; 
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