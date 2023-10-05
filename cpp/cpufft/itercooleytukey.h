#include "defines.h"

namespace cooleytukey{
    namespace iterative

    {
        unsigned bitreversechar(unsigned i, short order = ORDER)
        {
            return i = (i * 0x0202020202ULL & 0x010884422010ULL) % 1023;


        }
        unsigned bitreverse(unsigned x, short order = ORDER)
        //bit reverse scoped in order
        {
                x = (((x & 0xaaaaaaaa) >> 1) | ((x & 0x55555555) << 1));
                x = (((x & 0xcccccccc) >> 2) | ((x & 0x33333333) << 2));
                x = (((x & 0xf0f0f0f0) >> 4) | ((x & 0x0f0f0f0f) << 4));
                x = (((x & 0xff00ff00) >> 8) | ((x & 0x00ff00ff) << 8));
                return((x >> 16) | (x << 16)) >> ((sizeof(int)<<3)-order);

        }


        
        unsigned bitinverse(unsigned i, short order = ORDER)
        {
            return 
            !((UINT_MAX << order) & i) *  //sharp bitinverse
            (((~i) << ((sizeof(unsigned) << 3) - order)) >> ((sizeof(unsigned) * 8) - order));
        }



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
                        r.push_back(std::complex<double>(c[bitreverse(i, o)])); //fed with bitinverse index
                        //r.push_back(std::complex<double>(c[i])); //fed with bitinverse index
                    }
                    // for ( unsigned i = 0; i < r.size(); i++)
                    //     std::cout << r[i] << ","; 
                    // std::cout << std::endl;

                    //perform iterative fft
                    unsigned N_SIZE = c.size(); 
                    unsigned Partitionsize = 2; //2
                    for (unsigned i = 0; i < o; i++)
                    {
                        std::complex<double> ur = k_nth_root(1,Partitionsize);
                        for (unsigned Partition = 0; Partition<N_SIZE; Partition+= Partitionsize )
                            {
                                std::complex<double> multiplier(1,0); 
                                //extract uneven and even
                                for(unsigned k=0; k<Partitionsize/2; k++)  //2
                                {

                                    std::complex<double> t = multiplier * r[Partition + k + Partitionsize/2];
                                    std::complex<double> u = r[Partition + k];
                                    r[k+Partition] = u + t;  
                                    r[k+Partition+Partitionsize/2] = u - t; 
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
                            r.push_back(std::complex<double>(c[bitreverse(i, o)])); //fed with bitinverse index
                        }
                                         //perform iterative fft
                        unsigned N_SIZE = c.size(); 
                        unsigned Partitionsize = 2; //2
                        for (unsigned i = 0; i < o; i++)
                        {

                            std::complex<double> ur = inverse_k_nth_root(1,Partitionsize);
                                for (unsigned Partition = 0; Partition<N_SIZE; Partition+= Partitionsize )
                                    {
                                        std::complex<double> multiplier(1,0);  //twiddle
                                        //extract uneven and even
                                        for(unsigned k=0; k<Partitionsize/2; ++k)  //2
                                        {

                                            std::complex<double> t = multiplier * r[Partition + k + Partitionsize/2];
                                            std::complex<double> u = r[Partition + k];
                                            r[k+Partition] = u + t;  
                                            r[k+Partition+Partitionsize/2] = u - t; 
                                            multiplier *= ur; 
                                            
                                        }
                                    }
                                                                        
                            Partitionsize<<=1;
                        }


                        for (unsigned i = 0; i < r.size(); i++)
                        {
                            r[i] /= N_SIZE; 
                        }

                        return r; 

                    }
                    else {
                        throw WRONGSIZE; 
                    }
                    

            }

            void benchmark (const unsigned order) 
            {

                    double *v = new double[1 << order];
                    for ( unsigned i = 0; i < 1 << order; i++)
                    {
                        v[i] = (std::rand()%10000);
                    }
                    std::vector<double> values; 
                    values.insert(values.end(), &v[0], &v[1 << order]);
                    //std::copy(&v[0], &v[1<<order], back_inserter(values));
                    std::cout << "benchmark for ..iterative::vec::fft: " << std::endl;


                    Stopwatch sw; 
                    sw.reset(); 
                    sw.start(); 
                    std::vector<std::complex<double>> twiddles = cooleytukey::iterative::vec::fft(values);
                    sw.finish(); 
                    std::cout << "\telapsed time for cooleytukey::iterative::vec::fft with " << (1 << order) << " elements: "; 
                    sw.print_duration_in_milliseconds(); 
                    std::cout << std::endl; 


                    sw.reset(); 
                    sw.start(); 
                    std::vector<std::complex<double>> results = cooleytukey::iterative::vec::ifft(twiddles);
                    sw.finish(); 
                    std::cout << "\telapsed time for cooleytukey::iterative::vec::ifft with " << (1 << order) << " elements: "; 
                    sw.print_duration_in_milliseconds(); 
                    std::cout << std::endl; 

                    std::cout << "\tcomparing some elements: " << std::endl << "\t\toriginal sample: "; 

                    for ( unsigned i = 0; i < 4; ++i)
                    {
                        std::cout << v[i] << ",";
                    }
                    std::cout << std::endl; 
                    


                    std::cout << "\tcomparing some elements: " << std::endl << "\t\tifft sample: "; 

                    for ( unsigned i = 0; i < 4; ++i)
                    {
                        std::cout << results[i] << ",";
                    }
                    std::cout << std::endl; 
                    
            }

        }



    }
}