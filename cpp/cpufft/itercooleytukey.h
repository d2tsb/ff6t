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
                    for ( unsigned i = 0; i < r.size(); i++)
                        std::cout << r[i] << ","; 
                    std::cout << std::endl;

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
                            //r.push_back(std::complex<double>(c[bitinverse(i, o)])); //fed with bitinverse index
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
}