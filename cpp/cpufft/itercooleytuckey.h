#pragma once
#include "defines.h"

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