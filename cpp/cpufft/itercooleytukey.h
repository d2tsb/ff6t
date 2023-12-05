#include "defines.h"

namespace cooleytukey{
    namespace iterative

    {
        unsigned bitreversechar(unsigned i)
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

                    std::cout << "\telapsed time for cooleytukey::iterative::vec::fft with " << (1 << order) << " elements: "; 
                    sw.start(); 
                    std::vector<std::complex<double>> twiddles = cooleytukey::iterative::vec::fft(values);
                    sw.finish(); 
                    sw.print_duration_in_milliseconds(); 
                    std::cout << std::endl; 


                    sw.reset(); 

                    std::cout << "\telapsed time for cooleytukey::iterative::vec::ifft with " << (1 << order) << " elements: "; 
                    sw.start(); 
                    std::vector<std::complex<double>> results = cooleytukey::iterative::vec::ifft(twiddles);
                    sw.finish(); 
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
                    delete [] v; 
                    
            }

        }


        namespace farr
        {
            //full array ( operations on malloc )
            std::complex<double> * fft 
            (std::complex<double> const * c, const unsigned N_SIZE)
            { //DISCLAIMER: c is going to be deleted

                if ( IsPowerOfTwo(N_SIZE))
                {
                    unsigned o = std::log2(N_SIZE);
                    //conversion to complex<double> array
                    std::complex<double> * r = new std::complex<double>[N_SIZE]; 
                    for ( unsigned i = 0; i < N_SIZE; i++)
                    {
                        r[i] = c[bitreverse(i, o)]; //fed with bitinverse index
                        //r.push_back(std::complex<double>(c[i])); //fed with bitinverse index
                    }
                    // for ( unsigned i = 0; i < r.size(); i++)
                    //     std::cout << r[i] << ","; 
                    // std::cout << std::endl;

                    delete [] c;

                    //perform iterative fft
                    unsigned Partitionsize = 2; //2
                    for (unsigned i = 0; i < o; i++)
                    {
                        const std::complex<double> ur = k_nth_root(1,Partitionsize);
                        const unsigned Partitionsize_half = Partitionsize >> 1; 
                        for (unsigned Partition = 0; Partition<N_SIZE; Partition+= Partitionsize )
                            {
                                std::complex<double> multiplier(1,0); 
                                //extract uneven and even
                                for(unsigned k=0; k<Partitionsize_half; k++)  //2
                                {

                                    const std::complex<double> t = multiplier * r[Partition + k + Partitionsize_half];
                                    const std::complex<double> u = r[Partition + k];
                                    r[k+Partition] = u + t;  
                                    r[k+Partition+Partitionsize_half] = u - t; 
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

        std::complex<double> * ifft
                (std::complex<double> const * c, const unsigned N_SIZE)
                { //DISCLAIMER: c is going to be deleted
                    if ( IsPowerOfTwo(N_SIZE))
                    {
                        unsigned o = std::log2(N_SIZE); //calculate the order
                        //conversion to complex<double> array
                        std::complex<double> * r = new std::complex<double>[N_SIZE];  
                        for ( unsigned i = 0; i < N_SIZE; ++i)
                        {
                            r[i] = c[bitreverse(i, o)]; //fed with bitinverse index
                        }
                                         //perform iterative fft
                        delete [] c;
                        unsigned Partitionsize = 2; //2
                        for (unsigned i = 0; i < o; i++)
                        {

                            const std::complex<double> ur = inverse_k_nth_root(1,Partitionsize);
                            const unsigned Paritionsize_half = Partitionsize >> 1; 
                            for (unsigned Partition = 0; Partition<N_SIZE; Partition+= Partitionsize )
                                    {
                                        std::complex<double> multiplier(1,0);  //twiddle
                                        //extract uneven and even
                                        for(unsigned k=0; k<(Paritionsize_half); ++k)  //2
                                        {

                                            const std::complex<double> t = multiplier * r[Partition + k + (Paritionsize_half)];
                                            const std::complex<double> u = r[Partition + k];
                                            r[k+Partition] = u + t;  
                                            r[k+Partition+(Paritionsize_half)] = u - t; 
                                            multiplier *= ur; 
                                            
                                        }
                                    }
                                                                        
                            Partitionsize<<=1;
                        }


                        const std::complex<double> inv = 1. / N_SIZE; 
                        for (unsigned i = 0; i < N_SIZE; i++)
                        {
                            r[i] *= inv; 
                        }
                        return r;
                    }
                    else {
                        throw WRONGSIZE; 
                    }
                    

            }

#include <cassert>
            void benchmark (const unsigned order,
            
                    const bool test) 
            {

                    std::cout << std::endl; 


                    std::complex<double> *v = new std::complex<double>[1 << order];
                    std::complex<double> *original = new std::complex<double>[1 << order];
                    for ( unsigned i = 0; i < 1 << order; i++)
                    {
                        v[i] = (std::rand()%10000);
                        original[i] = v[i];
                    }
                    //std::copy(&v[0], &v[1<<order], back_inserter(values));
                    std::cout << "benchmark for ..iterative::farr::fft: " << std::endl;


                    Stopwatch sw; 
                    sw.reset(); 
                    std::cout << "\telapsed time for cooleytukey::iterative::farr::fft with " << (1 << order) << " elements: "; 
                    sw.start(); 
										std::complex<double> * tmp = v; 
                    v = cooleytukey::iterative::farr::fft(v, 1 << order);
										delete [] tmp; 

                    sw.finish(); 
                    sw.print_duration_in_milliseconds(); 
                    std::cout << std::endl; 


                    sw.reset(); 

                    std::cout << "\telapsed time for cooleytukey::iterative::farr::ifft with " << (1 << order) << " elements: "; 
                    sw.start(); 
                    v = cooleytukey::iterative::farr::ifft(v, 1 << order);
                    sw.finish(); 
                    sw.print_duration_in_milliseconds(); 
                    std::cout << std::endl; 

                    std::cout << "\tcomparing some elements: " << std::endl << "\t\toriginal sample: "; 

                    for ( unsigned i = 0; i < 4; ++i)
                    {
                        std::cout << original[i] << ",";
                    }
                    std::cout << std::endl; 
                    
                    std::cout << "\tcomparing some elements: " << std::endl << "\t\tifft sample: "; 

                    for ( unsigned i = 0; i < 4; ++i)
                    {
                        std::cout << v[i] << ",";
                    }


                    double error_real = 0, error_imag = 0; 
                    std::cout << std::endl; 
                    if ( test)
                    {   
                        for(int64_t i = 0; i<(1<<order); i++)
                        {
                            assert(original[i].real() == v[i].real());
                            error_real += std::abs(original[i].real() - v[i].real());
                            error_imag += std::abs(original[i].imag() - v[i].imag());
                        }

                    }

                    delete[] v;
                    delete[] original; 
                    
                    std::cout << "terminated with real error: " << error_real << " and imaginary residues: " << error_imag << std::endl; 
            }


            std::complex<double> * const fft_convolution(const unsigned order, 
                    std::complex<double> * A, std::complex<double> * B )
                    { //calculates A * B and stores in A.
                        std::complex<double> * A_ = cooleytukey::iterative::farr::fft(A, 1 << order);
                        std::complex<double> * B_ = cooleytukey::iterative::farr::fft(B, 1 << order);
                        for ( unsigned int i = 0; i < (1 << order); i++)
                        {
                            A_[i] *= B_[i];
                        }
                        delete [] B_; 
                        return cooleytukey::iterative::farr::ifft(A_, 1 << order);
                    }

            std::complex<double> * const fft_convolution_without_overflow(const unsigned order, 
                    std::complex<double> * A, std::complex<double> * B )
                    { //calculates A * B and stores in A. DISCLAIMER required a (2**(order + 1)) lot more space.
                        
                        std::complex<double> * A_ = new std::complex<double>[1 << (order + 1)];
                        std::copy(A, A + (1<<order), A_);
                        delete [] A; 
                        A_ = cooleytukey::iterative::farr::fft(A_, 1 << order); //A_ gets deleted and reassigned.

                        std::complex<double> * B_ = new std::complex<double>[1 << (order + 1)];
                        std::copy(B, B + (1<<order), B_);
                        delete [] B; 
                        B_ = cooleytukey::iterative::farr::fft(B_, 1 << order); //B_ gets deleted and reassigned.
                        for ( unsigned int i = 0; i < (1 << order); i++)
                        {
                            A_[i] *= B_[i];
                        }
                        delete [] B_; 
                        //A = cooleytukey::iterative::farr::ifft(A_, 1 << order); //A_ gets deleted
                        return cooleytukey::iterative::farr::ifft(A_, 1 << order); //A_ gets deleted

                    }

            void printasfullnumber (std::complex<double>*
                value,  unsigned order , unsigned base = 10
            ) 
            {
                double res = 0; 
                for (unsigned i = 0; i < (1 << order); i++)
                {
                    double product = round(std::round(value[i].real()) * std::pow(base, i));
                    res += (long) product; 

                }
                std::cout << std::endl; 
                std::cout << res << std::endl; 

            }



            void printasfullfloat (std::complex<double>*
                value,  unsigned order , unsigned base = 10
            ) 
            {
                double res = 0; 
                int exponent = -(1<<order - 1); 
                for (unsigned i = 0; i < (1 << order ); i++)
                {
                    double product = round(std::round(value[i].real()) * std::pow(base, exponent));
                    res += product; 
                    ++exponent;

                }
                std::cout << std::endl; 
                std::cout << res << std::endl; 

            }

            uint64_t approxNumOperationsNewtonRaphsonDivisionFFT(uint64_t dimension, unsigned costOfSimpleMultiplication = 1,
                            unsigned costOfSimpleAddition = 1)
            {
                uint64_t convolution_steps = log2(dimension) * dimension * costOfSimpleMultiplication * costOfSimpleAddition;
                return log2( dimension ) * (convolution_steps + convolution_steps + dimension) ;
            }
            uint64_t approxNumOperationsFFTConvolution(uint64_t dimension, unsigned costOfSimpleMultiplication = 1,
                            unsigned costOfSimpleAddition = 1)
            {
                uint64_t convolution_steps = log2(dimension) * dimension * costOfSimpleAddition * costOfSimpleMultiplication;
                return convolution_steps;
            }
            uint64_t approxNumOperationsFFTConvolutionInC(uint64_t dimension, unsigned costOfSimpleMultiplication = 1,
                            unsigned costOfSimpleAddition = 1)
            {
                uint64_t convolution_steps = log2(dimension * log2(dimension)) * dimension * costOfSimpleAddition * costOfSimpleMultiplication;
                return convolution_steps;
            }

            uint64_t approxNumOperationOfSchoolmath(uint64_t dimension, unsigned costOfSimpleMultiplication = 1,
                            unsigned costOfSimpleAddition = 1)
                            {

                                //approximates the cost of stereotypical schoolmath multiplication and division. which is O(n**2)
                                return dimension * dimension * costOfSimpleAddition * costOfSimpleMultiplication; 
                            }

            void benchmarkfftconvolution(unsigned order, bool verbose)
            {
                srand(time(NULL));
                std::complex<double> * A = new std::complex<double>[1 << order];
                std::complex<double> * B = new std::complex<double>[1 << order];
                
                for ( unsigned i = 0; i<(1 << order-1); i++)
                {
                    A[i] = rand() % 2; 
                    B[i] = rand() % 2; 
                }
                for ( unsigned i = (1 << order-1); i<(1 << order); i++)
                {
                    A[i] = 0; 
                    B[i] = 0; 
                }
                if ( verbose )
                {
                    std::cout << "print some values of A:" << std::endl << "\t "; 
                    for ( unsigned i = 0; i < ((1 << order) > 10 ? 10 : (1 << order)); i++)
                        std::cout << A[i] << ", ";
                    std::cout << std::endl; 
                    std::cout << "print some values of B:" << std::endl << "\t "; 
                    for ( unsigned i = 0; i < ((1 << order) > 10 ? 10 : (1 << order)); i++)
                        std::cout << B[i] << ", ";
                    std::cout << std::endl; 

                }

                Stopwatch sw; 
                sw.reset();
                sw.start();
                A = fft_convolution(order, A, B);
                sw.finish();


                std::cout << "Approx number of operations: " << approxNumOperationsFFTConvolutionInC((1 << order)); 
                std::cout << " ..vs stereotypical multiplication: " << approxNumOperationOfSchoolmath((1<<order)) << std::endl; 
                std::cout << "Time ellapsed for fftconvolution: (sizeof: " << (1 << order) << ")" << std::endl << "\t";
                sw.print_duration_in_microseconds(); 

                std::cout << "print some values of the convolution:" << std::endl << "\t "; 
                for ( unsigned i = 0; i < ((1 << order) > 10 ? 10 : (1 << order)); i++)
                {
                     std::cout << A[i] << ", ";
                }
                //printasfullnumber(A, order, 2);

								delete [] A; 

            }



    }}
}
