
#include "../cpufft/cpufft.h"

void execute_vec_cooleytuckey_fft()
{
    std::cout << "Test of (recursive) cooleytuckey::recursive::vec(tor)::fft/ifft" << std::endl;
    std::srand(std::time(0));
    //std::vector< double > v({ 1, 1 }); 
    //std::vector< double > v; 
    std::vector< double > v({ 1, 2, 3, 4, 5, 6, 7, 8 }); 

    // for ( unsigned i = 0; i < 1 << ORDER; i++)
    // {
    //     v.push_back(std::rand()%10000);
    // }
    
    
    /*

    v.resize(1<< ORDER);
     for (std::vector<double>::iterator it = v.begin(); it != v.end(); ++it)
    {
       *it = rand();
    }
    
    
    */
    //std::cout << fft(v).size(); o
    std::vector<std::complex<double>> cv = cooleytukey::recursive::vec::fft(v);
        for (unsigned i = 0; i < STATEMENT; i++)
                {
                    std::cout << cv[i] << "," ;
                }
        std::cout << "fft done" << std::endl;
    cooleytukey::recursive::vec::ifft(cv);
        
    std::cout << std::endl; 
    
   //std::cout << "Statement: " << STATEMENT; 
}

void execute_arr_cooleytuckey_fft()
{

    std::cout << "Test of (recursive) cooleytuckey::recursive::arr(std::array)::fft/ifft" << std::endl;
    std::srand(std::time(0));
    //std::vector< double > v({ 1, 1 }); 
    std::array< double, 1 << ORDER> v; 
    //std::vector< double > v({ 1, 2, 3, 4, 5, 6, 7, 8 }); 

    for ( unsigned i = 0; i < 1 << ORDER; i++)
    {
        v[i] = (std::rand()%10000);
    }
    
    
    /*

    v.resize(1<< ORDER);
     for (std::vector<double>::iterator it = v.begin(); it != v.end(); ++it)
    {
       *it = rand();
    }
    
    
    */
    //std::cout << fft(v).size(); o
    std::array<std::complex<double>, 1 << ORDER> cv = cooleytukey::recursive::arr::fft(v);
        for (unsigned i = 0; i < STATEMENT; i++)
                {
                    std::cout << cv[i] << "," ;
                }
        std::cout << "fft done" << std::endl;
    cooleytukey::recursive::arr::ifft(cv);
        
    std::cout << "ifft done\n";
    
   //std::cout << "Statement: " << STATEMENT; 
}


void execute_vec_cooleytuckey_fft_iterative()
{
    std::cout << "Test of (iterative) cooleytuckey::iterative::vec(std::vector)::fft/ifft" << std::endl;
    /*
    std::cout << "sizeof(unsinged):" << sizeof(unsigned) << std::endl; 
    std::cout << "bitinverse of 0:" << inverse(0, 4) << std::endl;
    std::cout << int_as_string(inverse(0,4)) << std::endl; 
    std::cout << int_as_string(7) << std::endl; 
    std::cout << int_as_string(8) << std::endl; 
    std::cout << int_as_string(inverse(8,4)) << std::endl; 
    */
    //std::vector<double> v({1,1});
    std::vector<double> v({1,2,3,4,5,6,7,8});
    auto r = cooleytukey::iterative::vec::fft(v);
    for (unsigned i = 0; i < (1 << ORDER); i++) 
    {
        std::cout << r[i] << ","; 
    }


    std::cout << "fft done\n";

    auto inverted = cooleytukey::iterative::vec::ifft(r);

    for (unsigned i = 0; i < (1 << ORDER); i++) 
    {
        std::cout << inverted[i] << ","; 
    }


}

void printasfullnumber (std::vector<std::complex<double>> 
    value
) 
{
    unsigned long res = 0; 
    for (unsigned i = 0; i < value.size(); i++)
    {
        double product = round(std::round(value[i].real()) * std::pow(10, i));
        res += (long) product; 


    }
    std::cout << std::endl; 
    std::cout << res << std::endl; 

}



void execute_convolution_recursive_multiplication_via_fft()
{
   std::cout << "first example" << std::endl;
   printasfullnumber(
   cooleytukey::convolution_recursive::multiplication_via_fft(
    std::vector<double>({  7, 1, 0, 0, 0, 0, 0, 0}),
    std::vector<double>({  7, 1, 0, 0, 0, 0, 0, 0}) //which is 17
   )); 

   std::cout << std::endl;
   std::cout << "second example" << std::endl;
   printasfullnumber(
    cooleytukey::convolution_recursive::multiplication_via_fft(
    std::vector<double>({  4, 1, 0, 0, 0, 0, 0, 0, 
                            0, 0, 0, 0, 0, 0, 0, 0}),
    std::vector<double>({  4, 1, 0, 0, 0, 0, 0, 0,
                            0, 0, 0, 0, 0, 0, 0, 0})
   ));

   std::cout << std::endl;
   std::cout << "third example (20 * 20 == 400)" << std::endl;
   printasfullnumber(
    cooleytukey::convolution_recursive::multiplication_via_fft(
    std::vector<double>({  0, 2, 0, 0, 0, 0, 0, 0}),
    std::vector<double>({  0, 2, 0, 0, 0, 0, 0, 0})
   ));
   std::cout << std::endl;
   std::cout << "forth example: (606*606 == 367236) " << std::endl;
   printasfullnumber(
    cooleytukey::convolution_recursive::multiplication_via_fft(
    std::vector<double>({  6, 0, 6, 0, 0, 0, 0, 0}),
    std::vector<double>({  6, 0, 6, 0, 0, 0, 0, 0})
   ));

    //actually numerically unstable.



}

void check_iterative_to_recursive_crosscompatibiliy () 
{

    std::cout << "Test of iterative to recursive (cross)compatibily" << std::endl;
    std::vector <double> values = { 1, 2, 3, 4, 5, 6, 7, 8};

    std::cout << "\nrecursive -> iterative:" << std::endl; 
    std::cout << "recursive::vec::fft" << std::endl; 
    auto res = cooleytukey::recursive::vec::fft(values); 
    std::cout << "iterative::vec::ifft" << std::endl; 
    auto inverted = cooleytukey::iterative::vec::ifft(res); 
    std::cout << "\n got: "<< std::endl; 
    for (unsigned i = 0; i < inverted.size(); i++)
    {
       std::cout << inverted[i] << ","; 
    }

    std::cout << "\nnow the other way around. iterative -> recursive" << std::endl; 
    std::cout << "iterative::vec::fft" << std::endl; 
    auto res2 = cooleytukey::iterative::vec::fft(values); 
    std::cout << "iterative::vec::ifft" << std::endl; 
    auto inverted2 = cooleytukey::recursive::vec::ifft(res2); 
    std::cout << "\n got: "<< std::endl; 
    for (unsigned i = 0; i < inverted2.size(); i++)
    {
       std::cout << inverted2[i] << ","; 
    }

    std::cout << std::endl; 
    std::cout << "test done." << std::endl; 
}



void execute_dft_ntt337 () 
{


    srand(time(NULL));
    unsigned size = 8;  
    std::vector < unsigned > v; 
    std::cout << std::endl; 
    std::cout << "testing dft::ntt337::ntt(number theortical transform) && intt of size(8): " << std::endl << "\tgot original: "; 
    for ( unsigned i = 0;   i<size;   i++)
    {
       v.push_back(rand() % 10);  //since the the maximum product of the coefficients need to be inside Z/337Z  and 10*10*log(10) is ruffly 200 as maximum product
       //v.push_back(i);
       std::cout << v[i] << ", "; 
    }
    std::cout << std::endl; 
    std::cout << std::endl; 

    v = dft::ntt337::ntt(v);
    std::cout << "\twith coefficients (after ntt): ";
    for ( unsigned i = 0;   i<size;   i++)
    {
       std::cout << v[i] << ", "; 
    }
    
    std::cout << std::endl; 
    std::cout << std::endl; 

    v = dft::ntt337::intt(v);
    std::cout << "\tand inverse (intt): ";
    for ( unsigned i = 0;   i<size;   i++)
    {
       std::cout << v[i] << ", "; 
    }
   
 
    std::cout << std::endl; 
    std::cout << std::endl; 


}

void testinverse(unsigned order = 4)
{


    std::vector<int> ordered = { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 27000};
    std::cout << "print ordered elements:" << " (order of "<< order << ")" << std::endl;
    for ( auto element : ordered)
        std::cout << element << ",";

    std::vector<int> bitinserse_indices;
    for ( int i = 0; i < ordered.size(); i++)
    {
        //being pedantic
        bitinserse_indices.push_back(  cooleytukey::iterative::bitinverse(ordered[i], order));
    }
    std::cout << std::endl; 
    std::cout << "print bitinversed elements:" << std::endl;
    for ( auto element : bitinserse_indices)
        std::cout << element << ",";

    std::cout << std::endl; 

    
}

void testbitreverse(const unsigned order)
{
    std::cout << "testing bitreverse: " << std::endl; 
    for ( unsigned i = 0; i < (1 << order); ++i)
    {
        std::cout << "\t\tgot " << i << " as: " << int_as_string(i) << std::endl; 
    }
    std::cout << "\tnow the corresponding bitreverse:" << std::endl;
    for ( unsigned i = 0; i < (1 << order); ++i)
    {
        std::cout << "\t\tgot " << i << " reversed (" << cooleytukey::iterative::bitreverse(i,order)\
        << ") as: " << int_as_string(cooleytukey::iterative::bitreverse(i, order)) << std::endl; 

    }

}


void execute_dft_modspace2k () 
{


    srand(time(NULL));
    unsigned order = 8; 

    //unsigned size = 1 << (order-1);  which is not prime, so it does no really work
    unsigned size = order*2; //make 2 order*2 nth root of unity in gf(2**n  + 1) 
    //2**8 + 1 is the last fermat prime
    std::vector < unsigned > v; 
    std::cout << std::endl; 
    std::cout << "testing dft::modspace2k::ntt(number theortical transform) && intt of size(8): " << std::endl << "\tgot original: "; 
    for ( unsigned i = 0;   i<size;   i++)
    {
       v.push_back(rand() % 3);  //since the the maximum product of the coefficients need to be inside Z/337Z  and 10*10*log(10) is ruffly 200 as maximum product
       //v.push_back(i);
       std::cout << v[i] << ", "; 
    }
    std::cout << std::endl; 
    std::cout << std::endl; 

    v = dft::modspace2k::ntt(v, order, 1);
    std::cout << "\twith coefficients (after ntt): ";
    for ( unsigned i = 0;   i<size;   i++)
    {
       std::cout << v[i] << ", "; 
    }
    
    std::cout << std::endl; 
    std::cout << std::endl; 

    v = dft::modspace2k::intt(v, order, 1);
    std::cout << "\tand inverse (intt): ";
    for ( unsigned i = 0;   i<size;   i++)
    {
       std::cout << v[i] << ", "; 
    }
   
 
    std::cout << std::endl; 
    std::cout << std::endl; 


}






int main()
{
    // Stopwatch sw; 
    // sw.start();
    //execute_vec_cooleytuckey_fft();
    // sw.finish();
    // sw.print_duration_in_milliseconds(); 

    // sw.reset(); 
    // sw.start();
    // execute_vec_cooleytuckey_fft();
    // sw.finish();
    // sw.print_duration_in_milliseconds(); 
    //testinverse(2);
    
    

    //execute_convolution_recursive_multiplication_via_fft(); 



    //execute_vec_cooleytuckey_fft_iterative();
    //check_iterative_to_recursive_crosscompatibiliy(); 
    //execute_dft_ntt337(); 

    spacepow2_k::executetestofurtest(); 
    execute_dft_modspace2k();  


    //cooleytukey::iterative::farr::benchmark(24);
    //cooleytukey::iterative::vec::benchmark(24);
    //powmod_recurcive_test();
    //testGcdAndBinaryGcd(); 
    //testmodinverse(5000000); 
    //gcdrace(5000000);:w

    
    //testbitreverse(5);
}

