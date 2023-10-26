
#include "defines.h"
#include "modops.h"


namespace dft
{
    namespace ntt337
    { 

        /**
         * works only with at max 2**n < 337 == 256 elements
        */
        const int m = 337; 
        const int ur = 85; 
        //template <typename T> 
        std::vector <unsigned> ntt(std::vector<unsigned>  v ) 
        {
            std::vector<unsigned> v_temp(v); 
                    //std::copy(&v.data()[0], &v.data()[v.size()], back_inserter(v_temp));
            
            for ( unsigned k = 0;  k < v_temp.size(); k++)
            {
                numtype res = 0; 
                for ( unsigned j = 0;  j < v_temp.size(); j++)
                {
                    res = (res + (v_temp[j] * iterativePowMod(ur, 
                                                        j * k,
                                                        m)) ) %m;

                    std::cout << res << ", "; 
                }
                    std::cout <<"k is :" << k << std::endl; 
                v[k] = res; 
            }
            return std::vector<unsigned> (v);
            //with time complexity n^2
        } 

        std::vector <unsigned> intt(std::vector<unsigned> & v ) 
        {
            const std::vector<unsigned> v_temp(v); 
            //std::copy(&v.data()[0], &v.data()[v.size()], back_inserter(v_temp));

            unsigned length = v_temp.size(); 

            for ( unsigned k = 0;  k < length; k++)
            {

                numtype res = 0; 

                for ( unsigned j = 0;  j < length; j++)
                {

                    res = (res + v_temp[j] * modInverseIT (iterativePowMod(ur,j*k, m), m))
                                                    %m;
                    std::cout << res << ", "; 
                }
                std::cout <<"k is :" << k << std::endl; 
                v[k] = res; 
            }
            std::cout << "now finally multiplying the inverse of length (dimension)" << std::endl;
            //with time complexity n^2
            unsigned m_inverse = modInverseIT(v_temp.size() % m, m );
            for ( unsigned t = 0;  t < v_temp.size(); t++)
            {

                v[t] = (m_inverse * v[t]) % m;
                std::cout <<"temp_inv * n^(-1) is :" << v[t] << std::endl; 

            }

            return std::vector<unsigned> (v);
        } 







    }
   
}