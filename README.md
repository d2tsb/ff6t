# ff6t

#### Description
Multiply numbers of size ~ 2^1600000 in a couple seconds!

Everything around the (fast) fourier transform implemented in C++ - lightweight, for different datastructures.
It's a bunch of headers which implement the all the checked features list bellow.
Note: all the not checked need to be done.


#### Features
                
+ modops ✅
    + modinverse ✅
         + iterative ✅
		 + recursive ✅
    + powmod ✅
         + iterative ✅
		 + recursive ✅
    + gcd ✅
         + binary ✅
		 + direct ✅
    + calculate/find primitive root ✅
	+ calculate unit root ✅
	+ modulo2 ✅
                
+ dct (I-VIII)
+ fft
	+ cooleytukey for different datastructures ✅
  		+ iterative ✅
		+ recursive ✅
	+ winograd
	+ pfft
	+ bluestein
   	+ chirp-z
 
   	+ dct (discrete cosine transform I-VII) as mod fft
   	+ dst (discrete sine transform as mod fft
 
   	+ functional: reindexing and bitinverse function
+ ntt
	+ dft
 		+ dft337 as example ✅
   	+ general dft over 2 related modspace (fermat- or mersenne-prime space)
   	+ cooleytukey in modspace
   		+ iterative
   	 	+ recursive
+ threaded ffts
  	+ for all dfts
  	+ for all ffts
  	+ gpu implementation
+ big number library
+ speedtests ✅
	+ gcdrace ✅
 	+ modInverseRace ✅
  	+ powmodrace ✅
  	  

+ dct (discrete cosine transform I-VII) as mod fft
+ dst (discrete sine transform) as mod fft

+ functional
	+ reindexing
 	+ bitinverse function

