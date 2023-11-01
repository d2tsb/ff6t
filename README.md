# ff6t
everything around the (fast) fourier transform implemented in c++

it's a bunch of headers which implement the following things:


modops
modInverse
    iterative
    recursive
powmod
    iterative
    recursive
gcd
    binary 
    direct 
    ..

calculate/find primitive root 
calculate unit root

modulo2


dct (I-VII)


fft cooleytuckey
        iterative   
        recursive
    winograd
    pfft
    bluestein
    chirp-z

dct (discrete cosine transform I-VII) as mod fft
dst (discrete sine transform) as mod fft

functional: reindexing and bitinverse function



ntt 
    dft
        dft337 as example
        general dft in power of 2 modspace

    cooleytuckey
        iterative
        recursive

        

threaded ffts
    all ntts 
    all cooleytuckeys
    using gpu


big number library


speedtests
    gcdrace 
    modInverseRace  
    powmodrace
