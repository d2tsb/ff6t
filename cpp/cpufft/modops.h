#include "defines.h"
/* operations in modular space


*/

#define numtype int
#define NO_SUCH_NUMBER_EXCEPTION 49

numtype mod(int a, unsigned m)
{
    return a % m;
}
numtype powmod_recursive(const int basis, unsigned exponent, const unsigned m)
{
    numtype res;
    if (exponent == 0)
    {
        res = 1;
    }
    else if (exponent == 1)
    {

        res = basis;
    }
    else
    {
        res = powmod_recursive(basis, exponent >> 1, m);
        res = (res * res) % m;
        if (exponent % 2)
            res = res * basis % m;
    }
    return res;
}
numtype iterativePowMod(numtype x, numtype y, numtype p)
{

    // Initialize answer
    unsigned long res = 1;

    // Check till the number becomes zero
    while (y > 0)
    {

        // If y is odd, multiply x with result
        if (y % 2 == 1)
            res = (res * x);

        // y = y/2
        y = y >> 1;

        // Change x to x^2
        x = (x * x) % p;
    }
    return (numtype) res;
}
void powmod_recurcive_test()
{
    std::cout << "Powmod 7**2%50: " << std::to_string(powmod_recursive(7, 2, 50)) << std::endl;
    std::cout << "Powmod 7**2%49: " << std::to_string(powmod_recursive(7, 2, 49)) << std::endl;
    std::cout << "Powmod 8**3%50: " << std::to_string(powmod_recursive(8, 3, 50)) << std::endl;
    std::cout << "Powmod 700**700%11: " << std::to_string(powmod_recursive(700, 700, 11)) << std::endl;
}
numtype gcd(numtype a, numtype b)
{

    if (b == 0)
    {
        return a;
    }
    else
    {
        return gcd(b, a % b);
    }
}
#define NOINVERSE (numtype) - 90
numtype modInverseWithPowModIT(numtype A, numtype M)
{
    if (gcd(A, M) != 1)
        return NOINVERSE;
    return iterativePowMod(A, M - 2, M);
}
numtype powerModrecursive(numtype x, unsigned numtype y, unsigned numtype M);
numtype modInverseWithPowMod(numtype A, numtype M)
{ // mod inverse with power
    numtype g = gcd(A, M);
    if (g != 1)
        return NOINVERSE; // actually it should throw.
    else
    {
        // If a and m are relatively prime (gcd(a,m) == 1), then modulo
        // inverse is a^(m-2) mode m
        return powerModrecursive(A, M - 2, M);
    }
}
numtype modInverseWPMOD_WOT(numtype A, numtype M)
{ // mod inverse with power
    // without any test
    //  If a and m are relatively prime (gcd(a,m) == 1), then modulo
    //  inverse is a^(m-2) mode m
    return powerModrecursive(A, M - 2, M);
}
numtype modInverseWPMOD_WOT_IT(numtype A, numtype M)
{ // mod inverse with power
    // without any test
    //  If a and m are relatively prime (gcd(a,m) == 1), then modulo
    //  inverse is a^(m-2) mode m
    return iterativePowMod(A, M - 2, M);
}
// To comput

// To compute x^y under modulo m
numtype powerModrecursive(numtype x, unsigned numtype y, unsigned numtype M)
{
    if (y == 0)
        return 1;

    numtype p = powerModrecursive(x, y / 2, M) % M;
    p = (p * p) % M;

    return (y % 2 == 0) ? p : (x * p) % M;
}
int modInverseIT(int A, int M)
{
    if (gcd(A, M) != 1)
        return NOINVERSE;
    int m0 = M;
    int y = 0, x = 1;

    if (M == 1)
        return 0;

    while (A > 1)
    {
        // q is quotient
        int q = A / M;
        int t = M;

        // m is remainder now, process same as
        // Euclid's algo
        M = A % M, A = t;
        t = y;

        // Update y and x
        y = x - q * y;
        x = t;
    }

    // Make x positive
    if (x < 0)
        x += m0;

    return x;
}
int modInverseWithoutTest(int A, int M)
{
    int m0 = M;
    int y = 0, x = 1;

    if (M == 1)
        return 0;

    while (A > 1)
    {
        // q is quotient
        int q = A / M;
        int t = M;

        // m is remainder now, process same as
        // Euclid's algo
        M = A % M, A = t;
        t = y;

        // Update y and x
        y = x - q * y;
        x = t;
    }

    // Make x positive
    if (x < 0)
        x += m0;

    return x;
}
numtype modInverse(numtype A, numtype M)
{
    for (int X = 1; X < M; X++)
        if (((A % M) * (X % M)) % M == 1)
            return X;
    return NOINVERSE;
}
static void modInverseRace(const int ITERATIONS_RACE, const int);
void testmodinverse(const int iterations = 100000)
{

    std::cout << modInverse(2, 8) << std::endl;
    std::cout << modInverseIT(2, 8) << std::endl;
    std::cout << modInverseIT(12378191, 39916801) << std::endl;

    modInverseRace(iterations, 3000000);
}
numtype modulo2(numtype a)
{
    return (a & 1);
}
numtype itbinarygcd(numtype x, numtype y)
{
    int shl = 0;

    while (x && y && x != y)
    {
        bool eu = !(x & 1);
        bool ev = !(y & 1);

        if (eu && ev)
        {
            ++shl;
            x >>= 1;
            y >>= 1;
        }
        else if (eu && !ev)
            x >>= 1;
        else if (!eu && ev)
            y >>= 1;
        else if (x >= y)
            x = (x - y) >> 1;
        else
        {
            int tmp = x;
            x = (y - x) >> 1;
            y = tmp;
        }
    }

    return !x ? y << shl : x << shl;
}
numtype binary_gcd(numtype x, numtype y)
{

    if (x == y)
        return x;
    if (x == 0)
        return y;
    if (y == 0)
        return x;

    if (~x & 1) // x is even
    {
        if (y & 1) // y is odd
            return gcd(x >> 1, y);
        else
            return gcd(x >> 1, y >> 1) << 1;
    }
    if (~y & 1) // x is odd, y even
        return gcd(x, y >> 1);

    if (x > y)
        return gcd((x - y) >> 1, y);
    return gcd((y - x) >> 1, x);
}
/*
numtype binary_gcd(numtype a, numtype b)
{
    if (a == b)
        return a;

    if (a == 0)

        return b;
    if (b == 0)
        return a;
    const short a_even = ~a & 1; // 0 is even, 1 is odd
    const short b_even = ~b & 1;

    if (a_even && b_even) // both are even
    {
        return binary_gcd(a >> 1, b >> 1) << 1; //  2*bgcd(x/2, y/2);
    }
    if (a_even && !b_even)
    { //  if a is even b is odd
        return binary_gcd(a >> 1, b);
    }
    if (b_even && !a_even)
    { // if a is odd and b even
        return binary_gcd(a, b >> 1);
    }

    // both are odd
    if (a > b)
    {
        return binary_gcd((a - b) >> 1, a);
    }
    else
    {
        return binary_gcd((b - a) >> 1, b);
    }

    // else  {
    //     //if both are odd
    //     return binary_gcd(( a <= b ? b - a : a - b)>>1, a <= b ? b : a);
    return 0;
}
*/
static void modInverseRace(const int ITERATIONS_RACE, const int LIMIT)
{

    const int prime = 39916801;

    // int numbers[ITERATIONS_RACE];
    // int * as = numbers;

    int *as = new int[ITERATIONS_RACE]; // using heap
    for (unsigned i = 0; i < ITERATIONS_RACE; i++)
    {
        as[i] = std::rand() % LIMIT;
    }

    Stopwatch sw;
    // sw.reset();
    // std::cout << "modinverse complexity m ("  << ITERATIONS_RACE << "): ";
    // sw.start();
    // for ( unsigned i = 0; i < ITERATIONS_RACE; i++)
    //     modInverse(as[i], prime);

    //     //modInverse(as[i], prime);
    // //std::cout << "gcd " << a << "," << b << ": " << std::to_string(gcd(a, b)) << std::endl;
    // sw.finish();
    // sw.print_duration_in_microseconds();
    sw.reset();

    std::cout << "modinverse iterative with gcd test (" << ITERATIONS_RACE << "): ";
    sw.start();
    for (unsigned i = 0; i < ITERATIONS_RACE; i++)
        modInverseIT(as[i], prime);
    // std::cout << "binarygcd " << a << "," << b << ": " << std::to_string(binary_gcd(a, b)) << std::endl;
    sw.finish();
    sw.print_duration_in_milliseconds();

    sw.reset();
    std::cout << "modinverse iterative without gcd test (" << ITERATIONS_RACE << "): ";
    sw.start();
    for (unsigned i = 0; i < ITERATIONS_RACE; i++)
        modInverseWithoutTest(as[i], prime);
    // std::cout << "binarygcd " << a << "," << b << ": " << std::to_string(binary_gcd(a, b)) << std::endl;
    sw.finish();
    sw.print_duration_in_milliseconds();

    sw.reset();
    std::cout << "modinverse with recursive powmod and gcd test (" << ITERATIONS_RACE << "): ";
    sw.start();
    for (unsigned i = 0; i < ITERATIONS_RACE; i++)
        modInverseWithPowMod(as[i], prime);
    // std::cout << "binarygcd " << a << "," << b << ": " << std::to_string(binary_gcd(a, b)) << std::endl;
    sw.finish();
    sw.print_duration_in_milliseconds();

    sw.reset();
    std::cout << "modinverse with iterative powmod and gcd test (" << ITERATIONS_RACE << "): ";
    sw.start();
    for (unsigned i = 0; i < ITERATIONS_RACE; i++)
        modInverseWithPowModIT(as[i], prime);
    // std::cout << "binarygcd " << a << "," << b << ": " << std::to_string(binary_gcd(a, b)) << std::endl;
    sw.finish();
    sw.print_duration_in_milliseconds();

    sw.reset();
    std::cout << "modinverse with recursive powmod, without gcd test (" << ITERATIONS_RACE << "): ";
    sw.start();
    for (unsigned i = 0; i < ITERATIONS_RACE; i++)
        modInverseWPMOD_WOT(as[i], prime);
    // std::cout << "binarygcd " << a << "," << b << ": " << std::to_string(binary_gcd(a, b)) << std::endl;
    sw.finish();
    sw.print_duration_in_milliseconds();

    sw.reset();
    std::cout << "modinverse with iterative powmod, without gcd test (" << ITERATIONS_RACE << "): ";
    sw.start();
    for (unsigned i = 0; i < ITERATIONS_RACE; i++)
        modInverseWPMOD_WOT_IT(as[i], prime);
    // std::cout << "binarygcd " << a << "," << b << ": " << std::to_string(binary_gcd(a, b)) << std::endl;
    sw.finish();
    sw.print_duration_in_milliseconds();
}
static void gcdrace(const int ITERATIONS_RACE = 100000)
{

    // int as[ITERATIONS_RACE];
    // int mods[ITERATIONS_RACE];
    int *as = new numtype[ITERATIONS_RACE]; // using heap
    int *mods = new numtype[ITERATIONS_RACE];

    for (unsigned i = 0; i < ITERATIONS_RACE; i++)
    {
        as[i] = std::rand() % 300000000;
        mods[i] = std::rand() % 300000000;
    }

    Stopwatch sw;
    sw.reset();
    std::cout << "gcd with modulo (" << ITERATIONS_RACE << "): ";
    sw.start();
    for (unsigned i = 0; i < ITERATIONS_RACE; i++)
        gcd(as[i], mods[i]);
    // std::cout << "gcd " << a << "," << b << ": " << std::to_string(gcd(a, b)) << std::endl;
    sw.finish();
    sw.print_duration_in_milliseconds();
    sw.reset();

    std::cout << "bin gcd with recursion (" << ITERATIONS_RACE << "): ";
    sw.start();
    for (unsigned i = 0; i < ITERATIONS_RACE; i++)
        binary_gcd(as[i], mods[i]);
    // std::cout << "binarygcd " << a << "," << b << ": " << std::to_string(binary_gcd(a, b)) << std::endl;
    sw.finish();
    sw.print_duration_in_milliseconds();

    sw.reset();
    std::cout << "iterative bin gcd (" << ITERATIONS_RACE << "): ";
    sw.start();
    for (unsigned i = 0; i < ITERATIONS_RACE; i++)
        itbinarygcd(as[i], mods[i]);
    // std::cout << "binarygcd " << a << "," << b << ": " << std::to_string(binary_gcd(a, b)) << std::endl;
    sw.finish();
    sw.print_duration_in_milliseconds();
}
void testGcdAndBinaryGcd()
{

    // std::cout << "gcd 6,4: " << std::to_string(gcd(6, 4)) << std::endl;
    // std::cout << "binarygcd 6,4: " << std::to_string(binary_gcd(6, 4)) << std::endl;
    // std::cout << "gcd 7000,800: " << std::to_string(gcd(7000, 800)) << std::endl;
    // std::cout << "binarygcd 7000,800: " << std::to_string(binary_gcd(7000, 800)) << std::endl;
    // std::cout << "gcd 7000000,800: " << std::to_string(gcd(7000000, 800)) << std::endl;
    // std::cout << "binarygcd 7000000,800: " << std::to_string(binary_gcd(7000000, 800)) << std::endl;
    // std::cout << "iterative binarygcd 7000000,800: " << std::to_string(itbinarygcd(7000000, 800)) << std::endl;
    // gcdrace(7123456, 800);

    gcdrace();
}
/*
numtype binary_gcd(numtype x, numtype y)
{

    if ( x == y )
        return x;
    if ( x == 0 )
        return y;
    if ( y == 0 )
        return x;

    if ( ~x & 1 ) // x is even
    {
        if (y & 1) // y is odd
            return gcd ( x >> 1, y);
        else
            return gcd(x>>1, y>>1) << 1;

    }
    if ( ~y & 1 ) // x is odd, y even
        return gcd(x,y >> 1);

    if ( x > y )
        return gcd((x-y) >> 1, y);
    return gcd((y - x) >> 1, x);
}
*/
numtype calculate_primitventhrootofunity(numtype n, numtype modulo)
{                                   // modulo has to be > 1      gcd ( n, modulo ) have to be atleast 1




    numtype d = gcd(n, modulo - 1); // number of possible solutions for roots of unity in moduloring
    if (d == 1)
        throw NO_SUCH_NUMBER_EXCEPTION;

    numtype solution = iterativePowMod(n, (modulo - 1) * modInverseIT(d, modulo), 213);
    return solution;
}






namespace spacepow2_k
{
    uint32_t cyclicshftl32 (uint32_t value, unsigned int count) {
        const unsigned int mask = CHAR_BIT * sizeof(value) - 1;
        count &= mask;
        return (value << count) | (value >> (-count & mask));
    }

    uint32_t cyclicshftr32 (uint32_t value, unsigned int count) {
        const unsigned int mask = CHAR_BIT * sizeof(value) - 1;
        count &= mask;
        return (value >> count) | (value << (-count & mask));
    }


    bool check_unit_root_property_in_2_space(int lot = 0, unsigned order = 8) 
    {
        unsigned zeuge = 0; 
        unsigned long modspace = (1 << order) + lot; 
        for ( unsigned l = 1; l < modspace; l++)
        {
            for (unsigned k = 0; k < modspace; k++)
            {
                zeuge += iterativePowMod(2, k*l, modspace);
                zeuge %= modspace; 
            }

        }
        return zeuge == 0;
    }
    void executetestofurtest() 
    {
        std::cout << "execute ur == 2 test for pow(2,8) + 1 mod space: " << std::endl; 
        std::cout << std::endl; 

        std::cout << "\tresult: " << (check_unit_root_property_in_2_space(1, 8) ? "true" : "false") << std::endl;  
        std::cout << "execute ur == 2 test for pow(2,8) - 1 mod space: " << std::endl; 
        std::cout << std::endl; 
        std::cout << "\tresult: " << (check_unit_root_property_in_2_space(-1, 8) ? "true" : "false") << std::endl;  



    }

     


    typedef unsigned long p2numtype; 
    constexpr p2numtype m2 = 65537<<1;     
    p2numtype m = 65537;     


}