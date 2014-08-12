*************************************************************************
*        Alternative (faster) implementation of the NIST STS tests      *
*                                                                       *
*            NIST STS Version 2.1.1 - Alt Version A                     *
*************************************************************************
*                                                                       *
*                           U S A G E                                   *
*                                                                       *
*************************************************************************

This is an alternative implementation of the NIST statistical randomness
tests. The source codes include both the original and the new alternative
implementation. Compile the source code the same way as you would the
original code. The new binary supports a new optional command line argument
which selects the code being executed: the original code or the new
alternative code.

Usage: assess <stream length> [-fast]
   <stream length> is the length of the individual bit stream(s) to be processed
   -fast           use the faster alternative implementation of tests


*************************************************************************
*                                                                       *
*                           Compilation                                 *
*                                                                       *
*************************************************************************

The source code can be compiled in one of the three modes. The modes can
be selected in config.h Uncoment exactlty one of the following lines:

//#define VERIFY_RESULTS 1
//#define SPEED 1
//#define FILE_OUTPUT 1

Uncommenting the definition of FILE_OUTPUT produces a binary with the
standard functionality. This is what you typically want and this is the
default. The usage of the binary is described in the previous section.

Uncommenting the definition of VERIFY_RESULTS produces a binary that
runs a series of tests which compare the results of the original code
with the results of the new code. The one and only coomand line 
parameter of resulting program is the number of the randomness test
(1-15). Bitsizes and secondary parameters can be modified in the source
code (search for funtion test() in main.c).

Uncommenting the definition of SPEED produces a binary that runs a series
of speed measurements. The commmand line arguments are as follows:
./a.out scale repeat test_from test_to
where scale is 0 (bitsizes are increades in steps) or 1 (one fixed bit size),
repeat is number of times each test is executed (minimum time is taken as the
result), test_from and test_to affect which randomness tests are executed
(1 to 15). Other parameters can be modified in the source code (search for
funtion speed in main.c).

*************************************************************************
*                                                                       *
*              Introduce support for m>25 (and m<=32)                   *
*              in Rank and Overlapping Template Matching tests          *
*                                                                       *
*************************************************************************

If you feel limited with the maximum value of m (m <=25) in the Rank
and Overlapping Template Matching tests then simply replace the calls
of get_nth_block4() with the calls of XXXXXXXXXX() in the relevant
parts of the source code. For the rank test search for function YYY()
in rank.c. For the Overlapping Template Matching test search for
the funtion OverlappingTemplateMatchings2() in 
overlappingTemplateMatchings.c
