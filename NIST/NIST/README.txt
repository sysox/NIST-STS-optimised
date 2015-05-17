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
improved variant/variants of tests. You can also compare performance of 
both implementations and correctness of new implementation.

Each test can be used in two versions v1, v2 (for instance Frequency_v1 
and Frequency_v2),  where v1 (Frequency_v1) defines original implementation
of the test and v2 (Frequency_v2) defines new implementation. Using main-speed.c 
you can compare their speed. Using main-correctness.c you can check 
(compare results) correctness of new version. 
 
You can add new implementation of the arbitrary test as follows: 
	1.Add the new function (Frequency_new(int n){...}) to file frequency.h.
	2.Add its prototype (Frequency_new(int n);) to stat-fncs.h
	3.Change #define Frequency_v2 Frequency3 
				to #define Frequency_v2 Frequency_new
 
*************************************************************************
*                                                                       *
*                           Randomness testing                          *
*                                                                       *
*************************************************************************
0.Create empty in Microsoft Visual Studio.
1. Add files (without main.c, main-speed.c, main-correctness.c) from directories src and include to project.
2. In config.h set the following:
    //#define VERIFY_RESULTS 1
    //#define SPEED 1
     #define FILE_OUTPUT 1
3. Add libfftw3-3.lib to project. 
4. Copy libfftw3-3.dll to .exe.
5. Usage: assess <stream length> [-fast] 
   <stream length> is the length of the individual bit stream(s) to be processed
   -fast           use the faster alternative implementation (version 2) of tests
   (set Command arguments: 10000 -fast if you want to test 10 000 bits with faster version of NIST tests)

*************************************************************************
*                                                                       *
*                           Performance testing                         *
*                                                                       *
*************************************************************************
0.Create empty in Microsoft Visual Studio.
1. Add files (without main.c, main-speed.c, main-correctness.c) from directories src and include to project.
2. In config.h set the following:
    //#define VERIFY_RESULTS 1
    //#define SPEED 1
     #define FILE_OUTPUT 1
3. Add libfftw3-3.lib to project. 
4. Copy libfftw3-3.dll to .exe .
5. Add main-speed.c to project.
6. Usage: main-speed.c scale repeat test_from test_to (for instance 0 10 1 2)
	scale      - 0 (one fixed bit size = 20MB) or 1 (bitsizes are increases in steps)
	repeat     - is number of times( 10) each test is executed (minimum time is taken as the result)
	test_from test_to  - which tests are measured ( 1 - 2 = Freqency and BlockFrequency) 

You can also compare speed of your implementation with some our implemenation.
If you want to compare speed of our Frequency2 and your Frequency_your function
it suffices to change defines to:
    1. #define Frequency_v1 Frequency2
    2. #define Frequency_v2 Frequency_your

*************************************************************************
*                                                                       *
*                           Correctness testing                         *
*                                                                       *
*************************************************************************
0.Create empty in Microsoft Visual Studio.
1. Add files (without main.c, main-speed.c, main-correctness.c) from directories src and include to project.
2. In config.h set the following:
    #define VERIFY_RESULTS 1
    //#define SPEED 1
    //#define FILE_OUTPUT 1
3. Add libfftw3-3.lib to project. 
4. Copy libfftw3-3.dll to .exe .
5. Add main-correctness.c to project.
6. Usage: main-speed.c scale repeat test_from test_to (for instance 0 10 1 2)
	scale      - 0 (one fixed bit size = 20MB) or 1 (bitsizes are increases in steps)
	repeat     - is number of times(10) each test is executed (minimum time is taken as the result)
	test_from test_to  - which tests are measured ( 1 - 2 = Frequency and BlockFrequency) 




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

 