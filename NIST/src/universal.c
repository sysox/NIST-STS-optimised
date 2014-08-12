#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include "../include/erf.h"
#include "../include/externs.h"
#include "../include/utilities.h"
#include "../include/cephes.h"
#include "../include/tools.h"

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
                         U N I V E R S A L  T E S T
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

void
Universal(int n)
{
	int		i, j, p, L, Q, K;
	double	arg, sqrt2, sigma, phi, sum, p_value, c;
	long	*T, decRep;
	double	expected_value[17] = { 0, 0, 0, 0, 0, 0, 5.2177052, 6.1962507, 7.1836656,
				8.1764248, 9.1723243, 10.170032, 11.168765,
				12.168070, 13.167693, 14.167488, 15.167379 };
	double   variance[17] = { 0, 0, 0, 0, 0, 0, 2.954, 3.125, 3.238, 3.311, 3.356, 3.384,
				3.401, 3.410, 3.416, 3.419, 3.421 };
	unsigned int mask;
	/* * * * * * * * * ** * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
	 * THE FOLLOWING REDEFINES L, SHOULD THE CONDITION:     n >= 1010*2^L*L       *
	 * NOT BE MET, FOR THE BLOCK LENGTH L.                                        *
	 * * * * * * * * * * ** * * * * * * * * * * * * * * * * * * * * * * * * * * * */
	L = 5;
	if ( n >= 387840 )     L = 6;
	if ( n >= 904960 )     L = 7;
	if ( n >= 2068480 )    L = 8;
	if ( n >= 4654080 )    L = 9;
	if ( n >= 10342400 )   L = 10;
	if ( n >= 22753280 )   L = 11;
	if ( n >= 49643520 )   L = 12;
	if ( n >= 107560960 )  L = 13;
	if ( n >= 231669760 )  L = 14;
	if ( n >= 496435200 )  L = 15;
	if ( n >= 1059061760 ) L = 16;
	
	Q = 10*(int)pow(2, L);
	K = (int) (floor(n/L) - (double)Q);	 		    /* BLOCKS TO TEST */
	mask = (1 << L)-1;
	p = (int)pow(2, L);
	if ( (L < 6) || (L > 16) || ((double)Q < 10*pow(2, L)) ||
		 ((T = (long *)calloc(p, sizeof(long))) == NULL) ) {
#ifdef FILE_OUTPUT		
		fprintf(stats[TEST_UNIVERSAL], "\t\tUNIVERSAL STATISTICAL TEST\n");
		fprintf(stats[TEST_UNIVERSAL], "\t\t---------------------------------------------\n");
		fprintf(stats[TEST_UNIVERSAL], "\t\tERROR:  L IS OUT OF RANGE.\n");
		fprintf(stats[TEST_UNIVERSAL], "\t\t-OR- :  Q IS LESS THAN %f.\n", 10*pow(2, L));
		fprintf(stats[TEST_UNIVERSAL], "\t\t-OR- :  Unable to allocate T.\n");
#endif
		return;
	}
	
	/* COMPUTE THE EXPECTED:  Formula 16, in Marsaglia's Paper */
	c = 0.7 - 0.8/(double)L + (4 + 32/(double)L)*pow(K, -3/(double)L)/15;
	sigma = c * sqrt(variance[L]/(double)K);
	sqrt2 = sqrt(2);
	sum = 0.0;
	for ( i=0; i<p; i++ )
		T[i] = 0;
	for ( i=1; i<=Q; i++ ) {		/* INITIALIZE TABLE */
		decRep = 0;
		for ( j=0; j<L; j++ )
			decRep += epsilon[(i-1)*L+j] * (long)pow(2, L-j-1);
		T[decRep] = i;
		//if( (mask & get_nth_block4(array,L*(i-1))) != decRep)printf("chyba %d",i-1);
		
	}
	for ( i=Q+1; i<=Q+K; i++ ) { 	/* PROCESS BLOCKS */
		decRep = 0;
		for ( j=0; j<L; j++ )
			decRep += epsilon[(i-1)*L+j] * (long)pow(2, L-j-1);
		sum += log(i - T[decRep])/log(2);
		T[decRep] = i;
	}
	//printf("%lf ",sum);
	phi = (double)(sum/(double)K);

#ifdef FILE_OUTPUT
	fprintf(stats[TEST_UNIVERSAL], "\t\tUNIVERSAL STATISTICAL TEST\n");
	fprintf(stats[TEST_UNIVERSAL], "\t\t--------------------------------------------\n");
	fprintf(stats[TEST_UNIVERSAL], "\t\tCOMPUTATIONAL INFORMATION:\n");
	fprintf(stats[TEST_UNIVERSAL], "\t\t--------------------------------------------\n");
	fprintf(stats[TEST_UNIVERSAL], "\t\t(a) L         = %d\n", L);
	fprintf(stats[TEST_UNIVERSAL], "\t\t(b) Q         = %d\n", Q);
	fprintf(stats[TEST_UNIVERSAL], "\t\t(c) K         = %d\n", K);
	fprintf(stats[TEST_UNIVERSAL], "\t\t(d) sum       = %f\n", sum);
	fprintf(stats[TEST_UNIVERSAL], "\t\t(e) sigma     = %f\n", sigma);
	fprintf(stats[TEST_UNIVERSAL], "\t\t(f) variance  = %f\n", variance[L]);
	fprintf(stats[TEST_UNIVERSAL], "\t\t(g) exp_value = %f\n", expected_value[L]);
	fprintf(stats[TEST_UNIVERSAL], "\t\t(h) phi       = %f\n", phi);
	fprintf(stats[TEST_UNIVERSAL], "\t\t(i) WARNING:  %d bits were discarded.\n", n-(Q+K)*L);
	fprintf(stats[TEST_UNIVERSAL], "\t\t-----------------------------------------\n");
#endif

	arg = fabs(phi-expected_value[L])/(sqrt2 * sigma);
	p_value = erfc(arg);

#ifdef VERIFY_RESULTS
	R1.universal.p_value=p_value;
	R1.universal.sum=sum;
	R1.universal.phi=phi;
#endif
#ifdef FILE_OUTPUT
	if ( isNegative(p_value) || isGreaterThanOne(p_value) )
		fprintf(stats[TEST_UNIVERSAL], "\t\tWARNING:  P_VALUE IS OUT OF RANGE\n");

	fprintf(stats[TEST_UNIVERSAL], "%s\t\tp_value = %f\n\n", p_value < ALPHA ? "FAILURE" : "SUCCESS", p_value); fflush(stats[TEST_UNIVERSAL]);
	fprintf(results[TEST_UNIVERSAL], "%f\n", p_value); fflush(results[TEST_UNIVERSAL]);
#endif
	free(T);
}

void
Universal2(int n)
{
	int		i, p, L, Q, K;
	double	arg, sqrt2, log2,sigma, phi, sum, p_value, c;
	long	*T;
	double	expected_value[17] = { 0, 0, 0, 0, 0, 0, 5.2177052, 6.1962507, 7.1836656,
				8.1764248, 9.1723243, 10.170032, 11.168765,
				12.168070, 13.167693, 14.167488, 15.167379 };
	double   variance[17] = { 0, 0, 0, 0, 0, 0, 2.954, 3.125, 3.238, 3.311, 3.356, 3.384,
				3.401, 3.410, 3.416, 3.419, 3.421 };
	unsigned int window,mask;

	/* * * * * * * * * ** * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
	 * THE FOLLOWING REDEFINES L, SHOULD THE CONDITION:     n >= 1010*2^L*L       *
	 * NOT BE MET, FOR THE BLOCK LENGTH L.                                        *
	 * * * * * * * * * * ** * * * * * * * * * * * * * * * * * * * * * * * * * * * */
	L = 5;
	if ( n >= 387840 )     L = 6;
	if ( n >= 904960 )     L = 7;
	if ( n >= 2068480 )    L = 8;
	if ( n >= 4654080 )    L = 9;
	if ( n >= 10342400 )   L = 10;
	if ( n >= 22753280 )   L = 11;
	if ( n >= 49643520 )   L = 12;
	if ( n >= 107560960 )  L = 13;
	if ( n >= 231669760 )  L = 14;
	if ( n >= 496435200 )  L = 15;
	if ( n >= 1059061760 ) L = 16;
	
	Q = 10*(int)pow(2, L);
	K = (int) (floor(n/L) - (double)Q);	 		    /* BLOCKS TO TEST */
	mask = (1 << L)-1;
	p = (int)pow(2, L);
	if ( (L < 6) || (L > 16) || ((double)Q < 10*pow(2, L)) ||
		 ((T = (long *)calloc(p, sizeof(long))) == NULL) ) {
#ifdef FILE_OUTPUT		
		fprintf(stats[TEST_UNIVERSAL], "\t\tUNIVERSAL STATISTICAL TEST\n");
		fprintf(stats[TEST_UNIVERSAL], "\t\t---------------------------------------------\n");
		fprintf(stats[TEST_UNIVERSAL], "\t\tERROR:  L IS OUT OF RANGE.\n");
		fprintf(stats[TEST_UNIVERSAL], "\t\t-OR- :  Q IS LESS THAN %f.\n", 10*pow(2, L));
		fprintf(stats[TEST_UNIVERSAL], "\t\t-OR- :  Unable to allocate T.\n");
#endif
		return;
	}
	
	/* COMPUTE THE EXPECTED:  Formula 16, in Marsaglia's Paper */
	c = 0.7 - 0.8/(double)L + (4 + 32/(double)L)*pow(K, -3/(double)L)/15;
	sigma = c * sqrt(variance[L]/(double)K);
	sqrt2 = sqrt(2);
	log2 = log(2);
	sum = 0.0;
	for ( i=0; i<p; i++ )
		T[i] = 0;
	for ( i=1; i<=Q; i++ ) {		/* INITIALIZE TABLE */
		window = get_nth_block4(array,L*(i-1));
		T[window & mask] = i;
		 
	}
	for ( i=Q+1; i<=Q+K; i++ ) { 	/* PROCESS BLOCKS */
		window = get_nth_block4(array,L*(i-1));
		sum += log(i - T[window & mask])/log2; //treba sem dat /log(2) a rovnake
		T[window & mask] = i;
	}
	//sum = sum / log(2);//nekoresponduje vdaka zaokruhlovaniu vid vyssie
	//printf("%lf ",sum);
	phi = (double)(sum/(double)K);
	
#ifdef FILE_OUTPUT
	fprintf(stats[TEST_UNIVERSAL], "\t\tUNIVERSAL STATISTICAL TEST\n");
	fprintf(stats[TEST_UNIVERSAL], "\t\t--------------------------------------------\n");
	fprintf(stats[TEST_UNIVERSAL], "\t\tCOMPUTATIONAL INFORMATION:\n");
	fprintf(stats[TEST_UNIVERSAL], "\t\t--------------------------------------------\n");
	fprintf(stats[TEST_UNIVERSAL], "\t\t(a) L         = %d\n", L);
	fprintf(stats[TEST_UNIVERSAL], "\t\t(b) Q         = %d\n", Q);
	fprintf(stats[TEST_UNIVERSAL], "\t\t(c) K         = %d\n", K);
	fprintf(stats[TEST_UNIVERSAL], "\t\t(d) sum       = %f\n", sum);
	fprintf(stats[TEST_UNIVERSAL], "\t\t(e) sigma     = %f\n", sigma);
	fprintf(stats[TEST_UNIVERSAL], "\t\t(f) variance  = %f\n", variance[L]);
	fprintf(stats[TEST_UNIVERSAL], "\t\t(g) exp_value = %f\n", expected_value[L]);
	fprintf(stats[TEST_UNIVERSAL], "\t\t(h) phi       = %f\n", phi);
	fprintf(stats[TEST_UNIVERSAL], "\t\t(i) WARNING:  %d bits were discarded.\n", n-(Q+K)*L);
	fprintf(stats[TEST_UNIVERSAL], "\t\t-----------------------------------------\n");
#endif

	arg = fabs(phi-expected_value[L])/(sqrt2 * sigma);
	p_value = erfc(arg);

#ifdef VERIFY_RESULTS
	R2.universal.p_value=p_value;
	R2.universal.sum=sum;
	R2.universal.phi=phi;
#endif

#ifdef FILE_OUTPUT
	if ( isNegative(p_value) || isGreaterThanOne(p_value) )
		fprintf(stats[TEST_UNIVERSAL], "\t\tWARNING:  P_VALUE IS OUT OF RANGE\n");

	fprintf(stats[TEST_UNIVERSAL], "%s\t\tp_value = %f\n\n", p_value < ALPHA ? "FAILURE" : "SUCCESS", p_value); fflush(stats[TEST_UNIVERSAL]);
	fprintf(results[TEST_UNIVERSAL], "%f\n", p_value); fflush(results[TEST_UNIVERSAL]);
#endif

	free(T);
}
