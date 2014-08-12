#include <stdio.h>
#include <math.h>
#include <string.h>
#include "../include/externs.h"
#include "../include/cephes.h"

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
                    B L O C K  F R E Q U E N C Y  T E S T
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

void
BlockFrequency(int M, int n)
{
	
	int		i, j, N, blockSum;
	double	p_value, sum, pi, v, chi_squared;
	
	

	N = n/M; 		/* # OF SUBSTRING BLOCKS      */
	sum = 0.0;
	
	for ( i=0; i<N; i++ ) {
		blockSum = 0;
		for ( j=0; j<M; j++ )
			blockSum += epsilon[j+i*M];
		
		pi = (double)blockSum/(double)M;
		v = pi - 0.5;
		sum += v*v;
		//printf("v*v: %lf\n",v*v);
		// added for printage printf("%d ",blockSum);
	}
	chi_squared = 4.0 * M * sum;
	p_value = cephes_igamc(N/2.0, chi_squared/2.0);
#ifdef VERIFY_RESULTS
	R1.blockfrequency.chi_squared=chi_squared;
	R1.blockfrequency.p_value=p_value;	
#endif
	//printf("%lf",sum);
#ifdef FILE_OUTPUT
	fprintf(stats[TEST_BLOCK_FREQUENCY], "\t\t\tBLOCK FREQUENCY TEST\n");
	fprintf(stats[TEST_BLOCK_FREQUENCY], "\t\t---------------------------------------------\n");
	fprintf(stats[TEST_BLOCK_FREQUENCY], "\t\tCOMPUTATIONAL INFORMATION:\n");
	fprintf(stats[TEST_BLOCK_FREQUENCY], "\t\t---------------------------------------------\n");
	fprintf(stats[TEST_BLOCK_FREQUENCY], "\t\t(a) Chi^2           = %f\n", chi_squared);
	fprintf(stats[TEST_BLOCK_FREQUENCY], "\t\t(b) # of substrings = %d\n", N);
	fprintf(stats[TEST_BLOCK_FREQUENCY], "\t\t(c) block length    = %d\n", M);
	fprintf(stats[TEST_BLOCK_FREQUENCY], "\t\t(d) Note: %d bits were discarded.\n", n % M);
	fprintf(stats[TEST_BLOCK_FREQUENCY], "\t\t---------------------------------------------\n");

	fprintf(stats[TEST_BLOCK_FREQUENCY], "%s\t\tp_value = %f\n\n", p_value < ALPHA ? "FAILURE" : "SUCCESS", p_value); fflush(stats[TEST_BLOCK_FREQUENCY]);
	fprintf(results[TEST_BLOCK_FREQUENCY], "%f\n", p_value); fflush(results[TEST_BLOCK_FREQUENCY]);
#endif
}


// works with global variables: array, R2
// supported rangle of n: 1 - limited by size of int
// supported range of M: 8 - limited by size of int
// not dependent on little/big endian
// bits in array filled from the lowest bits towards the highest bits

void
BlockFrequency2(int M, int n)
{
	/*
	short int LU_byte_weight[256] = 
	{0,1,1,2,1,2,2,3,1,2,2,3,2,3,3,4,1,2,2,3,2,3,3,4,2,3,3,4,3,4,4,5
	,1,2,2,3,2,3,3,4,2,3,3,4,3,4,4,5,2,3,3,4,3,4,4,5,3,4,4,5,4,5,5,6
	,1,2,2,3,2,3,3,4,2,3,3,4,3,4,4,5,2,3,3,4,3,4,4,5,3,4,4,5,4,5,5,6
	,2,3,3,4,3,4,4,5,3,4,4,5,4,5,5,6,3,4,4,5,4,5,5,6,4,5,5,6,5,6,6,7
	,1,2,2,3,2,3,3,4,2,3,3,4,3,4,4,5,2,3,3,4,3,4,4,5,3,4,4,5,4,5,5,6
	,2,3,3,4,3,4,4,5,3,4,4,5,4,5,5,6,3,4,4,5,4,5,5,6,4,5,5,6,5,6,6,7
	,2,3,3,4,3,4,4,5,3,4,4,5,4,5,5,6,3,4,4,5,4,5,5,6,4,5,5,6,5,6,6,7
	,3,4,4,5,4,5,5,6,4,5,5,6,5,6,6,7,4,5,5,6,5,6,6,7,5,6,6,7,6,7,7,8};
	*/
	int		i, N, blockSum,NM,counter;
	double	p_value, sum, pi, v, chi_squared;
	

	N = n/M; 		/* # OF SUBSTRING BLOCKS      */
	sum = 0.0;
	blockSum = 0;
	NM=(N*M / 8);
	if(NM*8<N*M)NM++; // new
	counter = 0;

	for ( i=0; i < NM /*+ 1*/; i++ ) { /// removed +1
		
		blockSum = blockSum + LU_byte_weight[array[i]];
		counter = counter + 8;
		if(counter >= M)
		{
			counter = counter - M ;
			blockSum -= LU_byte_weight[array[i] >> (8 - counter)];
			pi = (double)blockSum/(double)M;
			v = pi - 0.5;
			sum = sum + v*v;
			//printf("v*v: %lf\n",v*v);
			//printf("%d ",blockSum);
			blockSum = LU_byte_weight[array[i] >> (8 - counter)];
		}	
	}
	chi_squared = 4.0 * M * sum;
	p_value = cephes_igamc(N/2.0, chi_squared/2.0);
#ifdef VERIFY_RESULTS
	R2.blockfrequency.chi_squared=chi_squared;
	R2.blockfrequency.p_value=p_value;
#endif
	//printf("%lf ",sum);

#ifdef FILE_OUTPUT
	fprintf(stats[TEST_BLOCK_FREQUENCY], "\t\t\tBLOCK FREQUENCY TEST\n");
	fprintf(stats[TEST_BLOCK_FREQUENCY], "\t\t---------------------------------------------\n");
	fprintf(stats[TEST_BLOCK_FREQUENCY], "\t\tCOMPUTATIONAL INFORMATION:\n");
	fprintf(stats[TEST_BLOCK_FREQUENCY], "\t\t---------------------------------------------\n");
	fprintf(stats[TEST_BLOCK_FREQUENCY], "\t\t(a) Chi^2           = %f\n", chi_squared);
	fprintf(stats[TEST_BLOCK_FREQUENCY], "\t\t(b) # of substrings = %d\n", N);
	fprintf(stats[TEST_BLOCK_FREQUENCY], "\t\t(c) block length    = %d\n", M);
	fprintf(stats[TEST_BLOCK_FREQUENCY], "\t\t(d) Note: %d bits were discarded.\n", n % M);
	fprintf(stats[TEST_BLOCK_FREQUENCY], "\t\t---------------------------------------------\n");

	fprintf(stats[TEST_BLOCK_FREQUENCY], "%s\t\tp_value = %f\n\n", p_value < ALPHA ? "FAILURE" : "SUCCESS", p_value); fflush(stats[TEST_BLOCK_FREQUENCY]);
	fprintf(results[TEST_BLOCK_FREQUENCY], "%f\n", p_value); fflush(results[TEST_BLOCK_FREQUENCY]);
#endif
}
