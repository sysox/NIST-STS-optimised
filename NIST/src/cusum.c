#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "../include/externs.h"
#include "../include/cephes.h"
#include "../include/tools.h"  

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
		    C U M U L A T I V E  S U M S  T E S T
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

void
CumulativeSums(int n)
{
	int		S, sup, inf, z, zrev, k;
	double	sum1, sum2, p_value;
	

	S = 0;
	sup = 0;
	inf = 0;
	for ( k=0; k<n; k++ ) {
		epsilon[k] ? S++ : S--;
		if ( S > sup )
			sup++;
		if ( S < inf )
			inf--;
		z = (sup > -inf) ? sup : -inf;
		zrev = (sup-S > S-inf) ? sup-S : S-inf;
	}
	
	//printf("%d %d", z, zrev);
#ifdef VERIFY_RESULTS
	R1.cusum.z=z;
	R1.cusum.zrev=zrev;
#endif

	// forward
	sum1 = 0.0;
	for ( k=(-n/z+1)/4; k<=(n/z-1)/4; k++ ) {
		sum1 += cephes_normal(((4*k+1)*z)/sqrt(n));
		sum1 -= cephes_normal(((4*k-1)*z)/sqrt(n));
	}
	sum2 = 0.0;
	for ( k=(-n/z-3)/4; k<=(n/z-1)/4; k++ ) {
		sum2 += cephes_normal(((4*k+3)*z)/sqrt(n));
		sum2 -= cephes_normal(((4*k+1)*z)/sqrt(n));
	}

	p_value = 1.0 - sum1 + sum2;
#ifdef VERIFY_RESULTS
	R1.cusum.p_valueA=p_value;
	R1.cusum.sum1A=sum1;
	R1.cusum.sum2A=sum2;
#endif
	
#ifdef FILE_OUTPUT
	fprintf(stats[TEST_CUSUM], "\t\t      CUMULATIVE SUMS (FORWARD) TEST\n");
	fprintf(stats[TEST_CUSUM], "\t\t-------------------------------------------\n");
	fprintf(stats[TEST_CUSUM], "\t\tCOMPUTATIONAL INFORMATION:\n");
	fprintf(stats[TEST_CUSUM], "\t\t-------------------------------------------\n");
	fprintf(stats[TEST_CUSUM], "\t\t(a) The maximum partial sum = %d\n", z);
	fprintf(stats[TEST_CUSUM], "\t\t-------------------------------------------\n");

	if ( isNegative(p_value) || isGreaterThanOne(p_value) )
		fprintf(stats[TEST_CUSUM], "\t\tWARNING:  P_VALUE IS OUT OF RANGE\n");

	fprintf(stats[TEST_CUSUM], "%s\t\tp_value = %f\n\n", p_value < ALPHA ? "FAILURE" : "SUCCESS", p_value);
	fprintf(results[TEST_CUSUM], "%f\n", p_value);
#endif	
	// backwards
	sum1 = 0.0;
	for ( k=(-n/zrev+1)/4; k<=(n/zrev-1)/4; k++ ) {
		sum1 += cephes_normal(((4*k+1)*zrev)/sqrt(n));
		sum1 -= cephes_normal(((4*k-1)*zrev)/sqrt(n));
	}
	sum2 = 0.0;
	for ( k=(-n/zrev-3)/4; k<=(n/zrev-1)/4; k++ ) {
		sum2 += cephes_normal(((4*k+3)*zrev)/sqrt(n));
		sum2 -= cephes_normal(((4*k+1)*zrev)/sqrt(n));
	}
	p_value = 1.0 - sum1 + sum2;
#ifdef VERIFY_RESULTS
	R1.cusum.p_valueB=p_value;
	R1.cusum.sum1B=sum1;
	R1.cusum.sum2B=sum2;
#endif

#ifdef FILE_OUTPUT
	fprintf(stats[TEST_CUSUM], "\t\t      CUMULATIVE SUMS (REVERSE) TEST\n");
	fprintf(stats[TEST_CUSUM], "\t\t-------------------------------------------\n");
	fprintf(stats[TEST_CUSUM], "\t\tCOMPUTATIONAL INFORMATION:\n");
	fprintf(stats[TEST_CUSUM], "\t\t-------------------------------------------\n");
	fprintf(stats[TEST_CUSUM], "\t\t(a) The maximum partial sum = %d\n", zrev);
	fprintf(stats[TEST_CUSUM], "\t\t-------------------------------------------\n");

	if ( isNegative(p_value) || isGreaterThanOne(p_value) )
		fprintf(stats[TEST_CUSUM], "\t\tWARNING:  P_VALUE IS OUT OF RANGE\n");

	fprintf(stats[TEST_CUSUM], "%s\t\tp_value = %f\n\n", p_value < ALPHA ? "FAILURE" : "SUCCESS", p_value); fflush(stats[TEST_CUSUM]);
	fprintf(results[TEST_CUSUM], "%f\n", p_value); fflush(results[TEST_CUSUM]);
#endif
	//printf("sum1 %lf   sum2 %lf",sum1,sum2);
}

void
CumulativeSums2(int n)
{
	int		S, sup, inf, z, zrev, k, byte_ind;
	double	sum1, sum2, p_value;
	int max_minus[256] = {-8,-6,-6,-4,-6,-4,-4,-2,-6,-4,-4,-2,-4,-2,-2,0,-6,-4,-4,-2,-4,-2,-2,0,-4,-2,-2,0,-2,0,-1,0,-6,-4,-4,-2,-4,-2,-2,0,-4,-2,-2,0,-2,0,-1,0,-4,-2,-2,0,-2,0,-1,0,-3,-1,-1,0,-2,0,-1,0,-6,-4,-4,-2,-4,-2,-2,0,-4,-2,-2,0,-2,0,-1,0,-4,-2,-2,0,-2,0,-1,0,-3,-1,-1,0,-2,0,-1,0,-5,-3,-3,-1,-3,-1,-1,0,-3,-1,-1,0,-2,0,-1,0,-4,-2,-2,0,-2,0,-1,0,-3,-1,-1,0,-2,0,-1,0,-7,-5,-5,-3,-5,-3,-3,-1,-5,-3,-3,-1,-3,-1,-1,0,-5,-3,-3,-1,-3,-1,-1,0,-3,-1,-1,0,-2,0,-1,0,-5,-3,-3,-1,-3,-1,-1,0,-3,-1,-1,0,-2,0,-1,0,-4,-2,-2,0,-2,0,-1,0,-3,-1,-1,0,-2,0,-1,0,-6,-4,-4,-2,-4,-2,-2,0,-4,-2,-2,0,-2,0,-1,0,-4,-2,-2,0,-2,0,-1,0,-3,-1,-1,0,-2,0,-1,0,-5,-3,-3,-1,-3,-1,-1,0,-3,-1,-1,0,-2,0,-1,0,-4,-2,-2,0,-2,0,-1,0,-3,-1,-1,0,-2,0,-1,0};
	int max_plus[256] = {0,1,0,2,0,1,1,3,0,1,0,2,0,2,2,4,0,1,0,2,0,1,1,3,0,1,1,3,1,3,3,5,0,1,0,2,0,1,1,3,0,1,0,2,0,2,2,4,0,1,0,2,0,2,2,4,0,2,2,4,2,4,4,6,0,1,0,2,0,1,1,3,0,1,0,2,0,2,2,4,0,1,0,2,0,1,1,3,0,1,1,3,1,3,3,5,0,1,0,2,0,1,1,3,0,1,1,3,1,3,3,5,0,1,1,3,1,3,3,5,1,3,3,5,3,5,5,7,0,1,0,2,0,1,1,3,0,1,0,2,0,2,2,4,0,1,0,2,0,1,1,3,0,1,1,3,1,3,3,5,0,1,0,2,0,1,1,3,0,1,0,2,0,2,2,4,0,1,0,2,0,2,2,4,0,2,2,4,2,4,4,6,0,1,0,2,0,1,1,3,0,1,0,2,0,2,2,4,0,1,0,2,0,2,2,4,0,2,2,4,2,4,4,6,0,1,0,2,0,2,2,4,0,2,2,4,2,4,4,6,0,2,2,4,2,4,4,6,2,4,4,6,4,6,6,8};
	int byte_sum[256] = {-8,-6,-6,-4,-6,-4,-4,-2,-6,-4,-4,-2,-4,-2,-2,0,-6,-4,-4,-2,-4,-2,-2,0,-4,-2,-2,0,-2,0,0,2,-6,-4,-4,-2,-4,-2,-2,0,-4,-2,-2,0,-2,0,0,2,-4,-2,-2,0,-2,0,0,2,-2,0,0,2,0,2,2,4,-6,-4,-4,-2,-4,-2,-2,0,-4,-2,-2,0,-2,0,0,2,-4,-2,-2,0,-2,0,0,2,-2,0,0,2,0,2,2,4,-4,-2,-2,0,-2,0,0,2,-2,0,0,2,0,2,2,4,-2,0,0,2,0,2,2,4,0,2,2,4,2,4,4,6,-6,-4,-4,-2,-4,-2,-2,0,-4,-2,-2,0,-2,0,0,2,-4,-2,-2,0,-2,0,0,2,-2,0,0,2,0,2,2,4,-4,-2,-2,0,-2,0,0,2,-2,0,0,2,0,2,2,4,-2,0,0,2,0,2,2,4,0,2,2,4,2,4,4,6,-4,-2,-2,0,-2,0,0,2,-2,0,0,2,0,2,2,4,-2,0,0,2,0,2,2,4,0,2,2,4,2,4,4,6,-2,0,0,2,0,2,2,4,0,2,2,4,2,4,4,6,0,2,2,4,2,4,4,6,2,4,4,6,4,6,6,8};
unsigned int window;

	S = 0;
	sup = 0;
	inf = 0;
	for ( byte_ind = 0; byte_ind < n/8; byte_ind += 1 ) {

		if(S + max_plus[array[byte_ind]] > sup) sup = S + max_plus[array[byte_ind]];
		if(S + max_minus[array[byte_ind]] < inf) inf = S + max_minus[array[byte_ind]];
		S += byte_sum[array[byte_ind]];
	}
	//LAST BYTE
	if( n % 8 != 0)
	{
		for ( k=byte_ind*8; k<n; k++ ) {
			
			window = get_nth_block4(array,k);
			(window & 1) ? S++ : S--;

			if ( S > sup )
				sup++;
			if ( S < inf )
				inf--;
		}
	}


	z = (sup > -inf) ? sup : -inf;

	zrev = (sup-S > S-inf) ? sup-S : S-inf;

#ifdef VERIFY_RESULTS
	R2.cusum.z=z;
	R2.cusum.zrev=zrev;
#endif

	//printf("z = %d zrev=%d S= %d", z, zrev,S);
	
	// forward
	sum1 = 0.0;
	for ( k=(-n/z+1)/4; k<=(n/z-1)/4; k++ ) {
		sum1 += cephes_normal(((4*k+1)*z)/sqrt(n));
		sum1 -= cephes_normal(((4*k-1)*z)/sqrt(n));
	}
	sum2 = 0.0;
	for ( k=(-n/z-3)/4; k<=(n/z-1)/4; k++ ) {
		sum2 += cephes_normal(((4*k+3)*z)/sqrt(n));
		sum2 -= cephes_normal(((4*k+1)*z)/sqrt(n));
	}

	p_value = 1.0 - sum1 + sum2;
#ifdef VERIFY_RESULTS	
	R2.cusum.p_valueA=p_value;
	R2.cusum.sum1A=sum1;
	R2.cusum.sum2A=sum2;
#endif

#ifdef FILE_OUTPUT
	fprintf(stats[TEST_CUSUM], "\t\t      CUMULATIVE SUMS (FORWARD) TEST\n");
	fprintf(stats[TEST_CUSUM], "\t\t-------------------------------------------\n");
	fprintf(stats[TEST_CUSUM], "\t\tCOMPUTATIONAL INFORMATION:\n");
	fprintf(stats[TEST_CUSUM], "\t\t-------------------------------------------\n");
	fprintf(stats[TEST_CUSUM], "\t\t(a) The maximum partial sum = %d\n", z);
	fprintf(stats[TEST_CUSUM], "\t\t-------------------------------------------\n");

	if ( isNegative(p_value) || isGreaterThanOne(p_value) )
		fprintf(stats[TEST_CUSUM], "\t\tWARNING:  P_VALUE IS OUT OF RANGE\n");

	fprintf(stats[TEST_CUSUM], "%s\t\tp_value = %f\n\n", p_value < ALPHA ? "FAILURE" : "SUCCESS", p_value);
	fprintf(results[TEST_CUSUM], "%f\n", p_value);
#endif
	// backwards
	sum1 = 0.0;
	for ( k=(-n/zrev+1)/4; k<=(n/zrev-1)/4; k++ ) {
		sum1 += cephes_normal(((4*k+1)*zrev)/sqrt(n));
		sum1 -= cephes_normal(((4*k-1)*zrev)/sqrt(n));
	}
	sum2 = 0.0;
	for ( k=(-n/zrev-3)/4; k<=(n/zrev-1)/4; k++ ) {
		sum2 += cephes_normal(((4*k+3)*zrev)/sqrt(n));
		sum2 -= cephes_normal(((4*k+1)*zrev)/sqrt(n));
	}
	p_value = 1.0 - sum1 + sum2;
#ifdef VERIFY_RESULTS	
	R2.cusum.p_valueB=p_value;
	R2.cusum.sum1B=sum1;
	R2.cusum.sum2B=sum2;
#endif

#ifdef FILE_OUTPUT
	fprintf(stats[TEST_CUSUM], "\t\t      CUMULATIVE SUMS (REVERSE) TEST\n");
	fprintf(stats[TEST_CUSUM], "\t\t-------------------------------------------\n");
	fprintf(stats[TEST_CUSUM], "\t\tCOMPUTATIONAL INFORMATION:\n");
	fprintf(stats[TEST_CUSUM], "\t\t-------------------------------------------\n");
	fprintf(stats[TEST_CUSUM], "\t\t(a) The maximum partial sum = %d\n", zrev);
	fprintf(stats[TEST_CUSUM], "\t\t-------------------------------------------\n");

	if ( isNegative(p_value) || isGreaterThanOne(p_value) )
		fprintf(stats[TEST_CUSUM], "\t\tWARNING:  P_VALUE IS OUT OF RANGE\n");

	fprintf(stats[TEST_CUSUM], "%s\t\tp_value = %f\n\n", p_value < ALPHA ? "FAILURE" : "SUCCESS", p_value); fflush(stats[TEST_CUSUM]);
	fprintf(results[TEST_CUSUM], "%f\n", p_value); fflush(results[TEST_CUSUM]);
#endif
	//printf("sum1 %lf   sum2 %lf\n\n",sum1,sum2);
}