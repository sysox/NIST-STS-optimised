#include <stdio.h> 
#include <math.h> 
#include <string.h>
#include <stdlib.h>
#include "../include/externs.h"
#include "../include/cephes.h"
#include "../include/erf.h"
#include "../include/tools.h"


/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
            R A N D O M  E X C U R S I O N S  V A R I A N T  T E S T
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

void
RandomExcursionsVariant(int n)
{
	int		i, p, J, x, constraint, count, *S_k;
	int		stateX[18] = { -9, -8, -7, -6, -5, -4, -3, -2, -1, 1, 2, 3, 4, 5, 6, 7, 8, 9 };
	double	p_value;

#ifdef VERIFY_RESULTS
	R1.random_excursion_variant.valid=0;
#endif

	
	if ( (S_k = (int *)calloc(n, sizeof(int))) == NULL ) {
#ifdef FILE_OUTPUT
		fprintf(stats[TEST_RND_EXCURSION_VAR], "\t\tRANDOM EXCURSIONS VARIANT: Insufficent memory allocated.\n");
#endif
		return;
	}
	J = 0;
	S_k[0] = 2*(int)epsilon[0] - 1;
	for ( i=1; i<n; i++ ) {
		S_k[i] = S_k[i-1] + 2*epsilon[i] - 1;
		if ( S_k[i] == 0 )
			J++;
	}
	if ( S_k[n-1] != 0 )
		J++;

#ifdef FILE_OUTPUT
	fprintf(stats[TEST_RND_EXCURSION_VAR], "\t\t\tRANDOM EXCURSIONS VARIANT TEST\n");
	fprintf(stats[TEST_RND_EXCURSION_VAR], "\t\t--------------------------------------------\n");
	fprintf(stats[TEST_RND_EXCURSION_VAR], "\t\tCOMPUTATIONAL INFORMATION:\n");
	fprintf(stats[TEST_RND_EXCURSION_VAR], "\t\t--------------------------------------------\n");
	fprintf(stats[TEST_RND_EXCURSION_VAR], "\t\t(a) Number Of Cycles (J) = %d\n", J);
	fprintf(stats[TEST_RND_EXCURSION_VAR], "\t\t(b) Sequence Length (n)  = %d\n", n);
	fprintf(stats[TEST_RND_EXCURSION_VAR], "\t\t--------------------------------------------\n");
#endif

	constraint = (int)MAX(0.005*pow(n, 0.5), 500);
	if (J < constraint) {
#ifdef FILE_OUTPUT
		fprintf(stats[TEST_RND_EXCURSION_VAR], "\n\t\tWARNING:  TEST NOT APPLICABLE.  THERE ARE AN\n");
		fprintf(stats[TEST_RND_EXCURSION_VAR], "\t\t\t  INSUFFICIENT NUMBER OF CYCLES.\n");
		fprintf(stats[TEST_RND_EXCURSION_VAR], "\t\t---------------------------------------------\n");
		for ( i=0; i<18; i++ )
			fprintf(results[TEST_RND_EXCURSION_VAR], "%f\n", 0.0);
#endif
	}
	else {
		for ( p=0; p<=17; p++ ) {
			x = stateX[p];
			count = 0;
			for ( i=0; i<n; i++ )
				if ( S_k[i] == x )
					count++;
			//PRINT
			//printf("%d [%d]",count,J);
			p_value = erfc(fabs(count-J)/(sqrt(2.0*J*(4.0*fabs(x)-2))));

#ifdef FILE_OUTPUT
			if ( isNegative(p_value) || isGreaterThanOne(p_value) )
				fprintf(stats[TEST_RND_EXCURSION_VAR], "\t\t(b) WARNING: P_VALUE IS OUT OF RANGE.\n");
			fprintf(stats[TEST_RND_EXCURSION_VAR], "%s\t\t", p_value < ALPHA ? "FAILURE" : "SUCCESS");
			fprintf(stats[TEST_RND_EXCURSION_VAR], "(x = %2d) Total visits = %4d; p-value = %f\n", x, count, p_value);
			fprintf(results[TEST_RND_EXCURSION_VAR], "%f\n", p_value); fflush(results[TEST_RND_EXCURSION_VAR]);
#endif

#ifdef VERIFY_RESULTS
			R1.random_excursion_variant.valid=1;
			R1.random_excursion_variant.x[p]=x;
			R1.random_excursion_variant.count[p]=count;
			R1.random_excursion_variant.p_value[p]=p_value;
#endif

		}
	}
#ifdef FILE_OUTPUT
	fprintf(stats[TEST_RND_EXCURSION_VAR], "\n"); fflush(stats[TEST_RND_EXCURSION_VAR]);
#endif
	//printf("\n\n");
	free(S_k);
}

void
RandomExcursionsVariant2(int n)
{
	int		i, p, J = 0, x, constraint, count, S_k = 0, bit_ind,window;
	int		stateX[18] = { -9, -8, -7, -6, -5, -4, -3, -2, -1, 1, 2, 3, 4, 5, 6, 7, 8, 9 };
	double	p_value;
	int counter[19] = { 0, 0, 0, 0, 0, 0 , 0, 0, 0,  0, 0, 0, 0, 0, 0 , 0, 0, 0, 0 };

#ifdef VERIFY_RESULTS
	R2.random_excursion_variant.valid=0;
#endif

		for ( bit_ind = 0; bit_ind < n; bit_ind++ )
		{                          
			window = get_nth_block4(array,bit_ind);
			S_k += (window & 1)*2 - 1;
			if(S_k == 0)
			{
				++J;
			}
			if( (S_k >= -9) && (S_k <= 9)) counter[S_k+9]++;
		}
		////Last bit
		if(S_k)
		{
			J++;
		}
		for ( i=9; i<18; i++ )
		{
			counter[i] = counter[i+1];
		}
#ifdef FILE_OUTPUT
	fprintf(stats[TEST_RND_EXCURSION_VAR], "\t\t\tRANDOM EXCURSIONS VARIANT TEST\n");
	fprintf(stats[TEST_RND_EXCURSION_VAR], "\t\t--------------------------------------------\n");
	fprintf(stats[TEST_RND_EXCURSION_VAR], "\t\tCOMPUTATIONAL INFORMATION:\n");
	fprintf(stats[TEST_RND_EXCURSION_VAR], "\t\t--------------------------------------------\n");
	fprintf(stats[TEST_RND_EXCURSION_VAR], "\t\t(a) Number Of Cycles (J) = %d\n", J);
	fprintf(stats[TEST_RND_EXCURSION_VAR], "\t\t(b) Sequence Length (n)  = %d\n", n);
	fprintf(stats[TEST_RND_EXCURSION_VAR], "\t\t--------------------------------------------\n");
#endif

	constraint = (int)MAX(0.005*pow(n, 0.5), 500);
	if (J < constraint) {
#ifdef FILE_OUTPUT
		fprintf(stats[TEST_RND_EXCURSION_VAR], "\n\t\tWARNING:  TEST NOT APPLICABLE.  THERE ARE AN\n");
		fprintf(stats[TEST_RND_EXCURSION_VAR], "\t\t\t  INSUFFICIENT NUMBER OF CYCLES.\n");
		fprintf(stats[TEST_RND_EXCURSION_VAR], "\t\t---------------------------------------------\n");
		for ( i=0; i<18; i++ )
			fprintf(results[TEST_RND_EXCURSION_VAR], "%f\n", 0.0);
#endif
	}
	else {
		for ( p=0; p<=17; p++ ) {
			x = stateX[p];
			count = counter[p];
			//PRINT
			//printf("%d [%d]",count,J);
			p_value = erfc(fabs(count-J)/(sqrt(2.0*J*(4.0*fabs(x)-2))));

#ifdef VERIFY_RESULTS
			R2.random_excursion_variant.valid=1;
			R2.random_excursion_variant.x[p]=x;
			R2.random_excursion_variant.count[p]=count;
			R2.random_excursion_variant.p_value[p]=p_value;
#endif

#ifdef FILE_OUTPUT
			if ( isNegative(p_value) || isGreaterThanOne(p_value) )
				fprintf(stats[TEST_RND_EXCURSION_VAR], "\t\t(b) WARNING: P_VALUE IS OUT OF RANGE.\n");
			fprintf(stats[TEST_RND_EXCURSION_VAR], "%s\t\t", p_value < ALPHA ? "FAILURE" : "SUCCESS");
			fprintf(stats[TEST_RND_EXCURSION_VAR], "(x = %2d) Total visits = %4d; p-value = %f\n", x, count, p_value);
			fprintf(results[TEST_RND_EXCURSION_VAR], "%f\n", p_value); fflush(results[TEST_RND_EXCURSION_VAR]);
#endif
		}
	}
#ifdef FILE_OUTPUT
	fprintf(stats[TEST_RND_EXCURSION_VAR], "\n"); fflush(stats[TEST_RND_EXCURSION_VAR]);
#endif
}

