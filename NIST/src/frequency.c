#include <stdio.h>
#include <math.h>
#include <string.h>
#include "../include/externs.h"
#include "../include/erf.h"
#include "../include/tools.h"

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
                          F R E Q U E N C Y  T E S T
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

void
Frequency(int n)
{
	int		i;
	double	f, s_obs, p_value, sum, sqrt2 = 1.41421356237309504880;
	
	sum = 0.0;
	for ( i=0; i<n; i++ )
		sum += 2*(int)epsilon[i]-1;
		//sum += epsilon[i];
	s_obs = fabs(sum)/sqrt(n);
	f = s_obs/sqrt2;
	p_value = erfc(f);
	//printf("Pval: %lf\n",p_value);
	//printf("%lf ",sum);
#ifdef VERIFY_RESULTS
	R1.frequency.sum=sum;
	R1.frequency.sum_n=sum/n;
	R1.frequency.p_value=p_value;
#endif
	
#ifdef FILE_OUTPUT
	fprintf(stats[TEST_FREQUENCY], "\t\t\t      FREQUENCY TEST\n");
	fprintf(stats[TEST_FREQUENCY], "\t\t---------------------------------------------\n");
	fprintf(stats[TEST_FREQUENCY], "\t\tCOMPUTATIONAL INFORMATION:\n");
	fprintf(stats[TEST_FREQUENCY], "\t\t---------------------------------------------\n");
	fprintf(stats[TEST_FREQUENCY], "\t\t(a) The nth partial sum = %d\n", (int)sum);
	fprintf(stats[TEST_FREQUENCY], "\t\t(b) S_n/n               = %f\n", sum/n);
	fprintf(stats[TEST_FREQUENCY], "\t\t---------------------------------------------\n");

	fprintf(stats[TEST_FREQUENCY], "%s\t\tp_value = %f\n\n", p_value < ALPHA ? "FAILURE" : "SUCCESS", p_value); fflush(stats[TEST_FREQUENCY]);
	fprintf(results[TEST_FREQUENCY], "%f\n", p_value); fflush(results[TEST_FREQUENCY]);
#endif
}

/* --------------------------------------------------------------------------

The following code is distributed under the following BSD-style license:

Copyright � 2013-2014 Marek Sys (syso@fi.muni.cz) & Zdenek Riha (zriha@fi.muni.cz).
All Rights Reserved.

Redistribution and use in source and binary forms, with or without modification,
are permitted provided that the following conditions are met:

1. Redistributions of source code must retain the above copyright notice, this
list of conditions and the following disclaimer.

2. Redistributions in binary form must reproduce the above copyright notice, this
list of conditions and the following disclaimer in the documentation and/or other
materials provided with the distribution.

3. The name of the author may not be used to endorse or promote products derived
from this software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY AUTHORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES,
INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS
FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE AUTHORS BE LIABLE FOR
ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
(INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF
USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR
OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF
THE POSSIBILITY OF SUCH DAMAGE.

-------------------------------------------------------------------------- */

void
Frequency2(int n)
{
	int		i,int_sum;
	double	f, s_obs, p_value, sum, sqrt2 = 1.41421356237309504880;
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
	
	sum = 0.0;
	int_sum = 0;

	for ( i=0; i*8< n; i++ )
	{
		int_sum += LU_byte_weight[array[i]];
		//sum += LU_byte_weight[(get_nth_block4(array,i*8) & 255)];
	}
	if( i*8 > n)int_sum -= LU_byte_weight[array[i-1] >> (8-(i*8-n))];

	//sum = count_of_ones - (n-count_of_ones);
	sum = int_sum-(n-int_sum);
	s_obs = fabs(sum)/sqrt(n);
	f = s_obs/sqrt2;
	p_value = erfc(f);
	//printf("Pval: %lf\n",p_value);
	//printf("%lf ",(double)sum);
#ifdef VERIFY_RESULTS
	R2.frequency.sum=sum;
	R2.frequency.sum_n=sum/n;
	R2.frequency.p_value=p_value;
#endif

#ifdef FILE_OUTPUT
	fprintf(stats[TEST_FREQUENCY], "\t\t\t      FREQUENCY TEST\n");
	fprintf(stats[TEST_FREQUENCY], "\t\t---------------------------------------------\n");
	fprintf(stats[TEST_FREQUENCY], "\t\tCOMPUTATIONAL INFORMATION:\n");
	fprintf(stats[TEST_FREQUENCY], "\t\t---------------------------------------------\n");
	fprintf(stats[TEST_FREQUENCY], "\t\t(a) The nth partial sum = %d\n", (int)sum);
	fprintf(stats[TEST_FREQUENCY], "\t\t(b) S_n/n               = %f\n", sum/n);
	fprintf(stats[TEST_FREQUENCY], "\t\t---------------------------------------------\n");

	fprintf(stats[TEST_FREQUENCY], "%s\t\tp_value = %f\n\n", p_value < ALPHA ? "FAILURE" : "SUCCESS", p_value); fflush(stats[TEST_FREQUENCY]);
	fprintf(results[TEST_FREQUENCY], "%f\n", p_value); fflush(results[TEST_FREQUENCY]);
#endif
}

void
Frequency3(int n)
{
	int int_sum, mask;
	double	f, s_obs, p_value, sum, sqrt2 = 1.41421356237309504880;
	unsigned char *p_tmp, *p_end;

	int LUT_HW_size = 16;
	int LUT_HW_Bsize = 2;
	signed char *LUT_HW = LUT_HW_16;

	if(1)
	{
		LUT_HW_size = 8;
		LUT_HW_Bsize = 1;
		LUT_HW = LUT_HW_8;
	}

	sum = 0.0;
	int_sum = 0;
	mask = get_mask(LUT_HW_size);

	/*for (processed_bits = 0; processed_bits + LUT_HW_size < n + 1; processed_bits += LUT_HW_size){
		int_sum += LUT_HW[get_block_fast(array, Boffset) & mask];
		Boffset += LUT_HW_Bsize;
	}*/

	
	p_end = array + (n- (n%LUT_HW_size))/8;
	for (p_tmp = array; p_tmp < p_end; p_tmp += LUT_HW_Bsize){
		int_sum += LUT_HW[*((unsigned int*)p_tmp) & mask];
	}
	if (n % LUT_HW_size){
		int_sum += LUT_HW[ *((unsigned int*)p_tmp) & get_mask(n % LUT_HW_size)];
	}

	sum = int_sum - (n - int_sum);
	s_obs = fabs(sum) / sqrt(n);
	f = s_obs / sqrt2;
	p_value = erfc(f);

#ifdef VERIFY_RESULTS
	R2.frequency.sum = sum;
	R2.frequency.sum_n = sum / n;
	R2.frequency.p_value = p_value;
#endif

#ifdef FILE_OUTPUT
	fprintf(stats[TEST_FREQUENCY], "\t\t\t      FREQUENCY TEST\n");
	fprintf(stats[TEST_FREQUENCY], "\t\t---------------------------------------------\n");
	fprintf(stats[TEST_FREQUENCY], "\t\tCOMPUTATIONAL INFORMATION:\n");
	fprintf(stats[TEST_FREQUENCY], "\t\t---------------------------------------------\n");
	fprintf(stats[TEST_FREQUENCY], "\t\t(a) The nth partial sum = %d\n", (int)sum);
	fprintf(stats[TEST_FREQUENCY], "\t\t(b) S_n/n               = %f\n", sum / n);
	fprintf(stats[TEST_FREQUENCY], "\t\t---------------------------------------------\n");

	fprintf(stats[TEST_FREQUENCY], "%s\t\tp_value = %f\n\n", p_value < ALPHA ? "FAILURE" : "SUCCESS", p_value); fflush(stats[TEST_FREQUENCY]);
	fprintf(results[TEST_FREQUENCY], "%f\n", p_value); fflush(results[TEST_FREQUENCY]);
#endif
}
