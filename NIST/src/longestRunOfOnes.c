/* got rid of unused 'k' */

#include <stdio.h>
#include <math.h>
#include <string.h>
#include "../include/externs.h"
#include "../include/cephes.h"  
#include "../include/tools.h"

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
                      L O N G E S T  R U N S  T E S T
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

void
LongestRunOfOnes(int n)
{
	double			pval, chi2, pi[7];
	int				run, v_n_obs, N, i, j, K, M, V[7];
	unsigned int	nu[7] = { 0, 0, 0, 0, 0, 0, 0 };

	if ( n < 128 ) {
#ifdef FILE_OUTPUT
		fprintf(stats[TEST_LONGEST_RUN], "\t\t\t  LONGEST RUNS OF ONES TEST\n");
		fprintf(stats[TEST_LONGEST_RUN], "\t\t---------------------------------------------\n");
		fprintf(stats[TEST_LONGEST_RUN], "\t\t   n=%d is too short\n", n);
#endif
		return;
	}
	if ( n < 6272 ) {
		K = 3;
		M = 8;
		V[0] = 1; V[1] = 2; V[2] = 3; V[3] = 4;
		pi[0] = 0.21484375;
		pi[1] = 0.3671875;
		pi[2] = 0.23046875;
		pi[3] = 0.1875;
	}
	else if ( n < 750000 ) {
		K = 5;
		M = 128;
		V[0] = 4; V[1] = 5; V[2] = 6; V[3] = 7; V[4] = 8; V[5] = 9;
		pi[0] = 0.1174035788;
		pi[1] = 0.242955959;
		pi[2] = 0.249363483;
		pi[3] = 0.17517706;
		pi[4] = 0.102701071;
		pi[5] = 0.112398847;
	}
	else {
		K = 6;
		M = 10000;
			V[0] = 10; V[1] = 11; V[2] = 12; V[3] = 13; V[4] = 14; V[5] = 15; V[6] = 16;
		pi[0] = 0.0882;
		pi[1] = 0.2092;
		pi[2] = 0.2483;
		pi[3] = 0.1933;
		pi[4] = 0.1208;
		pi[5] = 0.0675;
		pi[6] = 0.0727;
	}
	
	N = n/M;
	for ( i=0; i<N; i++ ) {
		v_n_obs = 0;
		run = 0;
		for ( j=0; j<M; j++ ) {
			if ( epsilon[i*M+j] == 1 ) {
				run++;
				if ( run > v_n_obs )
					v_n_obs = run;
			}
			else
				run = 0;
		}
		if ( v_n_obs < V[0] )
			nu[0]++;
		for ( j=0; j<=K; j++ ) {
			if ( v_n_obs == V[j] )
				nu[j]++;
		}
		if ( v_n_obs > V[K] )
			nu[K]++;
		//printf("%d ", v_n_obs);
	}

	chi2 = 0.0;
	for ( i=0; i<=K; i++ )
		chi2 += ((nu[i] - N * pi[i]) * (nu[i] - N * pi[i])) / (N * pi[i]);

	pval = cephes_igamc((double)(K/2.0), chi2 / 2.0);

	//for(i = 0; i < K; i++)
	//	printf("%d ",nu[i]);
#ifdef VERIFY_RESULTS
	R1.longestrunofones.N=N;
	R1.longestrunofones.M=M;
	R1.longestrunofones.chi2=chi2;
	for(i = 0; i < 7; i++) R1.longestrunofones.nu[i]=nu[i];
	R1.longestrunofones.pval=pval;
#endif

#ifdef FILE_OUTPUT
	fprintf(stats[TEST_LONGEST_RUN], "\t\t\t  LONGEST RUNS OF ONES TEST\n");
	fprintf(stats[TEST_LONGEST_RUN], "\t\t---------------------------------------------\n");
	fprintf(stats[TEST_LONGEST_RUN], "\t\tCOMPUTATIONAL INFORMATION:\n");
	fprintf(stats[TEST_LONGEST_RUN], "\t\t---------------------------------------------\n");
	fprintf(stats[TEST_LONGEST_RUN], "\t\t(a) N (# of substrings)  = %d\n", N);
	fprintf(stats[TEST_LONGEST_RUN], "\t\t(b) M (Substring Length) = %d\n", M);
	fprintf(stats[TEST_LONGEST_RUN], "\t\t(c) Chi^2                = %f\n", chi2);
	fprintf(stats[TEST_LONGEST_RUN], "\t\t---------------------------------------------\n");
	fprintf(stats[TEST_LONGEST_RUN], "\t\t      F R E Q U E N C Y\n");
	fprintf(stats[TEST_LONGEST_RUN], "\t\t---------------------------------------------\n");

	if ( K == 3 ) {
		fprintf(stats[TEST_LONGEST_RUN], "\t\t  <=1     2     3    >=4   P-value  Assignment");
		fprintf(stats[TEST_LONGEST_RUN], "\n\t\t %3d %3d %3d  %3d ", nu[0], nu[1], nu[2], nu[3]);
	}
	else if ( K == 5 ) {
		fprintf(stats[TEST_LONGEST_RUN], "\t\t<=4  5  6  7  8  >=9 P-value  Assignment");
		fprintf(stats[TEST_LONGEST_RUN], "\n\t\t %3d %3d %3d %3d %3d  %3d ", nu[0], nu[1], nu[2],
				nu[3], nu[4], nu[5]);
	}
	else {
		fprintf(stats[TEST_LONGEST_RUN],"\t\t<=10  11  12  13  14  15 >=16 P-value  Assignment");
		fprintf(stats[TEST_LONGEST_RUN],"\n\t\t %3d %3d %3d %3d %3d %3d  %3d ", nu[0], nu[1], nu[2],
				nu[3], nu[4], nu[5], nu[6]);
	}
	if ( isNegative(pval) || isGreaterThanOne(pval) )
		fprintf(stats[TEST_LONGEST_RUN], "WARNING:  P_VALUE IS OUT OF RANGE.\n");

	fprintf(stats[TEST_LONGEST_RUN], "%s\t\tp_value = %f\n\n", pval < ALPHA ? "FAILURE" : "SUCCESS", pval); fflush(stats[TEST_LONGEST_RUN]);
	fprintf(results[TEST_LONGEST_RUN], "%f\n", pval); fflush(results[TEST_LONGEST_RUN]);
#endif
}

void
LongestRunOfOnes2(int n)
{
	unsigned char  start_run[256]={0,1,0,2,0,1,0,3,0,1,0,2,0,1,0,4,0,1,0,2,0,1,0,3,0,1,0,2,0,1,0,5,0,1,0,2,0,1,0,3,0,1,0,2,0,1,0,4,0,1,0,2,0,1,0,3,0,1,0,2,0,1,0,6,0,1,0,2,0,1,0,3,0,1,0,2,0,1,0,4,0,1,0,2,0,1,0,3,0,1,0,2,0,1,0,5,0,1,0,2,0,1,0,3,0,1,0,2,0,1,0,4,0,1,0,2,0,1,0,3,0,1,0,2,0,1,0,7,0,1,0,2,0,1,0,3,0,1,0,2,0,1,0,4,0,1,0,2,0,1,0,3,0,1,0,2,0,1,0,5,0,1,0,2,0,1,0,3,0,1,0,2,0,1,0,4,0,1,0,2,0,1,0,3,0,1,0,2,0,1,0,6,0,1,0,2,0,1,0,3,0,1,0,2,0,1,0,4,0,1,0,2,0,1,0,3,0,1,0,2,0,1,0,5,0,1,0,2,0,1,0,3,0,1,0,2,0,1,0,4,0,1,0,2,0,1,0,3,0,1,0,2,0,1,0,8};
	unsigned char  end_run[256]  ={0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,4,4,4,4,4,4,4,4,5,5,5,5,6,6,7,8};
	unsigned char  max_run[256]  ={0,1,1,2,1,1,2,3,1,1,1,2,2,2,3,4,1,1,1,2,1,1,2,3,2,2,2,2,3,3,4,5,1,1,1,2,1,1,2,3,1,1,1,2,2,2,3,4,2,2,2,2,2,2,2,3,3,3,3,3,4,4,5,6,1,1,1,2,1,1,2,3,1,1,1,2,2,2,3,4,1,1,1,2,1,1,2,3,2,2,2,2,3,3,4,5,2,2,2,2,2,2,2,3,2,2,2,2,2,2,3,4,3,3,3,3,3,3,3,3,4,4,4,4,5,5,6,7,1,1,1,2,1,1,2,3,1,1,1,2,2,2,3,4,1,1,1,2,1,1,2,3,2,2,2,2,3,3,4,5,1,1,1,2,1,1,2,3,1,1,1,2,2,2,3,4,2,2,2,2,2,2,2,3,3,3,3,3,4,4,5,6,2,2,2,2,2,2,2,3,2,2,2,2,2,2,3,4,2,2,2,2,2,2,2,3,2,2,2,2,3,3,4,5,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,4,4,4,4,4,4,4,4,4,5,5,5,5,6,6,7,8};


	double			pval, chi2, pi[7];
	int				run, v_n_obs, N, i, j, K, M,Mbytes, V[7],rest_bits;
	unsigned int	nu[7] = { 0, 0, 0, 0, 0, 0, 0 };
	unsigned int end_mask[8] = {255,254,252,248,240,224,192,128};
	if ( n < 128 ) {
#ifdef FILE_OUTPUT
		fprintf(stats[TEST_LONGEST_RUN], "\t\t\t  LONGEST RUNS OF ONES TEST\n");
		fprintf(stats[TEST_LONGEST_RUN], "\t\t---------------------------------------------\n");
		fprintf(stats[TEST_LONGEST_RUN], "\t\t   n=%d is too short\n", n);
#endif
		return;
	}
	if ( n < 6272 ) {
		K = 3;
		M = 8;
		V[0] = 1; V[1] = 2; V[2] = 3; V[3] = 4;
		pi[0] = 0.21484375;
		pi[1] = 0.3671875;
		pi[2] = 0.23046875;
		pi[3] = 0.1875;
	}
	else if ( n < 750000 ) {
		K = 5;
		M = 128;
		V[0] = 4; V[1] = 5; V[2] = 6; V[3] = 7; V[4] = 8; V[5] = 9;
		pi[0] = 0.1174035788;
		pi[1] = 0.242955959;
		pi[2] = 0.249363483;
		pi[3] = 0.17517706;
		pi[4] = 0.102701071;
		pi[5] = 0.112398847;
	}
	else {
		K = 6;
		M = 10000;
			V[0] = 10; V[1] = 11; V[2] = 12; V[3] = 13; V[4] = 14; V[5] = 15; V[6] = 16;
		pi[0] = 0.0882;
		pi[1] = 0.2092;
		pi[2] = 0.2483;
		pi[3] = 0.1933;
		pi[4] = 0.1208;
		pi[5] = 0.0675;
		pi[6] = 0.0727;
	}
	
	N = n/M;
	Mbytes = M/8;
	//bits(array, 8);
	for ( i=0; i<N; i++ ) {

		//start byte of block
		rest_bits = i*M % 8;
		if(rest_bits != 0 )
		{
			j = i*M/8 + 1;
			v_n_obs = max_run[array[j] >> (8 - rest_bits)];
			run = end_run[((array[j] >> (8 - rest_bits)) << (8 - rest_bits))];
			//run = end_run[ array[j] & end_mask[rest_bits] ];
		}
		else
		{
			j = i*M/8;
			v_n_obs = max_run[array[j]];
			run = end_run[array[j]];
		}
		//j = i*M/8 + (rest_bits != 0);
		//v_n_obs = max_run[array[j] >> (8 - rest_bits)];
		//run = end_run[ array[j] & end_mask[rest_bits] ];
		
		

		j++;
		//whole bytes
		for (; j < M*(i+1)/8; j++ ) {
			run += start_run[array[j]];
			
			if ( run > v_n_obs ){
				v_n_obs = run;
			}
			if (  max_run[array[j]] > v_n_obs ){
				v_n_obs =  max_run[array[j]];
			}
			
			//v_n_obs = v_n_obs+( (v_n_obs-run) >> 31)*(v_n_obs-run);
		
			
			if(array[j] != 255) run = end_run[array[j]];
			//if(v_n_obs > V[K])break;		
		

		}
		//enb byte of block
		rest_bits = (i+1)*M % 8;
		if(rest_bits != 0 )
		{
			run += start_run[array[j] & (255 >> (8 - rest_bits))];;
			if ( run > v_n_obs ){
				v_n_obs = run;
			}
			
		}
		
		
		//printf("%d ", v_n_obs);
		if(v_n_obs < V[0])nu[0]++;
		else if(v_n_obs > V[K])nu[K]++;
		else nu[v_n_obs-V[0]]++;
		


	}


	chi2 = 0.0;
	for ( i=0; i<=K; i++ )
		chi2 += ((nu[i] - N * pi[i]) * (nu[i] - N * pi[i])) / (N * pi[i]);

	pval = cephes_igamc((double)(K/2.0), chi2 / 2.0);
	//for(i = 0; i < K; i++)
	//	printf("%d ",nu[i]);

#ifdef VERIFY_RESULTS
	R2.longestrunofones.N=N;
	R2.longestrunofones.M=M;
	R2.longestrunofones.chi2=chi2;
	for(i = 0; i < 7; i++) R2.longestrunofones.nu[i]=nu[i];
	R2.longestrunofones.pval=pval;
#endif

#ifdef FILE_OUTPUT
	fprintf(stats[TEST_LONGEST_RUN], "\t\t\t  LONGEST RUNS OF ONES TEST\n");
	fprintf(stats[TEST_LONGEST_RUN], "\t\t---------------------------------------------\n");
	fprintf(stats[TEST_LONGEST_RUN], "\t\tCOMPUTATIONAL INFORMATION:\n");
	fprintf(stats[TEST_LONGEST_RUN], "\t\t---------------------------------------------\n");
	fprintf(stats[TEST_LONGEST_RUN], "\t\t(a) N (# of substrings)  = %d\n", N);
	fprintf(stats[TEST_LONGEST_RUN], "\t\t(b) M (Substring Length) = %d\n", M);
	fprintf(stats[TEST_LONGEST_RUN], "\t\t(c) Chi^2                = %f\n", chi2);
	fprintf(stats[TEST_LONGEST_RUN], "\t\t---------------------------------------------\n");
	fprintf(stats[TEST_LONGEST_RUN], "\t\t      F R E Q U E N C Y\n");
	fprintf(stats[TEST_LONGEST_RUN], "\t\t---------------------------------------------\n");

	if ( K == 3 ) {
		fprintf(stats[TEST_LONGEST_RUN], "\t\t  <=1     2     3    >=4   P-value  Assignment");
		fprintf(stats[TEST_LONGEST_RUN], "\n\t\t %3d %3d %3d  %3d ", nu[0], nu[1], nu[2], nu[3]);
	}
	else if ( K == 5 ) {
		fprintf(stats[TEST_LONGEST_RUN], "\t\t<=4  5  6  7  8  >=9 P-value  Assignment");
		fprintf(stats[TEST_LONGEST_RUN], "\n\t\t %3d %3d %3d %3d %3d  %3d ", nu[0], nu[1], nu[2],
				nu[3], nu[4], nu[5]);
	}
	else {
		fprintf(stats[TEST_LONGEST_RUN],"\t\t<=10  11  12  13  14  15 >=16 P-value  Assignment");
		fprintf(stats[TEST_LONGEST_RUN],"\n\t\t %3d %3d %3d %3d %3d %3d  %3d ", nu[0], nu[1], nu[2],
				nu[3], nu[4], nu[5], nu[6]);
	}
	if ( isNegative(pval) || isGreaterThanOne(pval) )
		fprintf(stats[TEST_LONGEST_RUN], "WARNING:  P_VALUE IS OUT OF RANGE.\n");

	fprintf(stats[TEST_LONGEST_RUN], "%s\t\tp_value = %f\n\n", pval < ALPHA ? "FAILURE" : "SUCCESS", pval); fflush(stats[TEST_LONGEST_RUN]);
	fprintf(results[TEST_LONGEST_RUN], "%f\n", pval); fflush(results[TEST_LONGEST_RUN]);
#endif
}
