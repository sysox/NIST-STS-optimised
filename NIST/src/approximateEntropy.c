#include <stdio.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>
#include "../include/externs.h"
#include "../include/utilities.h"
#include "../include/cephes.h"  
#include "../include/tools.h"  

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
                A P P R O X I M A T E  E N T R O P Y   T E S T
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

void
ApproximateEntropy(int m, int n)
{
	int				i, j, k, r, blockSize, seqLength, powLen, index;
#ifdef VERIFY_RESULTS
	int cc=0;
#endif
	double			sum, numOfBlocks, ApEn[2], apen, chi_squared, p_value;
	unsigned int	*P;
	
#ifdef FILE_OUTPUT
	fprintf(stats[TEST_APEN], "\t\t\tAPPROXIMATE ENTROPY TEST\n");
	fprintf(stats[TEST_APEN], "\t\t--------------------------------------------\n");
	fprintf(stats[TEST_APEN], "\t\tCOMPUTATIONAL INFORMATION:\n");
	fprintf(stats[TEST_APEN], "\t\t--------------------------------------------\n");
	fprintf(stats[TEST_APEN], "\t\t(a) m (block length)    = %d\n", m);
#endif

	seqLength = n;
	r = 0;

#ifdef VERIFY_RESULTS
	R1.approximate_entropy.P=malloc(sizeof(unsigned int)*(1<<(m+1))*2);
	if(R1.approximate_entropy.P==NULL) {printf("Approximate entropy test: Cannot allocate memory.\n"); return; }
#endif

	for ( blockSize=m; blockSize<=m+1; blockSize++ ) {
		if ( blockSize == 0 ) {
			ApEn[0] = 0.00;
			r++;
		}
		else {
			numOfBlocks = (double)seqLength;
			powLen = (int)pow(2, blockSize+1)-1;
			if ( (P = (unsigned int*)calloc(powLen,sizeof(unsigned int)))== NULL ) {
#ifdef FILE_OUTPUT
				fprintf(stats[TEST_APEN], "ApEn:  Insufficient memory available.\n");
#endif
				printf("ApEn:  Insufficient memory available.\n");
				return;
			}
			for ( i=1; i<powLen-1; i++ )
				P[i] = 0;
			for ( i=0; i<numOfBlocks; i++ ) { /* COMPUTE FREQUENCY */
				/*epsilon[0] = 1;
				epsilon[1] = 1;
				epsilon[2] = 1;*/
				k = 1; 
				for ( j=0; j<blockSize; j++ ) {
					k <<= 1;
					if ( (int)epsilon[(i+j) % seqLength] == 1 )
						k++;
				}
				P[k-1]++;
				//if(i < 100)printf(" %i ",k-1);
			}
			/* DISPLAY FREQUENCY */
			sum = 0.0;
			index = (int)pow(2, blockSize)-1;
			for ( i=0; i<(int)pow(2, blockSize); i++ ) {
				if ( P[index] > 0 )
					sum += P[index]*log(P[index]/numOfBlocks);
				//printf("[%i: %d] ",index,P[index]);
#ifdef VERIFY_RESULTS
				R1.approximate_entropy.P[cc++]=P[index];
#endif
				index++;				
			}
#ifdef VERIFY_RESULTS
			R1.approximate_entropy.pp=cc;
#endif

			sum /= numOfBlocks;
			ApEn[r] = sum;

			//printf("\n");
			//printf("SUM: %lf \n\n",sum);
			r++;
			free(P);
		}
	}
	apen = ApEn[0] - ApEn[1];
	
	chi_squared = 2.0*seqLength*(log(2) - apen);
	p_value = cephes_igamc(pow(2, m-1), chi_squared/2.0);

#ifdef VERIFY_RESULTS
	R1.approximate_entropy.p_value=p_value;
	R1.approximate_entropy.chi_squared=chi_squared;
	R1.approximate_entropy.ApEn[0]=ApEn[0];
	R1.approximate_entropy.ApEn[1]=ApEn[1];
#endif


	//printf("P-value %lf \n",p_value);
#ifdef FILE_OUTPUT
	fprintf(stats[TEST_APEN], "\t\t(b) n (sequence length) = %d\n", seqLength);
	fprintf(stats[TEST_APEN], "\t\t(c) Chi^2               = %f\n", chi_squared);
	fprintf(stats[TEST_APEN], "\t\t(d) Phi(m)	       = %f\n", ApEn[0]);
	fprintf(stats[TEST_APEN], "\t\t(e) Phi(m+1)	       = %f\n", ApEn[1]);
	fprintf(stats[TEST_APEN], "\t\t(f) ApEn                = %f\n", apen);
	fprintf(stats[TEST_APEN], "\t\t(g) Log(2)              = %f\n", log(2.0));
	fprintf(stats[TEST_APEN], "\t\t--------------------------------------------\n");
#endif

	if ( m > (int)(log(seqLength)/log(2)-5) ) {
		//printf("\t\tNote: The blockSize exceeds recommended value\n");
#ifdef FILE_OUTPUT
		fprintf(stats[TEST_APEN], "\t\tNote: The blockSize = %d exceeds recommended value of %d\n", m,
			MAX(1, (int)(log(seqLength)/log(2)-5)));
		fprintf(stats[TEST_APEN], "\t\tResults are inaccurate!\n");
		fprintf(stats[TEST_APEN], "\t\t--------------------------------------------\n");
#endif
	}
	
#ifdef FILE_OUTPUT
	fprintf(stats[TEST_APEN], "%s\t\tp_value = %f\n\n", p_value < ALPHA ? "FAILURE" : "SUCCESS", p_value); fflush(stats[TEST_APEN]);
	fprintf(results[TEST_APEN], "%f\n", p_value); fflush(results[TEST_APEN]);
#endif
}


void
ApproximateEntropy2(int m, int n)
{
	int				i, k , len,cc=0;
	double			sum, numOfBlocks, ApEn[2], apen, chi_squared, p_value;
	unsigned int	*P,mask,help;

#ifdef FILE_OUTPUT
	fprintf(stats[TEST_APEN], "\t\t\tAPPROXIMATE ENTROPY TEST\n");
	fprintf(stats[TEST_APEN], "\t\t--------------------------------------------\n");
	fprintf(stats[TEST_APEN], "\t\tCOMPUTATIONAL INFORMATION:\n");
	fprintf(stats[TEST_APEN], "\t\t--------------------------------------------\n");
	fprintf(stats[TEST_APEN], "\t\t(a) m (block length)    = %d\n", m);
#endif

	
	numOfBlocks = n;
	m++;
	mask = (1 << m)-1;
	
	len = (1 << m);

#ifdef VERIFY_RESULTS
	R2.approximate_entropy.P=malloc(sizeof(unsigned int)*len*2);
	if(R2.approximate_entropy.P==NULL) {printf("Approximate entropy test: Cannot allocate memory.\n"); return; }
#endif
		
	if ( (P = (unsigned int*)calloc(len,sizeof(unsigned int)))== NULL ) {
#ifdef FILE_OUTPUT
		fprintf(stats[TEST_APEN], "ApEn:  Insufficient memory available.\n");
#endif
		return;
	}
	for ( i=0; i < len; i++ )
		P[i] = 0;
	for ( i=0; i < n-m+1; i++ ) {		 /* COMPUTE FREQUENCY */
		++P[get_nth_block4(array,i)&mask];
		//if(i < 100)printf(" %i ",(1 << m)- 1 + Mirrored_int((get_nth_block4(array,i)&mask),m));
	}
	for ( i=1; i<m; i++ ) {		
		k = get_nth_block4(array,n-m+i)&(mask>>i);
		//printf("%d ",k);
		k ^= (((unsigned int*)array)[0] << (m-i));
		//printf("%d ",k);
		k &= mask;
		//printf("%d ",k);
		P[k]++;
	}

	//DISPLAY FREQUENCY
	sum =  0.0;
	for ( i=0; i < len/2; i++ ) {
		help = P[Mirrored_int(i,m-1)]+P[Mirrored_int(i,m-1)+len/2];
		if ( help > 0 )
			sum += help*log(help/numOfBlocks);
#ifdef VERIFY_RESULTS
		R2.approximate_entropy.P[cc++]=help;
#endif
		//printf("%i ",help);	
	}
	
	
	sum /= numOfBlocks;
	ApEn[0] = sum;

	sum =  0.0;
	for ( i=0; i < len; i++ ) {
		if ( P[i] > 0 )
			sum += P[i]*log(P[i]/numOfBlocks);
		//printf("[%d: %d] ",i,P[Mirrored_int(i,m)]);	
#ifdef VERIFY_RESULTS
		R2.approximate_entropy.P[cc++]=P[Mirrored_int(i,m)];
#endif
	}
#ifdef VERIFY_RESULTS
	R2.approximate_entropy.pp=cc;
#endif
	
	sum /= numOfBlocks;
	ApEn[1] = sum;

	

	//printf("\n");
	//printf("SUM: %lf \n\n",sum);

	//printf("%lf %lf",ApEn[0],ApEn[1]);
	apen = ApEn[0] - ApEn[1];
	chi_squared = 2.0*n*(log(2) - apen);
	p_value = cephes_igamc(pow(2, m-2), chi_squared/2.0);

#ifdef VERIFY_RESULTS
	R2.approximate_entropy.p_value=p_value;
	R2.approximate_entropy.chi_squared=chi_squared;
	R2.approximate_entropy.ApEn[0]=ApEn[0];
	R2.approximate_entropy.ApEn[1]=ApEn[1];
#endif

	//printf("P-value %lf \n",p_value);
#ifdef FILE_OUTPUT
	fprintf(stats[TEST_APEN], "\t\t(b) n (sequence length) = %d\n", n);
	fprintf(stats[TEST_APEN], "\t\t(c) Chi^2               = %f\n", chi_squared);
	fprintf(stats[TEST_APEN], "\t\t(d) Phi(m)	       = %f\n", ApEn[0]);
	fprintf(stats[TEST_APEN], "\t\t(e) Phi(m+1)	       = %f\n", ApEn[1]);
	fprintf(stats[TEST_APEN], "\t\t(f) ApEn                = %f\n", apen);
	fprintf(stats[TEST_APEN], "\t\t(g) Log(2)              = %f\n", log(2.0));
	fprintf(stats[TEST_APEN], "\t\t--------------------------------------------\n");

	if ( m > (int)(log(n)/log(2)-5) ) {
		fprintf(stats[TEST_APEN], "\t\tNote: The blockSize = %d exceeds recommended value of %d\n", m,
			MAX(1, (int)(log(n)/log(2)-5)));
		fprintf(stats[TEST_APEN], "\t\tResults are inaccurate!\n");
		fprintf(stats[TEST_APEN], "\t\t--------------------------------------------\n");
	}
	
	fprintf(stats[TEST_APEN], "%s\t\tp_value = %f\n\n", p_value < ALPHA ? "FAILURE" : "SUCCESS", p_value); fflush(stats[TEST_APEN]);
	fprintf(results[TEST_APEN], "%f\n", p_value); fflush(results[TEST_APEN]);
#endif
	free(P);
}