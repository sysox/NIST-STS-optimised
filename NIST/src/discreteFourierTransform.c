#include <stdio.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>
#include "../include/externs.h"
#include "../include/utilities.h"
#include "../include/cephes.h"
#include "../include/erf.h"
#include "../include/tools.h"

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
         D I S C R E T E  F O U R I E R  T R A N S F O R M  T E S T 
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
void  __ogg_fdrffti(int n, double *wsave, int *ifac);
void  __ogg_fdrfftf(int n, double *X, double *wsave, int *ifac);

void
DiscreteFourierTransform(int n)
{
	double	p_value, upperBound, percentile, N_l, N_o, d, *m=NULL, *X=NULL, *wsave=NULL;
	int		i, count, ifac[15];

	if ( ((X = (double*) calloc(n,sizeof(double))) == NULL) ||
		 ((wsave = (double *)calloc(2*n,sizeof(double))) == NULL) ||
		 ((m = (double*)calloc(n/2+1, sizeof(double))) == NULL) ) {
#ifdef FILE_OUTPUT
			fprintf(stats[7],"\t\tUnable to allocate working arrays for the DFT.\n");
#endif
			 printf("\t\tUnable to allocate working arrays for the DFT.\n");
			if( X != NULL )
				free(X);
			if( wsave != NULL )
				free(wsave);
			if( m != NULL )
				free(m);
			return;
	}
	for ( i=0; i<n; i++ ){
		X[i] = 2*(int)epsilon[i] - 1;
	}
	
	__ogg_fdrffti(n, wsave, ifac);		/* INITIALIZE WORK ARRAYS */
	__ogg_fdrfftf(n, X, wsave, ifac);	/* APPLY FORWARD FFT */
	
	m[0] = sqrt(X[0]*X[0]);	    /* COMPUTE MAGNITUDE */
	
	for ( i=0; i<n/2; i++ )
		m[i+1] = sqrt(pow(X[2*i+1],2)+pow(X[2*i+2],2)); 
	count = 0;				       /* CONFIDENCE INTERVAL */
	upperBound = sqrt(2.995732274*n);
	for ( i=0; i<n/2; i++ )
		if ( m[i] < upperBound )
			count++;
	percentile = (double)count/(n/2)*100;
	N_l = (double) count;       /* number of peaks less than h = sqrt(3*n) */
	N_o = (double) 0.95*n/2.0;
	d = (N_l - N_o)/sqrt(n/4.0*0.95*0.05);
	p_value = erfc(fabs(d)/sqrt(2.0));
	//printf("%lf ",p_value);

#ifdef FILE_OUTPUT
	fprintf(stats[TEST_FFT], "\t\t\t\tFFT TEST\n");
	fprintf(stats[TEST_FFT], "\t\t-------------------------------------------\n");
	fprintf(stats[TEST_FFT], "\t\tCOMPUTATIONAL INFORMATION:\n");
	fprintf(stats[TEST_FFT], "\t\t-------------------------------------------\n");
	fprintf(stats[TEST_FFT], "\t\t(a) Percentile = %f\n", percentile);
	fprintf(stats[TEST_FFT], "\t\t(b) N_l        = %f\n", N_l);
  //fprintf(stats[TEST_FFT], "\t\t(c) N_o        = r%f\n", N_o); replaced by
	fprintf(stats[TEST_FFT], "\t\t(c) N_o        = %f\n", N_o); 
	fprintf(stats[TEST_FFT], "\t\t(d) d          = %f\n", d);
	fprintf(stats[TEST_FFT], "\t\t-------------------------------------------\n");

	fprintf(stats[TEST_FFT], "%s\t\tp_value = %f\n\n", p_value < ALPHA ? "FAILURE" : "SUCCESS", p_value);
	fprintf(results[TEST_FFT], "%f\n", p_value);
#endif

#ifdef VERIFY_RESULTS
	R1.dft.p_value=p_value;
	R1.dft.percentile=percentile;
	R1.dft.N_l=N_l;
	R1.dft.N_o=N_o;
	R1.dft.d=d;
#endif

	free(X);
	free(wsave);
	free(m);
}


void
DiscreteFourierTransform2(int n)
{
	double	p_value, upperBound, percentile, N_l, N_o, d, *m=NULL, *X=NULL, *wsave=NULL;
	int		i, count, ifac[15];

	if ( ((X = (double*) calloc(n,sizeof(double))) == NULL) ||
		 ((wsave = (double *)calloc(2*n,sizeof(double))) == NULL) ||
		 ((m = (double*)calloc(n/2+1, sizeof(double))) == NULL) ) {
#ifdef FILE_OUTPUT
			fprintf(stats[7],"\t\tUnable to allocate working arrays for the DFT.\n");
#endif
			printf("\t\tUnable to allocate working arrays for the DFT.\n");
			if( X != NULL )
				free(X);
			if( wsave != NULL )
				free(wsave);
			if( m != NULL )
				free(m);
			return;
	}
	for ( i=0; i<n; i++ )
	{
		X[i] = 2*((double)(get_nth_block4(array,i) & 1))-1;
	}
	
	__ogg_fdrffti(n, wsave, ifac);		/* INITIALIZE WORK ARRAYS */
	__ogg_fdrfftf(n, X, wsave, ifac);	/* APPLY FORWARD FFT */
	
	m[0] = sqrt(X[0]*X[0]);	    /* COMPUTE MAGNITUDE */
	
	for ( i=0; i<n/2; i++ )
		m[i+1] = sqrt(pow(X[2*i+1],2)+pow(X[2*i+2],2)); 
	count = 0;				       /* CONFIDENCE INTERVAL */
	upperBound = sqrt(2.995732274*n);
	for ( i=0; i<n/2; i++ )
		if ( m[i] < upperBound )
			count++;
	percentile = (double)count/(n/2)*100;
	N_l = (double) count;       /* number of peaks less than h = sqrt(3*n) */
	N_o = (double) 0.95*n/2.0;
	d = (N_l - N_o)/sqrt(n/4.0*0.95*0.05);
	p_value = erfc(fabs(d)/sqrt(2.0));
	//printf("%lf\n",p_value);

#ifdef FILE_OUTPUT
	fprintf(stats[TEST_FFT], "\t\t\t\tFFT TEST\n");
	fprintf(stats[TEST_FFT], "\t\t-------------------------------------------\n");
	fprintf(stats[TEST_FFT], "\t\tCOMPUTATIONAL INFORMATION:\n");
	fprintf(stats[TEST_FFT], "\t\t-------------------------------------------\n");
	fprintf(stats[TEST_FFT], "\t\t(a) Percentile = %f\n", percentile);
	fprintf(stats[TEST_FFT], "\t\t(b) N_l        = %f\n", N_l);
	fprintf(stats[TEST_FFT], "\t\t(c) N_o        = %f\n", N_o);
	fprintf(stats[TEST_FFT], "\t\t(d) d          = %f\n", d);
	fprintf(stats[TEST_FFT], "\t\t-------------------------------------------\n");

	fprintf(stats[TEST_FFT], "%s\t\tp_value = %f\n\n", p_value < ALPHA ? "FAILURE" : "SUCCESS", p_value);
	fprintf(results[TEST_FFT], "%f\n", p_value);
#endif

#ifdef VERIFY_RESULTS
	R2.dft.p_value=p_value;
	R2.dft.percentile=percentile;
	R2.dft.N_l=N_l;
	R2.dft.N_o=N_o;
	R2.dft.d=d;
#endif

	free(X);
	free(wsave);
	free(m);
}
