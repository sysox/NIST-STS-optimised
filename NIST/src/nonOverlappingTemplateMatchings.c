#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "../include/externs.h"
#include "../include/utilities.h"
#include "../include/cephes.h"  
#include "../include/tools.h"  

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
          N O N O V E R L A P P I N G  T E M P L A T E  T E S T
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

void
NonOverlappingTemplateMatchings(int m, int n)
{
	int		numOfTemplates[100] = {0, 0, 2, 4, 6, 12, 20, 40, 74, 148, 284, 568, 1116,
						2232, 4424, 8848, 17622, 35244, 70340, 140680, 281076, 562152};
	/*----------------------------------------------------------------------------
	NOTE:  Should additional templates lengths beyond 21 be desired, they must 
	first be constructed, saved into files and then the corresponding 
	number of nonperiodic templates for that file be stored in the m-th 
	position in the numOfTemplates variable.
	----------------------------------------------------------------------------*/
	unsigned int	bit, W_obs, nu[6], *Wj = NULL; 
	FILE			*fp=NULL;
	double			sum, chi2, p_value, lambda, pi[6], varWj;
	int				i, j, jj, k, match, SKIP, M, N, K = 5;
	char			directory[100];
	BitSequence		*sequence = NULL;

	N = 8;
	M = n/N;

	if ( (Wj = (unsigned int*)calloc(N, sizeof(unsigned int))) == NULL ) {
#ifdef FILE_OUTPUT
		fprintf(stats[TEST_NONPERIODIC], "\tNONOVERLAPPING TEMPLATES TESTS ABORTED DUE TO ONE OF THE FOLLOWING : \n");
		fprintf(stats[TEST_NONPERIODIC], "\tInsufficient memory for required work space.\n");
#endif
		printf("\tNONOVERLAPPING TEMPLATES TESTS: Insufficient memory for required work space.\n");
		return;
	}
	lambda = (M-m+1)/pow(2, m);
	varWj = M*(1.0/pow(2.0, m) - (2.0*m-1.0)/pow(2.0, 2.0*m));
	sprintf(directory, "templates/template%d", m);

	if ( ((isNegative(lambda)) || (isZero(lambda))) ||
		 ((fp = fopen(directory, "r")) == NULL) ||
		 ((sequence = (BitSequence *) calloc(m, sizeof(BitSequence))) == NULL) ) {
#ifdef FILE_OUTPUT
		fprintf(stats[TEST_NONPERIODIC], "\tNONOVERLAPPING TEMPLATES TESTS ABORTED DUE TO ONE OF THE FOLLOWING : \n");
		fprintf(stats[TEST_NONPERIODIC], "\tLambda (%f) not being positive!\n", lambda);
		fprintf(stats[TEST_NONPERIODIC], "\tTemplate file <%s> not existing\n", directory);
		fprintf(stats[TEST_NONPERIODIC], "\tInsufficient memory for required work space.\n");
#endif
		printf("\tNONOVERLAPPING TEMPLATES TESTS ABORTED.\n");
		if ( sequence != NULL )
			free(sequence);

#ifdef VERIFY_RESULTS
		R1.nonoverlapping.templates=0;
		R1.nonoverlapping.W=NULL;
		R1.nonoverlapping.chi2=NULL;
		R1.nonoverlapping.p_value=NULL;
#endif
	}
	else {
#ifdef FILE_OUTPUT
		fprintf(stats[TEST_NONPERIODIC], "\t\t  NONPERIODIC TEMPLATES TEST\n");
		fprintf(stats[TEST_NONPERIODIC], "-------------------------------------------------------------------------------------\n");
		fprintf(stats[TEST_NONPERIODIC], "\t\t  COMPUTATIONAL INFORMATION\n");
		fprintf(stats[TEST_NONPERIODIC], "-------------------------------------------------------------------------------------\n");
		fprintf(stats[TEST_NONPERIODIC], "\tLAMBDA = %f\tM = %d\tN = %d\tm = %d\tn = %d\n", lambda, M, N, m, n);
		fprintf(stats[TEST_NONPERIODIC], "-------------------------------------------------------------------------------------\n");
		fprintf(stats[TEST_NONPERIODIC], "\t\tF R E Q U E N C Y\n");
		fprintf(stats[TEST_NONPERIODIC], "Template   W_1  W_2  W_3  W_4  W_5  W_6  W_7  W_8    Chi^2   P_value Assignment Index\n");
		fprintf(stats[TEST_NONPERIODIC], "-------------------------------------------------------------------------------------\n");
#endif

		if ( numOfTemplates[m] < MAXNUMOFTEMPLATES )
			SKIP = 1;
		else
			SKIP = (int)(numOfTemplates[m]/MAXNUMOFTEMPLATES);
		numOfTemplates[m] = (int)numOfTemplates[m]/SKIP;

#ifdef VERIFY_RESULTS
		R1.nonoverlapping.templates=MIN(MAXNUMOFTEMPLATES, numOfTemplates[m]);
		R1.nonoverlapping.W=(unsigned int*)malloc(sizeof(unsigned int)*R1.nonoverlapping.templates*8);
		R1.nonoverlapping.chi2=(double *)malloc(sizeof(double)*R1.nonoverlapping.templates);
		R1.nonoverlapping.p_value=(double *)malloc(sizeof(double)*R1.nonoverlapping.templates);
		if(R1.nonoverlapping.W==NULL||R1.nonoverlapping.chi2==NULL||R1.nonoverlapping.p_value==NULL)
		{ printf("NONOVERLAPPING TEMPLATES TEST: Cannot allocate memory"); return; }
#endif
		
		
		sum = 0.0;
		for ( i=0; i<2; i++ ) {                      /* Compute Probabilities */
			pi[i] = exp(-lambda+i*log(lambda)-cephes_lgam(i+1));
			sum += pi[i];
		}
		pi[0] = sum;
		for ( i=2; i<=K; i++ ) {                      /* Compute Probabilities */
			pi[i-1] = exp(-lambda+i*log(lambda)-cephes_lgam(i+1));
			sum += pi[i-1];
		}
		pi[K] = 1 - sum;

		for( jj=0; jj<MIN(MAXNUMOFTEMPLATES, numOfTemplates[m]); jj++ ) {
			sum = 0;

			for ( k=0; k<m; k++ ) {
				fscanf(fp, "%d", &bit);
				sequence[k] = bit;
				//printf("%d", sequence[k]);
#ifdef FILE_OUTPUT
				fprintf(stats[TEST_NONPERIODIC], "%d", sequence[k]);
#endif
			}
			//printf(" ");
#ifdef FILE_OUTPUT
			fprintf(stats[TEST_NONPERIODIC], " ");
#endif
			for ( k=0; k<=K; k++ )
				nu[k] = 0;
			for ( i=0; i<N; i++ ) {
				W_obs = 0;
				for ( j=0; j<M-m+1; j++ ) {
					match = 1;
					for ( k=0; k<m; k++ ) {
						if ( (int)sequence[k] != (int)epsilon[i*M+j+k] ) {
							match = 0;
							break;
						}
					}
					if ( match == 1 )
						W_obs++;
				}
				Wj[i] = W_obs;
			}
			sum = 0;
			chi2 = 0.0;                                   /* Compute Chi Square */
			for ( i=0; i<N; i++ ) {
#ifdef FILE_OUTPUT
				if ( m == 10 )
					fprintf(stats[TEST_NONPERIODIC], "%3d  ", Wj[i]);
				else
					fprintf(stats[TEST_NONPERIODIC], "%4d ", Wj[i]);
#endif
				/*if ( m == 10 )
					printf("%3d  ", Wj[i]);
				else
					printf("%4d ", Wj[i]);*/
#ifdef VERIFY_RESULTS
				R1.nonoverlapping.W[jj*N+i]=Wj[i];
#endif
				//
				chi2 += pow(((double)Wj[i] - lambda)/pow(varWj, 0.5), 2);
			}
			p_value = cephes_igamc(N/2.0, chi2/2.0);
			//printf("(1) chi2:%lf value: %lf\n",chi2,p_value);
#ifdef VERIFY_RESULTS
			R1.nonoverlapping.chi2[jj]=chi2;
			R1.nonoverlapping.p_value[jj]=p_value;
#endif
			
#ifdef FILE_OUTPUT
			if ( isNegative(p_value) || isGreaterThanOne(p_value) )
				fprintf(stats[TEST_NONPERIODIC], "\t\tWARNING:  P_VALUE IS OUT OF RANGE.\n");

			fprintf(stats[TEST_NONPERIODIC], "%9.6f %f %s %3d\n", chi2, p_value, p_value < ALPHA ? "FAILURE" : "SUCCESS", jj);
#endif
			if ( SKIP > 1 )
				fseek(fp, (long)(SKIP-1)*2*m, SEEK_CUR);
#ifdef FILE_OUTPUT
			fprintf(results[TEST_NONPERIODIC], "%f\n", p_value); fflush(results[TEST_NONPERIODIC]);
#endif
		}
	}
	
#ifdef FILE_OUTPUT
	fprintf(stats[TEST_NONPERIODIC], "\n"); fflush(stats[TEST_NONPERIODIC]);
#endif
	if ( sequence != NULL )
		free(sequence);

	free(Wj);
	if(fp)fclose(fp);
}

//best
void
NonOverlappingTemplateMatchings2(int m, int n)
{
	int		numOfTemplates[100] = {0, 0, 2, 4, 6, 12, 20, 40, 74, 148, 284, 568, 1116,
						2232, 4424, 8848, 17622, 35244, 70340, 140680, 281076, 562152};
	/*----------------------------------------------------------------------------
	NOTE:  Should additional templates lengths beyond 21 be desired, they must 
	first be constructed, saved into files and then the corresponding 
	number of nonperiodic templates for that file be stored in the m-th 
	position in the numOfTemplates variable.
	----------------------------------------------------------------------------*/ 
	FILE			*fp=NULL;
	double			sum, chi2, p_value, lambda, pi[6], varWj;
	int	i, j, k,  M, N, K = 5,SKIP;
	char			directory[100];

	unsigned int sequence,numoftemplates,window,one_template, **Wj = NULL, *templates,mask;
	int bit,bit_ind;
	int block;

	N = 8;
	M = n/N;

	if ( numOfTemplates[m] < MAXNUMOFTEMPLATES )
			SKIP = 1;
	else	
		SKIP = (int)(numOfTemplates[m]/MAXNUMOFTEMPLATES);
	
	
	numoftemplates = numOfTemplates[m];

	
	templates = (unsigned int*)malloc(numoftemplates*sizeof(unsigned int));
	if(templates==NULL) 
	{ 
#ifdef FILE_OUTPUT
		fprintf(stats[TEST_NONPERIODIC], "\tNONOVERLAPPING TEMPLATES TESTS: CANNOT ALLOCATE MEMORY\n");
#endif
		printf("NONOVERLAPPING TEMPLATES TESTS: Cannot allocate memory.\n"); 
		return; 
	}
	//Wj = (unsigned int**)malloc((1 << m)*sizeof(unsigned int*));
	Wj = (unsigned int**)malloc(N*sizeof(unsigned int*));
	if(Wj==NULL) 
	{ 
#ifdef FILE_OUTPUT
		fprintf(stats[TEST_NONPERIODIC], "\tNONOVERLAPPING TEMPLATES TESTS: CANNOT ALLOCATE MEMORY\n");
#endif
		printf("NONOVERLAPPING TEMPLATES TESTS: Cannot allocate memory.\n"); 
		return; 
	}	
	for(i = 0; i < N; i++)
	{
		Wj[i] = (unsigned int*)malloc((1 << m)*sizeof(unsigned int));
		if(Wj[i]==NULL) 
		{
#ifdef FILE_OUTPUT
			fprintf(stats[TEST_NONPERIODIC], "\tNONOVERLAPPING TEMPLATES TESTS: CANNOT ALLOCATE MEMORY\n");
#endif
			printf("NONOVERLAPPING TEMPLATES TESTS: Cannot allocate memory.\n"); 

			for(j = 0; j < i; j++)
			{
				free(Wj[j]);
			}
			free(Wj);
			return; 
		}
		for(j = 0; j < (1 << m); j++)
		{
			Wj[i][j] = 0;
		}
	}
	
	mask = (1 << m ) - 1;
	lambda = (M-m+1)/pow(2, m);
	varWj = M*(1.0/pow(2.0, m) - (2.0*m-1.0)/pow(2.0, 2.0*m));
	
	sprintf(directory, "templates/template%d", m);
	//fp=fopen(directory, "r");
	
	if ( ((isNegative(lambda)) || (isZero(lambda))) ||
		 ((fp = fopen(directory, "r")) == NULL) 
	    ) 
	{
#ifdef FILE_OUTPUT
		fprintf(stats[TEST_NONPERIODIC], "\tNONOVERLAPPING TEMPLATES TESTS ABORTED DUE TO ONE OF THE FOLLOWING : \n");
		fprintf(stats[TEST_NONPERIODIC], "\tLambda (%f) not being positive!\n", lambda);
		fprintf(stats[TEST_NONPERIODIC], "\tTemplate file <%s> not existing\n", directory);
		fprintf(stats[TEST_NONPERIODIC], "\tInsufficient memory for required work space.\n");
#endif
		//if(fp==NULL)
		//	printf("NONOVERLAPPING TEMPLATES TESTS: Cannot open templates file.\n"); 

#ifdef VERIFY_RESULTS
		R2.nonoverlapping.templates=0;
		R2.nonoverlapping.W=(unsigned int*)malloc(sizeof(unsigned int)*R2.nonoverlapping.templates*8);
		R2.nonoverlapping.chi2=(double *)malloc(sizeof(double)*R2.nonoverlapping.templates);
		R2.nonoverlapping.p_value=(double *)malloc(sizeof(double)*R2.nonoverlapping.templates);
		if(R2.nonoverlapping.W==NULL||R2.nonoverlapping.chi2==NULL||R2.nonoverlapping.p_value==NULL)
		{ printf("NONOVERLAPPING TEMPLATES TEST: Cannot allocate memory"); return; }
#endif
	}
	else
	{	
#ifdef VERIFY_RESULTS
		R2.nonoverlapping.templates=MIN(numoftemplates,MAXNUMOFTEMPLATES);
		R2.nonoverlapping.W=(unsigned int*)malloc(sizeof(unsigned int)*R2.nonoverlapping.templates*8);
		if(R2.nonoverlapping.W==NULL)
		{
			printf("NONOVERLAPPING TEMPLATES TESTS: Cannot allocate memory.\n"); 
			return; 
		}
		R2.nonoverlapping.chi2=(double *)malloc(sizeof(double)*R2.nonoverlapping.templates);
		if(R2.nonoverlapping.chi2==NULL)
		{ 
			printf("NONOVERLAPPING TEMPLATES TESTS: Cannot allocate memory.\n"); 
			return; 
		}
		R2.nonoverlapping.p_value=(double *)malloc(sizeof(double)*R2.nonoverlapping.templates);
		if(R2.nonoverlapping.p_value==NULL) 
		{ 
			printf("NONOVERLAPPING TEMPLATES TESTS: Cannot allocate memory.\n"); 
			return; 
		}
#endif

#ifdef FILE_OUTPUT
		fprintf(stats[TEST_NONPERIODIC], "\t\t  NONPERIODIC TEMPLATES TEST\n");
		fprintf(stats[TEST_NONPERIODIC], "-------------------------------------------------------------------------------------\n");
		fprintf(stats[TEST_NONPERIODIC], "\t\t  COMPUTATIONAL INFORMATION\n");
		fprintf(stats[TEST_NONPERIODIC], "-------------------------------------------------------------------------------------\n");
		fprintf(stats[TEST_NONPERIODIC], "\tLAMBDA = %f\tM = %d\tN = %d\tm = %d\tn = %d\n", lambda, M, N, m, n);
		fprintf(stats[TEST_NONPERIODIC], "-------------------------------------------------------------------------------------\n");
		fprintf(stats[TEST_NONPERIODIC], "\t\tF R E Q U E N C Y\n");
		fprintf(stats[TEST_NONPERIODIC], "Template   W_1  W_2  W_3  W_4  W_5  W_6  W_7  W_8    Chi^2   P_value Assignment Index\n");
		fprintf(stats[TEST_NONPERIODIC], "-------------------------------------------------------------------------------------\n");
#endif
	
		for(i = 0; i < (int)numoftemplates; i++)
		{
			sequence = 0;	
			for ( k=0; k< m; k++ ) {
				fscanf(fp, "%d", &bit);
				//printf("%i",bit);
			//delete
			/*
#ifdef FILE_OUTPUT
				fprintf(stats[TEST_NONPERIODIC], "%d", bit);
#endif
			*/
				sequence <<= 1;
				sequence ^= bit;
			}
			templates[i] = Mirrored_int(sequence,m); //sequence;
			//printf("%d ",sequence%296);
		}
		//printf("\n\n\n");
		
		//delete
		/*
#ifdef FILE_OUTPUT
			fprintf(stats[TEST_NONPERIODIC], " ");
#endif
			*/
	
		//timings();
		for(block = 0; block < N; block++)
		{
			for(bit_ind = block*M; bit_ind < (block+1)*M-m+1; bit_ind++)
			{
				window = get_nth_block4(array,bit_ind);
				//bits(&window,4);
				Wj[block][window & mask]++;
				//printf("%d",window & mask);
			}
		}
		
		sum = 0.0;
		for ( i=0; i<2; i++ ) {                      /* Compute Probabilities */
			pi[i] = exp(-lambda+i*log(lambda)-cephes_lgam(i+1));
			sum += pi[i];
		}
		
		pi[0] = sum;
		for ( i=2; i<=K; i++ ) {                      /* Compute Probabilities */
			pi[i-1] = exp(-lambda+i*log(lambda)-cephes_lgam(i+1));
			sum += pi[i-1];
		}
		pi[K] = 1 - sum;
		
		///new
#ifdef FILE_OUTPUT
			fclose(fp);
			fp = fopen(directory, "r");
#endif
		///
		for(j = 0; j/SKIP < (int)MIN(numoftemplates,MAXNUMOFTEMPLATES); j += SKIP)
		{
			///new
#ifdef FILE_OUTPUT
			for ( k=0; k< m; k++ ) {
				fscanf(fp, "%d", &bit);
				fprintf(stats[TEST_NONPERIODIC], "%d", bit);
			
			}
			fprintf(stats[TEST_NONPERIODIC], " ");
#endif
			///
			one_template = templates[j];
			//printf("%3d ",one_template);
			chi2 = 0.0;                                   /* Compute Chi Square */
			for ( i=0; i<N; i++ ) {
				chi2 += pow(((double)Wj[i][one_template] - lambda)/pow(varWj, 0.5), 2);

#ifdef FILE_OUTPUT
				if ( m == 10 )
					fprintf(stats[TEST_NONPERIODIC], "%3d  ", Wj[i][one_template]);
				else
					fprintf(stats[TEST_NONPERIODIC], "%4d ", Wj[i][one_template]);
#endif
				/*if ( m == 10 )
					printf("%3d  ", Wj[i][one_template]);
				else
					printf("%4d ", Wj[i][one_template]);*/

#ifdef VERIFY_RESULTS
				//if((j/SKIP)>=R2.nonoverlapping.templates) printf("!!!\n");
				R2.nonoverlapping.W[(j/SKIP)*N+i]=Wj[i][one_template];
#endif
			
			}
			p_value = cephes_igamc(N/2.0, chi2/2.0);
			//printf("(2) chi2:%lf value: %lf\n",chi2,p_value);
			
#ifdef VERIFY_RESULTS
			//if((j/SKIP)>=R2.nonoverlapping.templates) printf("!!!\n");
			R2.nonoverlapping.chi2[j/SKIP]=chi2;
			R2.nonoverlapping.p_value[j/SKIP]=p_value;
#endif

#ifdef FILE_OUTPUT
			if ( isNegative(p_value) || isGreaterThanOne(p_value) )
				fprintf(stats[TEST_NONPERIODIC], "\t\tWARNING:  P_VALUE IS OUT OF RANGE.\n");

			fprintf(stats[TEST_NONPERIODIC], "%9.6f %f %s %3d\n", chi2, p_value, p_value < ALPHA ? "FAILURE" : "SUCCESS", j/SKIP);
#endif

#ifdef FILE_OUTPUT
			fprintf(results[TEST_NONPERIODIC], "%f\n", p_value); fflush(results[TEST_NONPERIODIC]);
#endif

		}
		
	}

#ifdef FILE_OUTPUT
	fprintf(stats[TEST_NONPERIODIC], "\n"); fflush(stats[TEST_NONPERIODIC]);
#endif

	for(i = 0; i < N; i++)
	{
		if(Wj[i])free(Wj[i]);
	}
	free(Wj);
	//TODO free(Wj);
    if(fp)fclose(fp); //PADA PRI M == 16
	
}
