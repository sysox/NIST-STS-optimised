#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <time.h>
#include <assert.h>
//#include "../include/decls.h"
#include "../include/externs.h"
#include "../include/cephes.h"  
#include "../include/utilities.h"
#include "../include/tools.h"
#include "../include/stat_fncs.h"

#ifdef _WIN32
#define M_LOG2E 1.44269504088896340736 //log2(e)
long double log2(const long double x){
    return  log(x) * M_LOG2E;
}
#endif


// Temp hack!!!
int n;

void timings(){
	
	static clock_t start = 0,end;
	static FILE *times;
	if(times == NULL)times=fopen("times.txt","w");
	//fprintf(times,"time");
	end = clock();
	if(start != 0)fprintf(times,"time %lf \n ",1.0*(end - start)/CLOCKS_PER_SEC);
	printf("time %lf \n",1.0*(end - start)/CLOCKS_PER_SEC);
	start = end;
}
unsigned char* load_array(FILE* f,int bit_size){
	unsigned char *array;
	unsigned int byte_size;

	if(bit_size < 1){
		fseek(f,0,SEEK_END);
		bit_size = ftell(f)*8;
		fseek (f,0,SEEK_SET);
	}
	n = bit_size;
	byte_size = (unsigned int)ceil(n/8)+4; 
	array = (unsigned char*)malloc(sizeof( unsigned char)*(byte_size));
	if(array==NULL) { printf("Cannot allocate memory.\n"); exit(10); }
	fread(array,1,byte_size,f);
	//after_byte = &array[byte_size]+1;
	array[byte_size-1]=array[byte_size-2]=array[byte_size-3]=array[byte_size-4]=0;
	return array;
}
void transform(BitSequence	*epsilon, unsigned char* array, int n){
	int i;

	for(i = 0; i < n; i++)
	{
		//printf("%d %d %d\n",array[i >> 3],(1 << (i & 7)), array[i >> 3] & (1 << (i & 7)));
		if( (array[i >> 3] & (1 << (i & 7))) != 0)epsilon[i]=1;
		else epsilon[i] = 0;
	}
}

void data_all_zeros(int bit_size)
{
	int i;
	int byte_size=bit_size/8;
	if(byte_size*8<bit_size)byte_size++;
	epsilon = (unsigned char*)malloc(sizeof(unsigned char)*(byte_size*8));
	array = (unsigned char*)malloc(sizeof( unsigned char)*(byte_size)+4);
	if(epsilon == NULL || array==NULL) { printf("Cannot allocate memory.\n"); exit(10); }
	for (i=0;i<byte_size;i++) array[i]=0;
	transform(epsilon,array,bit_size);
}

void data_all_ones_padded(int bit_size)
{
	int i;
	int byte_size=bit_size/8;
	if(byte_size*8<bit_size)byte_size++;
	epsilon = (unsigned char*)malloc(sizeof(unsigned char)*(byte_size*8));
	array = (unsigned char*)malloc(sizeof( unsigned char)*(byte_size)+4);
	if(epsilon == NULL || array==NULL) { printf("Cannot allocate memory.\n"); exit(10); }
	for (i=0;i<byte_size;i++) array[i]=0xff;
	transform(epsilon,array,bit_size);
}

void data_all_ones(int bit_size)
{
	int i, j, pad;
	unsigned char mask=0x80;
	unsigned char padding=0x00;
	int byte_size=bit_size/8;
	if(byte_size*8<bit_size)byte_size++;
	pad = byte_size*8 - bit_size;
	epsilon = (unsigned char*)malloc(sizeof(unsigned char)*(byte_size*8));
	array = (unsigned char*)malloc(sizeof( unsigned char)*(byte_size)+4);
	if(epsilon == NULL || array==NULL) { printf("Cannot allocate memory.\n"); exit(10); }
	for (i=0;i<byte_size;i++) array[i]=0xff;

	for(j=0;j<pad;j++)
	{
		padding+=mask;
		mask >>=1;
	}

	array[i-1]^=padding;
	transform(epsilon,array,bit_size);
}

void data_all_zeros_padded(int bit_size)
{
	int i, j, pad;
	unsigned char mask=0x80;
	unsigned char padding=0x00;
	int byte_size=bit_size/8;
	if(byte_size*8<bit_size)byte_size++;
	pad = byte_size*8 - bit_size;
	epsilon = (unsigned char*)malloc(sizeof(unsigned char)*(byte_size*8));
	array = (unsigned char*)malloc(sizeof( unsigned char)*(byte_size)+4);
	if(epsilon == NULL || array==NULL) { printf("Cannot allocate memory.\n"); exit(10); }
	for (i=0;i<byte_size;i++) array[i]=0;

	for(j=0;j<pad;j++)
	{
		padding+=mask;
		mask >>=1;
	}

	array[i-1]|=padding;
	transform(epsilon,array,bit_size);
}

void data_pattern_01_padded(int bit_size)
{
	int i, j, pad;
	unsigned char mask=0x80;
	unsigned char padding=0x00;
	int byte_size=bit_size/8;
	if(byte_size*8<bit_size)byte_size++;
	pad = byte_size*8 - bit_size;
	epsilon = (unsigned char*)malloc(sizeof(unsigned char)*(byte_size*8));
	array = (unsigned char*)malloc(sizeof( unsigned char)*(byte_size)+4);
	if(epsilon == NULL || array==NULL) { printf("Cannot allocate memory.\n"); exit(10); }
	for (i=0;i<byte_size;i++) array[i]=0xAA;

	for(j=0;j<pad;j++)
	{
		padding+=mask;
		mask >>=1;
	}

	array[i-1]|=padding;
	transform(epsilon,array,bit_size);
}

void data_pattern_01(int bit_size)
{
	int i, j, pad;
	unsigned char mask=0x80;
	unsigned char padding=0x00;
	int byte_size=bit_size/8;
	if(byte_size*8<bit_size)byte_size++;
	pad = byte_size*8 - bit_size;
	epsilon = (unsigned char*)malloc(sizeof(unsigned char)*(byte_size*8));
	array = (unsigned char*)malloc(sizeof( unsigned char)*(byte_size)+4);
	if(epsilon == NULL || array==NULL) { printf("Cannot allocate memory.\n"); exit(10); }
	for (i=0;i<byte_size;i++) array[i]=0xAA;

	for(j=0;j<pad;j++)
	{
		padding+=mask;
		mask >>=1;
	}

	padding=0xFF^padding;
	array[i-1]&=padding;
	transform(epsilon,array,bit_size);
}

void data_pattern_10_padded(int bit_size)
{
	int i, j, pad;
	unsigned char mask=0x80;
	unsigned char padding=0x00;
	int byte_size=bit_size/8;
	if(byte_size*8<bit_size)byte_size++;
	pad = byte_size*8 - bit_size;
	epsilon = (unsigned char*)malloc(sizeof(unsigned char)*(byte_size*8));
	array = (unsigned char*)malloc(sizeof( unsigned char)*(byte_size)+4);
	if(epsilon == NULL || array==NULL) { printf("Cannot allocate memory.\n"); exit(10); }
	for (i=0;i<byte_size;i++) array[i]=0x55;

	for(j=0;j<pad;j++)
	{
		padding+=mask;
		mask >>=1;
	}

	array[i-1]|=padding;
	transform(epsilon,array,bit_size);
}

void data_pattern_10(int bit_size)
{
	int i, j, pad;
	unsigned char mask=0x80;
	unsigned char padding=0x00;
	int byte_size=bit_size/8;
	if(byte_size*8<bit_size)byte_size++;
	pad = byte_size*8 - bit_size;
	epsilon = (unsigned char*)malloc(sizeof(unsigned char)*(byte_size*8));
	array = (unsigned char*)malloc(sizeof( unsigned char)*(byte_size)+4);
	if(epsilon == NULL || array==NULL) { printf("Cannot allocate memory.\n"); exit(10); }
	for (i=0;i<byte_size;i++) array[i]=0x55;

	for(j=0;j<pad;j++)
	{
		padding+=mask;
		mask >>=1;
	}

	padding=0xFF^padding;
	array[i-1]&=padding;
	transform(epsilon,array,bit_size);
}

void data_prandom_padded(int bit_size)
{
	int i, j, pad;
	unsigned char mask=0x80;
	unsigned char padding=0x00;
	int byte_size=bit_size/8;
	if(byte_size*8<bit_size)byte_size++;
	pad = byte_size*8 - bit_size;
	srand((unsigned)0/*time(NULL)*/);
	epsilon = (unsigned char*)malloc(sizeof(unsigned char)*(byte_size*8));
	array = (unsigned char*)malloc(sizeof( unsigned char)*(byte_size)+4);
	if(epsilon == NULL || array==NULL) { printf("Cannot allocate memory.\n"); exit(10); }
	for (i=0;i<byte_size;i++) array[i]=(0xFF&((unsigned)rand()));

	for(j=0;j<pad;j++)
	{
		padding+=mask;
		mask >>=1;
	}

	array[i-1]|=padding;
	transform(epsilon,array,bit_size);
}

void data_prandom(int bit_size)
{
	int i, j, pad;
	unsigned char mask=0x80;
	unsigned char padding=0x00;
	int byte_size=bit_size/8;
	if(byte_size*8<bit_size)byte_size++;
	pad = byte_size*8 - bit_size;
	srand((unsigned)0/*time(NULL)*/);
	epsilon = (unsigned char*)malloc(sizeof(unsigned char)*(byte_size*8));
	array = (unsigned char*)malloc(sizeof( unsigned char)*(byte_size)+4);
	if(epsilon == NULL || array==NULL) { printf("Cannot allocate memory.\n"); exit(10); }
	for (i=0;i<byte_size;i++) array[i]=(0xFF&((unsigned)rand()));

	for(j=0;j<pad;j++)
	{
		padding+=mask;
		mask >>=1;
	}

	padding=0xFF^padding;
	array[i-1]&=padding;
	transform(epsilon,array,bit_size);
}

void data_prandom_padded_random(int bit_size)
{
	int i, pad;
	unsigned char mask=0x80;
	unsigned char padding=0x00;
	int byte_size=bit_size/8;
	if(byte_size*8<bit_size)byte_size++;
	pad = byte_size*8 - bit_size;
	srand((unsigned)time(NULL));
	epsilon = (unsigned char*)malloc(sizeof(unsigned char)*(byte_size*8));
	array = (unsigned char*)malloc(sizeof( unsigned char)*(byte_size)+4);
	if(epsilon == NULL || array==NULL) { printf("Cannot allocate memory.\n"); exit(10); }
	for (i=0;i<byte_size;i++) array[i]=(0xFF&((unsigned)rand()));

	transform(epsilon,array,bit_size);
}

void data_bad_prandom_padded(int bit_size)
{
	int i, j, pad;
	unsigned char mask=0x80;
	unsigned char padding=0x00;
	int byte_size=bit_size/8;
	if(byte_size*8<bit_size)byte_size++;
	pad = byte_size*8 - bit_size;
	srand((unsigned)time(NULL));
	epsilon = (unsigned char*)malloc(sizeof(unsigned char)*(byte_size*8));
	array = (unsigned char*)malloc(sizeof( unsigned char)*(byte_size)+4);
	if(epsilon == NULL || array==NULL) { printf("Cannot allocate memory.\n"); exit(10); }
	for (i=0;i<byte_size;i++) array[i]=(0xFF&((unsigned)rand())) | (0xFF&((unsigned)rand()));

	for(j=0;j<pad;j++)
	{
		padding+=mask;
		mask >>=1;
	}

	array[i-1]|=padding;
	transform(epsilon,array,bit_size);
}

void data_bad_prandom(int bit_size)
{
	int i, j, pad;
	unsigned char mask=0x80;
	unsigned char padding=0x00;
	int byte_size=bit_size/8;
	if(byte_size*8<bit_size)byte_size++;
	pad = byte_size*8 - bit_size;
	srand((unsigned)time(NULL));
	epsilon = (unsigned char*)malloc(sizeof(unsigned char)*(byte_size*8));
	array = (unsigned char*)malloc(sizeof( unsigned char)*(byte_size)+4);
	if(epsilon == NULL || array==NULL) { printf("Cannot allocate memory.\n"); exit(10); }
	for (i=0;i<byte_size;i++) array[i]=(0xFF&((unsigned)rand()))  | (0xFF&((unsigned)rand())) ;

	for(j=0;j<pad;j++)
	{
		padding+=mask;
		mask >>=1;
	}

	padding=0xFF^padding;
	array[i-1]&=padding;
	transform(epsilon,array,bit_size);
}

void data_bad_prandom_padded_random(int bit_size)
{
	int i, pad;
	unsigned char mask=0x80;
	unsigned char padding=0x00;
	int byte_size=bit_size/8;
	if(byte_size*8<bit_size)byte_size++;
	pad = byte_size*8 - bit_size;
	srand((unsigned)time(NULL));
	epsilon = (unsigned char*)malloc(sizeof(unsigned char)*(byte_size*8));
	array = (unsigned char*)malloc(sizeof( unsigned char)*(byte_size)+4);
	if(epsilon == NULL || array==NULL) { printf("Cannot allocate memory.\n"); exit(10); }
	for (i=0;i<byte_size;i++) array[i]=(0xFF&((unsigned)rand()))  | (0xFF&((unsigned)rand()));

	transform(epsilon,array,bit_size);
}

#include <float.h>

union Double_Int
{
#ifdef _WIN32
    __int64 i;
#else
	int64_t i;
#endif
    double d;
};
 
int ReasonablyEqualDoubles(double D1, double D2, int Max)
{
    union Double_Int DD1;
    union Double_Int DD2;
#ifdef _WIN32
	__int64 Difference;
#else
	int64_t Difference;
#endif

	if(D1==D2) return 1;

	assert(sizeof(double) == 8);

	DD1.d=D1;
	DD2.d=D2;
	// Compare the sign
    if (((DD1.i >> 63) != 0) != ((DD2.i >> 63) != 0))
    {
        //if(D1 == D2) return 1;
        return 0;
    }
    // Compare the integers made from doubles
    Difference = DD1.i - DD2.i;
	if(Difference<0) Difference = -Difference;
    if(Difference<=Max) return 1; 
    return 0;
} 

#ifdef VERIFY_RESULTS
int equal_frequency_results()
{
	if(R1.frequency.sum!=R2.frequency.sum)return 0;
	if(R1.frequency.sum_n!=R2.frequency.sum_n)return 0;
	if(R1.frequency.p_value!=R2.frequency.p_value)return 0;
	return 1;
}

int equal_blockfrequency_results()
{
	if(R1.blockfrequency.chi_squared!=R2.blockfrequency.chi_squared)return 0;
	if(R1.blockfrequency.p_value!=R2.blockfrequency.p_value)return 0;
	return 1;
}



int equal_runs_results()
{
	if(!ReasonablyEqualDoubles(R1.runs.p_value,R2.runs.p_value,1)) return 0;
	//if(R1.runs.p_value!=R2.runs.p_value) 
	//{ printf("Comparing %lf and %lf\n",R1.runs.p_value,R2.runs.p_value);return 0;}
	if(R1.runs.V!=R2.runs.V) return 0;
	if(R1.runs.pi!=R2.runs.pi) return 0;
	//if(R1.runs.erfc_arg!=R2.runs.erfc_arg) return 0;
	//if(fabs(R1.runs.erfc_arg-R2.runs.erfc_arg)>FLT_EPSILON) return 0;
	if(!ReasonablyEqualDoubles(R1.runs.erfc_arg,R2.runs.erfc_arg,1)) return 0;
	printf("#");
	return 1;
}

int equal_longestrunofones_results()
{
	int i;
	if(R1.longestrunofones.pval!=R2.longestrunofones.pval)return 0;
	if(R1.longestrunofones.chi2!=R2.longestrunofones.chi2)return 0;
	if(R1.longestrunofones.N!=R2.longestrunofones.N)return 0;
	if(R1.longestrunofones.M!=R2.longestrunofones.M)return 0;
	for(i=0;i<7;i++) if(R1.longestrunofones.nu[i]!=R2.longestrunofones.nu[i])return 0;
	return 1;
}

int equal_rank_results()
{
	if(R1.rank.p_30!=R2.rank.p_30)return 0;
	if(R1.rank.p_31!=R2.rank.p_31)return 0;
	if(R1.rank.p_32!=R2.rank.p_32)return 0;
	if(R1.rank.F_30!=R2.rank.F_30)return 0;
	if(R1.rank.F_31!=R2.rank.F_31)return 0;
	if(R1.rank.F_32!=R2.rank.F_32)return 0;
	if(R1.rank.chi_squared!=R2.rank.chi_squared)return 0;
	if(R1.rank.p_value!=R2.rank.p_value)return 0;
	if(R1.rank.N!=R2.rank.N)return 0;
	printf("#");
	return 1;
}
int equal_serial_results()
{
	if(R1.serial.psim0!=R2.serial.psim0)return 0;
	if(R1.serial.psim1!=R2.serial.psim1)return 0;
	if(R1.serial.psim2!=R2.serial.psim2)return 0;
	if(R1.serial.del1!=R2.serial.del1)return 0;
	if(R1.serial.del2!=R2.serial.del2)return 0;
	if(R1.serial.p_value1!=R2.serial.p_value1)return 0;
	if(R1.serial.p_value2!=R2.serial.p_value2)return 0;	
	return 1;
}

int equal_nonoverlapping_results()
{
	unsigned int i,j;
	if(R1.nonoverlapping.templates!=R2.nonoverlapping.templates) 
		return 0;
	for(i=0;i<R1.nonoverlapping.templates;i++)
	{
		if(R1.nonoverlapping.chi2[i]!=R2.nonoverlapping.chi2[i]) 
			return 0;
		if(R1.nonoverlapping.p_value[i]!=R2.nonoverlapping.p_value[i]) return 0;
		for(j=0;j<8;j++)
			if(R1.nonoverlapping.W[i*8+j]!=R2.nonoverlapping.W[i*8+j]) return 0;
	}
	if(R1.nonoverlapping.chi2) free(R1.nonoverlapping.chi2);
	if(R1.nonoverlapping.p_value) free(R1.nonoverlapping.p_value);
	if(R1.nonoverlapping.W) free(R1.nonoverlapping.W);

	if(R2.nonoverlapping.chi2) free(R2.nonoverlapping.chi2);
	if(R2.nonoverlapping.p_value) free(R2.nonoverlapping.p_value);
	if(R2.nonoverlapping.W) free(R2.nonoverlapping.W);
	printf("#");
	return 1;
}

int equal_overlapping_results()
{
	int i;
	if(R1.overlapping.chi2==R1.overlapping.chi2)
	{
		if(R1.overlapping.chi2!=R2.overlapping.chi2) return 0;
	}
	else
	{
		if(R2.overlapping.chi2==R2.overlapping.chi2) return 0;
	}
	if(R1.overlapping.p_value==R1.overlapping.p_value)
	{
		if(R1.overlapping.p_value!=R2.overlapping.p_value) return 0;
	}
	else
	{
		if(R2.overlapping.p_value==R2.overlapping.p_value) return 0;
	}
	for(i=0;i<6;i++)
	{
		if(R1.overlapping.nu[i]!=R2.overlapping.nu[i]) return 0;
	}
	printf("#");
	return 1;
}

int equal_universal_results()
{
	if(R1.universal.phi!=R2.universal.phi) 
		return 0;
	if(R1.universal.p_value!=R2.universal.p_value) 
		return 0;
	if(R1.universal.sum!=R2.universal.sum) return 0;
	return 1;
}

int equal_apen_results()
{
	unsigned int i;
	if(R1.approximate_entropy.pp!=R2.approximate_entropy.pp) return 0;
	for(i=0;i<R1.approximate_entropy.pp;i++)
		if(R1.approximate_entropy.P[i]!=R2.approximate_entropy.P[i]) return 0;
	if(R1.approximate_entropy.P)free(R1.approximate_entropy.P);
	if(R2.approximate_entropy.P)free(R2.approximate_entropy.P);

	/*
	if(!ReasonablyEqualDoubles(R1.approximate_entropy.chi_squared, R2.approximate_entropy.chi_squared,16)) return 0;
	if(!ReasonablyEqualDoubles(R1.approximate_entropy.p_value, R2.approximate_entropy.p_value,16)) return 0;
	if(!ReasonablyEqualDoubles(R1.approximate_entropy.ApEn[0], R2.approximate_entropy.ApEn[0],16)) return 0;
	if(!ReasonablyEqualDoubles(R1.approximate_entropy.ApEn[1], R2.approximate_entropy.ApEn[1],16)) return 0;
	*/

	if(fabs(R1.approximate_entropy.chi_squared-R2.approximate_entropy.chi_squared)>FLT_EPSILON) return 0;
	if(fabs(R1.approximate_entropy.p_value-R2.approximate_entropy.p_value)>FLT_EPSILON) return 0;
	if(fabs(R1.approximate_entropy.ApEn[0]-R2.approximate_entropy.ApEn[0])>FLT_EPSILON) return 0;
	if(fabs(R1.approximate_entropy.ApEn[1]-R2.approximate_entropy.ApEn[1])>FLT_EPSILON) return 0;

	/*
	if(R1.approximate_entropy.chi_squared!=R2.approximate_entropy.chi_squared) 
		return 0;
	if(R1.approximate_entropy.p_value!=R2.approximate_entropy.p_value) 
		return 0;
	if(R1.approximate_entropy.ApEn[0]!=R2.approximate_entropy.ApEn[0]) 
		return 0;
	if(R1.approximate_entropy.ApEn[1]!=R2.approximate_entropy.ApEn[1]) 
		return 0;
		*/
	printf("#");
	return 1;
}

int equal_cusum_results()
{
	if(R1.cusum.z!=R2.cusum.z) return 0;
	if(R1.cusum.zrev!=R2.cusum.zrev) return 0;
	if(R1.cusum.sum1A!=R2.cusum.sum1A) return 0;
	if(R1.cusum.sum1B!=R2.cusum.sum1B) return 0;
	if(R1.cusum.sum2A!=R2.cusum.sum2A) return 0;
	if(R1.cusum.sum2B!=R2.cusum.sum2B) return 0;
	if(R1.cusum.p_valueA!=R2.cusum.p_valueA) return 0;
	if(R1.cusum.p_valueB!=R2.cusum.p_valueB) return 0;
	return 1;
}

int equal_random_excursion_results()
{
	int i;
	// temp hack
	if(!R1.random_excursion.valid) return 1;
	//
	if(R1.random_excursion.valid!=R2.random_excursion.valid) return 0;
	if(R1.random_excursion.valid)
	{
		for(i=0;i<8;i++)
		{
			if(R1.random_excursion.p_value[i]!=R2.random_excursion.p_value[i]) return 0;
			if(R1.random_excursion.sum[i]!=R2.random_excursion.sum[i]) return 0;
			if(R1.random_excursion.J[i]!=R2.random_excursion.J[i]) return 0;
			if(R1.random_excursion.x[i]!=R2.random_excursion.x[i]) return 0;
		}
	}
	printf("#");
	return 1;
}

int equal_random_excursion_var_results()
{
	int i;
	// temp hack
	if(!R1.random_excursion_variant.valid) return 1;
	//
	if(R1.random_excursion_variant.valid!=R2.random_excursion_variant.valid) return 0;
	if(R1.random_excursion_variant.valid)
	{
		for(i=0;i<18;i++)
		{
			if(R1.random_excursion_variant.p_value[i]!=R2.random_excursion_variant.p_value[i]) 
				return 0;
			if(R1.random_excursion_variant.count[i]!=R2.random_excursion_variant.count[i]) return 0;
			if(R1.random_excursion_variant.x[i]!=R2.random_excursion_variant.x[i]) return 0;
		}
	}
	printf("#");
	return 1;
}

int equal_linear_complexity_results()
{
	int i;

	if(R1.linear_complexity.p_value!=R2.linear_complexity.p_value) return 0;
	//printf("[1chi2]: %f [2chi2]: %f\n",R1.linear_complexity.chi2, R2.linear_complexity.chi2);
	if(R1.linear_complexity.chi2!=R2.linear_complexity.chi2) return 0;
	for(i=0;i<7;i++)
	{
		if(R1.linear_complexity.nu[i]!=R2.linear_complexity.nu[i]) return 0;
	}
	printf("#");
	return 1;
}

int equal_dft_results()
{
	if(R1.dft.p_value!=R2.dft.p_value) return 0;
	// works also with NaN
	if(R1.dft.percentile==R1.dft.percentile)
	{
		if(R1.dft.percentile!=R2.dft.percentile) return 0;
	}
	else
	{
		if(R2.dft.percentile==R2.dft.percentile) return 0;
	}
	if(R1.dft.N_l!=R2.dft.N_l) return 0;
	if(R1.dft.N_o!=R2.dft.N_o) return 0;
	if(R1.dft.d!=R2.dft.d) return 0;
	return 1;
}

int compare_results(int what)
{
	int result=0;
	switch(what)
	{
	case TEST_FREQUENCY: result=equal_frequency_results(); break;
	case TEST_BLOCK_FREQUENCY: result=equal_blockfrequency_results(); break;
	case TEST_RUNS: result=equal_runs_results(); break;
	case TEST_LONGEST_RUN: result=equal_longestrunofones_results(); break;
	case TEST_RANK: result=equal_rank_results(); break;
	case TEST_SERIAL: result=equal_serial_results(); break;
	case TEST_NONPERIODIC: result=equal_nonoverlapping_results(); break;
	case TEST_OVERLAPPING: result=equal_overlapping_results(); break;
	case TEST_UNIVERSAL: result=equal_universal_results(); break;
	case TEST_APEN: result=equal_apen_results(); break;
	case TEST_CUSUM: result=equal_cusum_results(); break;
	case TEST_RND_EXCURSION: result=equal_random_excursion_results(); break;
	case TEST_RND_EXCURSION_VAR: result=equal_random_excursion_var_results(); break;
	case TEST_LINEARCOMPLEXITY: result=equal_linear_complexity_results(); break;
	case TEST_FFT: result=equal_dft_results(); break;
	}

	if(result)
		printf("");
	else
		printf("!");

	return result;
}

#define DATA_ALL_ZEROS 0
#define DATA_ALL_ZEROS_PADDED 1
#define DATA_ALL_ONES 2
#define DATA_ALL_ONES_PADDED 3
#define DATA_PATTERN_01 4
#define DATA_PATTERN_01_PADDED 5
#define DATA_PATTERN_10 6
#define DATA_PATTERN_10_PADDED 7
#define DATA_PRANDOM 8
#define DATA_PRANDOM_PADDED 9
#define DATA_PRANDOM_PADDED_RANDOM 10
#define DATA_BAD_PRANDOM 11
#define DATA_BAD_PRANDOM_PADDED 12
#define DATA_BAD_PRANDOM_PADDED_RANDOM 13
#define DATA_MAX_NUMBER 13

void prepare_data(int j,int i)
{
	switch(j)
	{
	case DATA_ALL_ZEROS: data_all_zeros(i); break;
	case DATA_ALL_ZEROS_PADDED: data_all_zeros_padded(i); break;
	case DATA_ALL_ONES: data_all_ones(i); break;
	case DATA_ALL_ONES_PADDED: data_all_ones_padded(i); break;
	case DATA_PATTERN_01: data_pattern_01(i); break;
	case DATA_PATTERN_01_PADDED: data_pattern_01_padded(i); break;
	case DATA_PATTERN_10: data_pattern_10(i); break;
	case DATA_PATTERN_10_PADDED: data_pattern_10_padded(i); break;
	case DATA_PRANDOM: data_prandom(i); break;
	case DATA_PRANDOM_PADDED: data_prandom_padded(i); break;
	case DATA_PRANDOM_PADDED_RANDOM: data_prandom_padded_random(i); break;
	case DATA_BAD_PRANDOM: data_bad_prandom(i); break;
	case DATA_BAD_PRANDOM_PADDED: data_bad_prandom_padded(i); break;
	case DATA_BAD_PRANDOM_PADDED_RANDOM: data_bad_prandom_padded_random(i); break;
	}
}


void test(int testcase,int smallnumbers)
{
	unsigned int i,j,k,r;
	unsigned int from, to;

	if(smallnumbers)
	{
		from=1;
		to=1000000;
	}
	else
	{
		from=1;
		to=1000000000;
	}

	for (i=from;i<=to;i=(smallnumbers)?i+1:i*10+rand()%8)
	{
		if(!smallnumbers || i%1000==0)
		{
			printf("Testing %s %li bit.\n", testNames[testcase], i);
			fflush(stdout);
		}

		for (j=0;j<=DATA_MAX_NUMBER;j++) 
		{
			prepare_data(j,i);
			
			if(testcase == TEST_FREQUENCY)
			{
                Frequency(i); Frequency2(i);
				r=compare_results(TEST_FREQUENCY);
				if(!r) 
				{
					printf("i: %i, j: %i\n",i,j);
					exit(1);	
				}
			}
			
			if(testcase == TEST_BLOCK_FREQUENCY)
			{
				for(k=8;k<=i/10;k++)
				{
					if(i>=k)
					{
						BlockFrequency(k,i); BlockFrequency2(k,i);
						r=compare_results(TEST_BLOCK_FREQUENCY);
						if(!r) 
						{
							printf("i: %i, j: %i\n",i,j);
							exit(1);
							
						}
					}
				}
			}

			if(testcase == TEST_CUSUM)
			{
				CumulativeSums(i); CumulativeSums2(i);
				r=compare_results(TEST_CUSUM);
				if(!r) 
				{
					printf("i: %i, j: %i\n",i,j);
					getchar();
					exit(1);
					
				}
			}
			
			if(testcase == TEST_RUNS)
			{
				Runs(i); Runs2(i);
				r=compare_results(TEST_RUNS);
				if(!r) 
				{
					printf("i: %i, j: %i\n",i,j);
					exit(1);
					
				}
			}
			
			if(testcase == TEST_LONGEST_RUN)
			{
				if(i>=128)
				{
					LongestRunOfOnes(i); LongestRunOfOnes2(i);
					r=compare_results(TEST_LONGEST_RUN);
					if(!r) 
					{
						printf("i: %i, j: %i\n",i,j);
						exit(1);
						
					}
				}
			}

			if(testcase == TEST_RANK)
			{
				if(i>32*32)
				{
					Rank(i); Rank2(i);
					r=compare_results(TEST_RANK);
					if(!r) 
					{
						printf("i: %i, j: %i\n",i,j);
						exit(1);
						
					}
				}
			}

			if(testcase == TEST_FFT)
			{
				DiscreteFourierTransform(i); DiscreteFourierTransform2(i);
				r=compare_results(TEST_FFT);
				if(!r) 
				{
					printf("i: %i, j: %i\n",i,j);
					exit(1);
					
				}
			}

			if(testcase == TEST_NONPERIODIC)
			{
				for(k=2;k<=15;k++) // 2-15
				{
					if(i>=k)
					{
						NonOverlappingTemplateMatchings(k,i); NonOverlappingTemplateMatchings2(k,i);
						r=compare_results(TEST_NONPERIODIC);
						if(!r) 
						{
							printf("i: %i, j: %i, m: %i\n",i,j,k);
							exit(1);
							
						}
					}
				}
			}


			if(testcase == TEST_OVERLAPPING)
			{
				for(k=9;k<=10;k++)
				{
					if(i>=k)
					{
						OverlappingTemplateMatchings(k,i); OverlappingTemplateMatchings2(k,i);
						r=compare_results(TEST_OVERLAPPING);
						if(!r) 
						{
							printf("i: %i, j: %i, m: %i\n",i,j,k);
							exit(1);
							
						}
					}
				}
			}

			if(testcase == TEST_UNIVERSAL)
			{
				Universal(i); Universal2(i);
				r=compare_results(TEST_UNIVERSAL);
				if(!r) 
				{
					printf("i: %i, j: %i\n",i,j);
					exit(1);
					
				}
			}

			if(testcase == TEST_APEN)
			{
				for(k=1;k<=23;k++) // 1-23
				{
					if(i>=k)
					{
						ApproximateEntropy(k,i); ApproximateEntropy2(k,i);
						r=compare_results(TEST_APEN);
						if(!r) 
						{
							printf("i: %i, j: %i, m: %i\n",i,j,k);
							exit(1);
						}
					}
				}
			}

			if(testcase == TEST_RND_EXCURSION)
			{
                RandomExcursions(i); RandomExcursions2(i);
				r=compare_results(TEST_RND_EXCURSION);
				if(!r) 
				{
					printf("i: %i, j: %i\n",i,j);
					exit(1);
					
				}
			}


			if(testcase == TEST_RND_EXCURSION_VAR)
			{
                RandomExcursionsVariant(i); RandomExcursionsVariant2(i);
				r=compare_results(TEST_RND_EXCURSION_VAR);
				if(!r) 
				{
					printf("i: %i, j: %i\n",i,j);
					exit(1);
					
				}
			}
		
			if(testcase == TEST_SERIAL)
			{
				for(k=1;k<=(log2(i)+1);k++)
				{
					if(i>=k)
					{
						Serial(k,i); Serial2(k,i);
						r=compare_results(TEST_SERIAL);
						if(!r) 
						{
							printf("i: %i, j: %i, k: %i\n",i,j,k);
							exit(1);
							
						}
					}
				}
			}


			if(testcase == TEST_LINEARCOMPLEXITY)
			{
				if(i>=1000000)
				{
					for(k=500;k<=5000;k++)
					{
						LinearComplexity(k,i); LinearComplexity2(k,i);
						r=compare_results(TEST_LINEARCOMPLEXITY);
						if(!r) 
						{
							printf("i: %i, j: %i, m: %i\n",i,j,k);
							exit(1);
							
						}
					}
				}
			}

			free(epsilon); free(array);
		}
	}
}
#endif 

#ifdef _WIN64
// In x64 Windows systems use __rdtsc function (asm is not available)
#include <intrin.h>
__int64 GetCpuClocks()
{
	return __rdtsc();
}
#else
#ifdef _WIN32
// In win32 systems use asm
__int64 GetCpuClocks()
{
    struct { __int32 low, high; } clocks;
    __asm push EAX
    __asm push EDX
    __asm __emit 0fh 
	__asm __emit 031h 
    __asm mov clocks.low, EAX
    __asm mov clocks.high, EDX
    __asm pop EDX
    __asm pop EAX
    return *(__int64 *)(&clocks);
}
#else
// On Linux use asm
#include <stdint.h>
uint64_t GetCpuClocks(){
    unsigned int low, high;
    __asm__ __volatile__("rdtsc":"=a"(low),"=d"(high));
    return ((uint64_t)high << 32) | low;
}
#endif
#endif

#ifdef _WIN32
#include <windows.h>
#else
#include <time.h>
clock_t GetTickCount()
{
	return clock();
}
#endif

#include <locale.h>

void speed(int scale,int repeat, int test_from, int test_to)
{
	int n,t,j,from,to;
	FILE *f;
#ifdef VERIFY_RESULTS
	int r;
#endif
#ifdef _WIN32
	__int64 Astart, Amiddle, Aend, Atime1, Atime2, Amintime1,Amintime2;
	unsigned long Bstart, Bmiddle, Bend, Btime1, Btime2;
//	LARGE_INTEGER Cstart, Cmiddle, Cend;
//	__int64 Ctime1, Ctime2, Dtime1, Dtime2;
//	FILETIME Dstart,Dmiddle,Dend,Dtest1,Dtest2,Dtest3,Dvoid;
	DWORD affinity;
#else
	int64_t Astart, Amiddle, Aend, Atime1, Atime2, Amintime1,Amintime2;
	clock_t Bstart, Bmiddle, Bend, Btime1, Btime2;
#endif

#ifdef _WIN32
	affinity=1;
	SetThreadAffinityMask(GetCurrentThread(),affinity);
	SetPriorityClass(GetCurrentProcess(),HIGH_PRIORITY_CLASS);
#endif

	setlocale(LC_NUMERIC,"Czech");
	f=fopen("speed.csv","wt");
	if(!f) {printf("Cannot open file speed.csv for writing.\n"); exit(3); } ;

	if(scale)
	{
		from=1024*1*8;
		to=1024*1024*8*100;
	}
	else
	{
		from=1024*1024*8*200;
		to=1024*1024*8*200;
	}

	for(n=from;n<=to;n=(scale)?n*10+rand()%8:n+1)
	{
        for(t=test_from;t<=test_to;t++)
		{
			data_prandom(n);

			for(j=1;j<=repeat;j++)
			{
				Astart=GetCpuClocks();
				Bstart=GetTickCount();
#ifdef _WIN32
				/*
				GetProcessTimes(GetCurrentProcess(),&Dvoid,&Dvoid,&Dtest1,&Dstart);
				QueryPerformanceCounter(&Cstart);
				*/
#endif

				switch(t)
				{
				case TEST_FREQUENCY: Frequency(n); break;
				case TEST_BLOCK_FREQUENCY: if(n>=20) BlockFrequency(20,n); break;
				case TEST_CUSUM: CumulativeSums(n); break;
				case TEST_RUNS: Runs(n); break;
				case TEST_LONGEST_RUN: LongestRunOfOnes(n); break;
				case TEST_RANK: if(n>32*32) Rank(n); break;
				case TEST_FFT: DiscreteFourierTransform(n); break;
				case TEST_NONPERIODIC: NonOverlappingTemplateMatchings(2,n); break;
				case TEST_OVERLAPPING: OverlappingTemplateMatchings(10,n); break;
				case TEST_UNIVERSAL: Universal(n); break;
				case TEST_APEN: ApproximateEntropy(2,n); break;
				case TEST_RND_EXCURSION: RandomExcursions(n); break;
				case TEST_RND_EXCURSION_VAR: RandomExcursionsVariant(n); break;
				case TEST_SERIAL: if(log2(n)>=8) Serial(2,n); break;
				case TEST_LINEARCOMPLEXITY: LinearComplexity(5000,n); break;
				}
				 
				Amiddle=GetCpuClocks();
				Bmiddle=GetTickCount();
#ifdef _WIN32
				/*
				QueryPerformanceCounter(&Cmiddle);
				GetProcessTimes(GetCurrentProcess(),&Dvoid,&Dvoid,&Dtest2,&Dmiddle);
				*/
#endif

				switch(t)
				{
				case TEST_FREQUENCY: Frequency2(n); break;
				case TEST_BLOCK_FREQUENCY: if(n>=20) BlockFrequency2(20,n); break;
				case TEST_CUSUM: CumulativeSums2(n); break;
				case TEST_RUNS: Runs2(n); break;
				case TEST_LONGEST_RUN: LongestRunOfOnes2(n); break;
				case TEST_RANK: if(n>32*32) Rank2(n); break;
				case TEST_FFT: DiscreteFourierTransform2(n);  break;
				case TEST_NONPERIODIC: NonOverlappingTemplateMatchings2(2,n); break;
				case TEST_OVERLAPPING: OverlappingTemplateMatchings2(10,n); break;
				case TEST_UNIVERSAL: Universal2(n); break;
				case TEST_APEN: ApproximateEntropy2(2,n); break;
				case TEST_RND_EXCURSION: RandomExcursions2(n); break;
				case TEST_RND_EXCURSION_VAR: RandomExcursionsVariant2(n); break;
				case TEST_SERIAL: if(log2(n)>=2) Serial2(2,n); break;
				case TEST_LINEARCOMPLEXITY: LinearComplexity2(5000,n);  break;
				}

				Aend=GetCpuClocks();
				Bend=GetTickCount();
#ifdef _WIN32
				/*
				QueryPerformanceCounter(&Cend);
				GetProcessTimes(GetCurrentProcess(),&Dvoid,&Dvoid,&Dtest3,&Dend);
				*/
#endif

				Atime1=Amiddle-Astart;
				Atime2=Aend-Amiddle;
				if(j==1)
				{
					Amintime1=Atime1;
					Amintime2=Atime2;
				}
				if(Atime1<Amintime1) Amintime1=Atime1;
				if(Atime2<Amintime2) Amintime2=Atime2;

				Btime1=Bmiddle-Bstart;
				Btime2=Bend-Bmiddle;

#ifdef _WIN32
				/*
				Ctime1=Cmiddle.QuadPart-Cstart.QuadPart;
				Ctime2=Cend.QuadPart-Cmiddle.QuadPart;
				Dtime1=((((__int64)Dmiddle.dwHighDateTime)<<32) | (__int64)(Dmiddle.dwLowDateTime)) - ((((__int64)Dstart.dwHighDateTime)<<32) | (__int64)Dstart.dwLowDateTime);
				Dtime2=((((__int64)Dend.dwHighDateTime)<<32) | (__int64)(Dend.dwLowDateTime)) - ((((__int64)Dmiddle.dwHighDateTime)<<32) | (__int64)Dmiddle.dwLowDateTime);
				*/
#endif

				//printf("%s (%i bits): %f x faster [CPU clocks], %f x faster [GetTickCount - %ul ms vs. %ul ms], %f x faster [QueryPerformanceCounter], %f x faster [GetProcessTimes]\n",testNames[t],n,(float)Atime1/Atime2, (float)Btime1/Btime2, Btime1, Btime2, (float)Ctime1/Ctime2, (float)Dtime1/Dtime2);
				printf("%s (%i bits): %f x faster [CPU clocks], %f x faster [GetTickCount - %lu ms vs. %lu ms].\n",testNames[t],n,(float)Atime1/Atime2, (float)Btime1/Btime2, Btime1, Btime2);
				//printf("%s (%i bits): %f x faster\n",testNames[t],n,(float)Atime1/Atime2);
				//printf("%s (%i bits): %f x faster {GetProcessTimes}\n",testNames[t],n,(float)Dtime1/Dtime2);
				fprintf(f, "%s;%I64i;%I64i;%f;%lu;%lu;%f\n", testNames[t], Atime1, Atime2, (float)Atime1 / Atime2, Btime1, Btime2, (float)Btime1 / Btime2);
				fflush(stdout);
				fflush(f);

#ifdef VERIFY_RESULTS
				r=compare_results(t);
				if(!r)
				{
					printf("%s (%i bits): %f x faster [CPU clocks], %f x faster [GetTickCount - %ul ms vs. %ul ms].\n",testNames[t],n,(float)Atime1/Atime2, (float)Btime1/Btime2, Btime1, Btime2);
					printf("Results DO NOT MATCH\n");
					//exit(2);
				}
#endif

				}
				//printf("MINIMUM - %s (%i bits): %f x faster [CPU clocks]\n",testNames[t],n,(float)Amintime1/Amintime2);
				free(epsilon); free(array);
		}
		printf("Done n=%i.\n",n);
	}
	fclose(f);
}

#ifdef VERIFY_RESULTS
int main(int argc, char **argv)
{
	int testcase;

	if(argc>=2)
		testcase=atoi(argv[1]);
	else
		testcase=1;

	test(testcase,0);
	test(testcase,1);

	return 0;
}
#endif

#ifdef SPEED
int main(int argc, char **argv)
{
	int scale,repeat, test_from, test_to;

	if(argc>=2)
		scale=atoi(argv[1]);
	else
		scale=0;
	if(argc>=3)
		repeat=atoi(argv[2]);
	else
		repeat=1;
	if (argc >= 4)
		test_from = atoi(argv[3]);
	else
		test_from = 1;
	if (argc >= 5)
		test_to = atoi(argv[4]);
	else
		test_to = 15;
	
	speed(scale,repeat,test_from,test_to);
	
	return 0;
}
#endif