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


char *test_name(int t)
{
	if (t <= 15) return testNames[t];
	switch (t)
	{
	case 16: 
	case 29: return "BlockFrequency";
	case 17:
	case 18: 
	case 30: return "NonOverlappingTemplate";
	case 19:
	case 20: 
	case 31: return "OverlappingTemplate";
	case 21:
	case 22:
	case 23: 
	case 34:
	case 35:
	case 36: return "ApproximateEntropy";
	case 24:
	case 25: 
	case 32:
	case 26:
	case 33: return "Serial";
	case 27:
	case 28: return "LinearComplexity";
	}
	return "Unknown";
}

#include <locale.h>

void speed(int scale,int repeat, int test_from, int test_to)
{
	int n,t,j,from,to,param;
	FILE *f;

#ifdef _WIN32
	__int64 Astart, Amiddle, Aend, Atime1, Atime2, Amintime1,Amintime2;
	unsigned long Bstart, Bmiddle, Bend, Btime1, Btime2, Bmintime1, Bmintime2;
//	LARGE_INTEGER Cstart, Cmiddle, Cend;
//	__int64 Ctime1, Ctime2, Dtime1, Dtime2;
//	FILETIME Dstart,Dmiddle,Dend,Dtest1,Dtest2,Dtest3,Dvoid;
	DWORD affinity;
#else
	int64_t Astart, Amiddle, Aend, Atime1, Atime2, Amintime1,Amintime2;
	clock_t Bstart, Bmiddle, Bend, Btime1, Btime2, Bmintime1, Bmintime2;
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
		from=1024*1024*8*20;
		to=1024*1024*8*20;
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

				param = 0;

				switch(t)
				{
				case TEST_FREQUENCY: Frequency_v1(n); break;
				case TEST_BLOCK_FREQUENCY: if (n >= 100) BlockFrequency_v1(n / 100, n); param = n / 100;  break;
				case TEST_CUSUM: CumulativeSums_v1(n); break;
				case TEST_RUNS: Runs_v1(n); break;
				case TEST_LONGEST_RUN: LongestRunOfOnes_v1(n); break;
				case TEST_RANK: if(n>32*32) Rank_v1(n); break;
				case TEST_FFT: DiscreteFourierTransform_v1(n); break;
				case TEST_NONPERIODIC: NonOverlappingTemplateMatchings_v1(10, n); param = 10;  break;
				case TEST_OVERLAPPING: OverlappingTemplateMatchings_v1(10, n); param = 10;  break;
				case TEST_UNIVERSAL: Universal_v1(n); break;
				case TEST_APEN: ApproximateEntropy_v1(9, n); param = 9;  break;
				case TEST_RND_EXCURSION: RandomExcursions_v1(n); break;
				case TEST_RND_EXCURSION_VAR: RandomExcursionsVariant_v1(n); break;
				case TEST_SERIAL: if (log2(n) >= 9) Serial_v1(9, n); param = 9;  break;
				case TEST_LINEARCOMPLEXITY: LinearComplexity_v1(5000, n); param = 5000;  break;
				case 16: BlockFrequency_v1(20, n); param = 20; break;
				case 17: NonOverlappingTemplateMatchings_v1(2, n); param = 2; break;
				case 18: NonOverlappingTemplateMatchings_v1(9, n); param = 9; break;
				case 19: OverlappingTemplateMatchings_v1(2, n); param = 2;  break;
				case 20: OverlappingTemplateMatchings_v1(9, n); param = 9;  break;
				case 21: ApproximateEntropy_v1(2, n); param = 2;  break; 
				case 22: ApproximateEntropy_v1(5, n); param = 5;  break;
				case 23: ApproximateEntropy_v1(24, n); param = 24;  break;
				case 24: if (log2(n) >= 2) Serial_v1(2, n); param = 2;  break;
				case 25: if (log2(n) >= 24) Serial_v1(24, n); param = 24;  break;
				case 26: if (log2(n) >= 5) Serial_v1(5, n); param = 5;  break;
				case 27: LinearComplexity_v1(500, n); param = 500;   break;
				case 28: LinearComplexity_v1(1000, n); param = 1000;   break;
				case 29: BlockFrequency_v1(128, n); param = 128; break;
				case 30: NonOverlappingTemplateMatchings_v1(21, n); param = 21; break;
				case 31: OverlappingTemplateMatchings_v1(24, n); param = 24;  break;
				case 32: if (log2(n) >= 15) Serial_v1(13, n); param = 13;  break;
				case 33: if (log2(n) >= 17) Serial_v1(14, n); param = 14;  break;
				case 34: ApproximateEntropy_v1(27, n); param = 27;  break;
				case 35: ApproximateEntropy_v1(8, n); param = 8;  break;
				case 36: ApproximateEntropy_v1(10, n); param = 10;  break;
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
				case TEST_FREQUENCY: Frequency_v2(n); break;
				case TEST_BLOCK_FREQUENCY: if (n >= 100) BlockFrequency_v2(n/100, n); break;
				case TEST_CUSUM: CumulativeSums_v2(n); break;
				case TEST_RUNS: Runs_v2(n); break;
				case TEST_LONGEST_RUN: LongestRunOfOnes_v2(n); break;
				case TEST_RANK: if(n>32*32) Rank_v2(n); break;
				case TEST_FFT: DiscreteFourierTransform_v2(n);  break;
				case TEST_NONPERIODIC: NonOverlappingTemplateMatchings_v2(10,n); break;
				case TEST_OVERLAPPING: OverlappingTemplateMatchings_v2(10,n); break;
				case TEST_UNIVERSAL: Universal_v2(n); break;
				case TEST_APEN: ApproximateEntropy_v2(9,n); break;
				case TEST_RND_EXCURSION: RandomExcursions_v2(n); break;
				case TEST_RND_EXCURSION_VAR: RandomExcursionsVariant_v2(n); break;
				case TEST_SERIAL: if(log2(n)>=9) Serial_v2(9,n); break;
				case TEST_LINEARCOMPLEXITY: LinearComplexity_v2(5000,n);  break;
				case 16: BlockFrequency_v2(20, n); break;
				case 17: NonOverlappingTemplateMatchings_v2(2, n); break;
				case 18: NonOverlappingTemplateMatchings_v2(9, n); break;
				case 19: OverlappingTemplateMatchings_v2(2, n); break;
				case 20: OverlappingTemplateMatchings_v2(9, n); break;
				case 21: ApproximateEntropy_v2(2, n); break;
				case 22: ApproximateEntropy_v2(5, n); break;
				case 23: ApproximateEntropy_v2(24, n); break;
				case 24: if (log2(n) >= 2) Serial_v2(2, n); break;
				case 25: if (log2(n) >= 24) Serial_v2(24, n); break;
				case 26: if (log2(n) >= 5) Serial_v2(5, n); break;
				case 27: LinearComplexity_v2(500, n); break;
				case 28: LinearComplexity_v2(1000, n); break;
				case 29: BlockFrequency_v2(128, n); break;
				case 30: NonOverlappingTemplateMatchings_v2(21, n); break;
				case 31: OverlappingTemplateMatchings_v2(24, n); break;
				case 32: if (log2(n) >= 15) Serial_v2(13, n); break;
				case 33: if (log2(n) >= 17) Serial_v2(14, n); break;
				case 34: ApproximateEntropy_v2(27, n); break;
				case 35: ApproximateEntropy_v2(8, n); break;
				case 36: ApproximateEntropy_v2(10, n); break;
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
				Btime1 = Bmiddle - Bstart;
				Btime2 = Bend - Bmiddle;

				if(j==1)
				{
					Amintime1=Atime1;
					Amintime2=Atime2;
					Bmintime1 = Btime1;
					Bmintime2 = Btime2;
				}

				if(Atime1<Amintime1) Amintime1=Atime1;
				if(Atime2<Amintime2) Amintime2=Atime2;
				if (Btime1<Bmintime1) Bmintime1 = Btime1;
				if (Btime2<Bmintime2) Bmintime2 = Btime2;


#ifdef _WIN32
				/*
				Ctime1=Cmiddle.QuadPart-Cstart.QuadPart;
				Ctime2=Cend.QuadPart-Cmiddle.QuadPart;
				Dtime1=((((__int64)Dmiddle.dwHighDateTime)<<32) | (__int64)(Dmiddle.dwLowDateTime)) - ((((__int64)Dstart.dwHighDateTime)<<32) | (__int64)Dstart.dwLowDateTime);
				Dtime2=((((__int64)Dend.dwHighDateTime)<<32) | (__int64)(Dend.dwLowDateTime)) - ((((__int64)Dmiddle.dwHighDateTime)<<32) | (__int64)Dmiddle.dwLowDateTime);
				*/
#endif

				//printf("%s (%i bits): %f x faster [CPU clocks], %f x faster [GetTickCount - %ul ms vs. %ul ms], %f x faster [QueryPerformanceCounter], %f x faster [GetProcessTimes]\n",testNames[t],n,(float)Atime1/Atime2, (float)Btime1/Btime2, Btime1, Btime2, (float)Ctime1/Ctime2, (float)Dtime1/Dtime2);
				printf("%s [%i] (%i bits): %f x faster [CPU clocks], %f x faster [GetTickCount - %lu ms vs. %lu ms].\n", test_name(t), param,n, (float)Atime1 / Atime2, (float)Btime1 / Btime2, Btime1, Btime2);
				//printf("%s (%i bits): %f x faster\n",testNames[t],n,(float)Atime1/Atime2);
				//printf("%s (%i bits): %f x faster {GetProcessTimes}\n",testNames[t],n,(float)Dtime1/Dtime2);
				//fprintf(f, "%s [%i];%I64i;%I64i;%f;%lu;%lu;%f\n", test_name(t),param, Atime1, Atime2, (float)Atime1 / Atime2, Btime1, Btime2, (float)Btime1 / Btime2);
				fflush(stdout);
				//fflush(f);
				}
				printf("MINIMUM - %s (%i bits) [%i]: %f x faster [CPU clocks: %I64i vs. %I64i], %f x faster [ms: %i vs. %i]\n", test_name(t), n, param, (float)Amintime1 / Amintime2, Amintime1, Amintime2, (float)Bmintime1 / Bmintime2, Bmintime1, Bmintime2);
				fprintf(f, "%s [%i];%I64i;%I64i;%f;%lu;%lu;%f\n", test_name(t), param, Amintime1, Amintime2, (float)Amintime1 / Amintime2, Bmintime1, Bmintime2, (float)Bmintime1 / Bmintime2);
				fflush(f);
				free(epsilon); free(array);
		}
		printf("Done n=%i.\n",n);
	}
	fclose(f);
}




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
		repeat=10;
	if (argc >= 4)
		test_from = atoi(argv[3]);
	else
		test_from = 1;
	if (argc >= 5)
		test_to = atoi(argv[4]);
	else
		test_to = 36;
	
	speed(scale,repeat,test_from,test_to);
	
	return 0;
}
#endif
