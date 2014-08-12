#include <stdio.h>
#include "../include/tools.h"


void bits(unsigned char*array, int byte_size)
{
	int i,j;
	unsigned char val;
	for(i = 0; i < byte_size; i++)
	{
		val = array[i];
		for(j = 0; j < 8; j++)
		{
			printf("%d",val&1);
			val >>= 1;
		}
		printf(" ");
	}
	printf("\n");
}

unsigned int get_nth_block4(unsigned char* array, int offset)
{
	return (*(unsigned int*)(array+(offset>>3))) >> (offset & 7);//(array2[(offset >> 3)&3][(offset >> 3)] >> (offset & 7));
}
unsigned int get_nth_block_effect(unsigned char* array, int offset)
{	
	int shift = (offset & 7);
	int byte = (offset >> 3);
	if (shift == 0) return (*(unsigned int*)(array + byte) >> shift);
	else return (*(unsigned int*)(array + byte) >> shift)^(*(unsigned int*)(array + byte+4) << (32 - shift));
}

/*
void test_blocks(unsigned char* array, int size){
	int i,sum = 0,value = 0,counter = 7;
	unsigned char *copy_array = array - 1;
	

	timings();
	for(i = 0; i < size; i++)
	{
		
		sum  += get_nth_block4(array,i);
		
	}
	timings();
	printf("%d",sum);
}
*/

int Mirrored_int(unsigned int val, int m){
	int res = 0,i;
	for(i=0; i < m; i++)
	{
		if(val & (1 << i)) res += (1 << (m - 1 - i));
	}
	return res;
}
/*
int test(unsigned char*array, int n)
{
	int sum, i;

	for(i = 0; i < n; i++)
	{
		sum += get_nth_block4(array,i);
	}
	return sum;
}
*/