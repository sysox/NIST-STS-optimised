/* --------------------------------------------------------------------------

The following code is distributed under the following BSD-style license:

Copyright © 2013-2014 Marek Sys (syso@fi.muni.cz) & Zdenek Riha (zriha@fi.muni.cz).
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

unsigned int get_block_fast(unsigned char* array, int byte_offset)
{
	unsigned int res = (*((unsigned int*)(array + byte_offset)));

	return res;
}

unsigned int get_2bytes(unsigned char* array, int byte_offset)
{

	return array[byte_offset] ^ (array[byte_offset + 1] << 8);
}

unsigned int get_mask(int size){
	return (1 << size) - 1;
}
