#ifndef _BM_H_
#define _BM_H_
#include "../include/tools.h"


/*typedef unsigned int type;


typedef struct mybitset{
	type *array;
	int first, last;
	int size,array_size;
} mybitset;



void set_bit(mybitset* bitset,int index);
void clear(mybitset* bitset);
void resize(mybitset* bitset, int new_size);
void left_shift(mybitset* bitset, int bits, int from, int to);

void XOR(mybitset* a, mybitset* b, int from, int to);
void copy(mybitset* a, mybitset* b, int from, int to);

void rand_array(mybitset* bitset);
int BM(mybitset bitstream, int N);
*/
typedef unsigned int type;
void print_bits(unsigned char c);
void array_as_bits(unsigned char* c,int byte_size);

void left_shift(type* array, int array_size, int from, int to);
void copy(type* array, type* b, int array_size, int from, int to);
void XOR(type* array, type* b, int array_size, int from, int to);
int BM_c(type* bitstream,int N,type* c, type *b,type* t);

#endif