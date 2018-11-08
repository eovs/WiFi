#ifndef _MATRIX_H_
#define _MATRIX_H_

typedef struct
{
	int nrow;
	int ncol;
	int M;
	int *matrix;
} CODE_CFG;

CODE_CFG select_code( int rate, int M );
#endif//_MATRIX_H_
