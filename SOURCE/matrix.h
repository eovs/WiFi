#ifndef _MATRIX_H_
#define _MATRIX_H_
/*
extern int H1x2m27[12][24];
extern int H2x3m27[8][24];
extern int H3x4m27[6][24];
extern int H5x6m27[4][24];
extern int H1x2m54[12][24];
extern int H2x3m54[8][24];
extern int H3x4m54[6][24];
extern int H5x6m54[4][24];
extern int H1x2m81[12][24];
extern int H2x3m81[6][24];
extern int H3x4m81[6][24];
extern int H5x6m81[4][24];
*/
typedef struct
{
	int nrow;
	int ncol;
	int M;
	int *matrix;
} CODE_CFG;

CODE_CFG select_code( int rate, int M );
#endif//_MATRIX_H_
