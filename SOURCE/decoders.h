#ifndef _DECODERS_H_
#define _DECODERS_H_

#include "matrix.h"

#define IL_SOFT_FPP	1
#define ONE_IL_SOFT ( 1 << IL_SOFT_FPP )

#define MAX_STATE 10000		// max number of decoder instances

enum  DEC_ID
{
	LMS_DEC,
	ILMS_DEC,
	LCHE_DEC,
	ILCHE_DEC,
};

extern char const * const DEC_FULL_NAME[];


#define SKIP  -1


typedef unsigned int        ui32;
typedef unsigned short      ui16;
typedef short               i16;
typedef int                 i32;
typedef unsigned char       ui8;
typedef signed char         i8;


#if 01
typedef double MS_DATA;
#else
typedef char MS_DATA;
#endif

typedef int IMS_DATA;


typedef struct  
{
	MS_DATA min1; 
	MS_DATA min2; 
	int pos;
	int sign;
} MS_DEC_STATE;

typedef struct  
{
	IMS_DATA min1; 
	IMS_DATA min2; 
	int pos;
	int sign;

	int sum;
	int cnt;
} IMS_DEC_STATE;

typedef struct
{
	int iter;
	int dec_index;
}DEC_RES;

typedef struct
{
	int nh; 
	int rh;
	int m;
    int codelen;
	int synd_len;
	int maxiter;
    int codec_id;
	int   **hd; 
	int   *syndr;                     


	// Layered Min-Sum Decoder
	MS_DATA *lms_soft;
	int	*lms_BnNS;
	MS_DEC_STATE *lms_dcs;
	MS_DEC_STATE *lms_tmps;
	MS_DATA *lms_buffer;
	MS_DATA *lms_rbuffer;
	MS_DATA *lms_rsoft;

	// Integer Layered Min-Sum Decoder
	IMS_DATA *ilms_y;
	IMS_DATA *ilms_decword;
	IMS_DATA *ilms_soft;
	int	*ilms_BnNS;
	IMS_DEC_STATE *ilms_dcs;
	IMS_DEC_STATE *ilms_tmps;
	IMS_DATA *ilms_buffer;
	IMS_DATA *ilms_rbuffer;
	IMS_DATA *ilms_rsoft;

	// Low complexity-high efficienty Decoder
	double  *lche_data0;
	double	*lche_tmp;
	double  *lche_soft_out;
	double  **lche_state;

	// integer Low complexity-high efficienty Decoder
	int  *ilche_data0;
	int	 *ilche_tmp;
	int  *ilche_soft_out;
	int  **ilche_state;

}DEC_STATE;


// functions prototypes
DEC_STATE* decod_open( int decoder_id, int mh, int nh, int M );
int decod_init( void* st );
void decod_close( DEC_STATE* st );
//int min_sum_decod_qc_lm( DEC_STATE* st, int soft[], int decword[], int maxiter, double alpha );    
int lmin_sum_decod_qc_lm( DEC_STATE* st, int soft[], int decword[], int maxiter, double alpha, double beta, int pre_shift );    
int il_min_sum_decod_qc_lm( DEC_STATE* st, int soft[], int decword[], int maxiter, double alpha, double beta,  int inner_data_bits, int pre_shift );    
int lche_decod( DEC_STATE* st, int soft[], int decword[], int maxiter );    
int ilche_decod( DEC_STATE* st, int soft[], int decword[], int maxsteps );

double** Alloc2d_double( int b, int c );

int icheck_syndrome( int **matr, int rh, int nh, IMS_DATA *soft, IMS_DATA *rsoft, int m_ldpc, int *synd );
void il_min_sum_reset( DEC_STATE *st );
int il_min_sum_iterate( DEC_STATE* st, int inner_data_bits );

void open_ext_il_minsum( int irate, int M, int num );
CODE_CFG open_ext_il_minsum( char *file_name, int M );
DEC_RES ext_il_min_sum( int *dec_input, int *dec_output, int n_iter, double alpha, double beta, int inner_data_bits );
void close_ext_il_minsum( void );


#endif	//_DECODERS_H_
