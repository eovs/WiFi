#ifndef _DECODERS_H_
#define _DECODERS_H_


//#define ORIG_TABLES

#define DEC_DECISION 0	// 0 - hard decision, 1 - soft decision

//#define BP_USE_EPS	// be careful!!!

#define IASP_FIXED_POINT

#define MS_MUL_CORRECTION
#define MS_ALPHA_FPP 4

enum  DEC_ID
{
	BP_DEC, 
	SP_DEC, 
	ASP_DEC, 
	MS_DEC,
	IMS_DEC,
	IASP_DEC,
	FHT_DEC,
	TASP_DEC,
	LMS_DEC,
	LCHE_DEC,
	IL_MS_DEC,
	ILCHE_DEC,

};

extern char const * const DEC_FULL_NAME[];


#define SKIP  -1
#define TRUE_CIRCULANT

/*
#ifdef MS_MUL_CORRECTION
#define MS_ALPHA 0.8
#else
#define MS_ALPHA 0.02//0.075//0.02
#endif
*/
#define MS_ALPHA 0.8//0.8
#define MS_BETA  0.4

#define MS_THR   1.4
#define MS_QBITS 6//6//5//6  6 - for qantenne
#define MS_DBITS (MS_QBITS + 2) //16


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
} IMS_DEC_STATE;


typedef struct
{
	int nh; 
	int rh;
	int m;
    int n;
	int maxiter;
    int codec_id;
	int   **hd; 
	int   *syndr;                     

	
//#ifdef KEEP_STATISTIC
	double *sign_counter;
	double *min_abs_llr;
	double *prev_soft;
//#endif

	// Min-Sum Decoder
	MS_DATA *ms_soft;
	int	*ms_BnNS;
	MS_DEC_STATE *ms_dcs;
	MS_DEC_STATE *ms_tmps;
	MS_DATA *ms_buffer;
	MS_DATA *ms_rbuffer;
	MS_DATA *ms_rsoft;

	// Integer Min-Sum Decoder
	IMS_DATA *ims_soft;
	int	 *ims_BnNS;
	IMS_DATA *ims_y;
	IMS_DEC_STATE *ims_dcs;
	IMS_DEC_STATE *ims_tmps;
	IMS_DATA *ims_buffer;
	IMS_DATA *ims_rbuffer;
	IMS_DATA *ims_rsoft;

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
	double	*ilche_tmp;
	int  *ilche_soft_out;
	int  **ilche_state;

}DEC_STATE;


// functions prototypes
DEC_STATE* decod_open( int decoder_id, int q_bits, int mh, int nh, int M );
int decod_init( void* st );
void decod_close( DEC_STATE* st );
int min_sum_decod_qc_lm( DEC_STATE* st, int soft[], int decword[], int maxiter, double alpha );    
int imin_sum_decod_qc_lm( DEC_STATE* st, int soft[], int decword[], int maxiter, double alpha, double thr, int qbits, int dbits );    
int lmin_sum_decod_qc_lm( DEC_STATE* st, int soft[], int decword[], int maxiter, double alpha, double beta );    
int il_min_sum_decod_qc_lm( DEC_STATE* st, int soft[], int decword[], int maxiter, double alpha, double beta,  int inner_data_bits );    
int lche_decod( DEC_STATE* st, int soft[], int decword[], int maxiter );    
int ilche_decod( DEC_STATE* st, int soft[], int decword[], int maxsteps );


void icheck_syndrome( int **matr, int rh, int nh, IMS_DATA *soft, IMS_DATA *rsoft, int m_ldpc, int *synd );
void il_min_sum_init( IMS_DEC_STATE *prev, int r, int *signs, int n );




#endif	//_DECODERS_H_
