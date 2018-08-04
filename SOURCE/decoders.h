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
	LCHE_DEC
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
#define MS_QBITS 7//6//5//6  7 - for qantenne
#define MS_DBITS (MS_QBITS + 2) //16


#define UFLT_MNT_16
#define FLT_MNT_16
//#define FLT_POW_16

#ifdef _MSC_VER
typedef unsigned __int64    ui64;
typedef __int64             i64;
#else
typedef unsigned long long  ui64;
typedef long long           i64;
#endif

typedef unsigned int        ui32;
typedef unsigned short      ui16;
typedef short               i16;
typedef int                 i32;
typedef unsigned char       ui8;
typedef signed char         i8;


#if 01
typedef double MS_DATA;
//typedef unsigned short UDATA;
#else
typedef char MS_DATA;
typedef unsigned char UDATA;
#endif

typedef short IMS_DATA;


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
	int q_bits;
	int q;
	int nh; 
	int rh;
	int m;
    int n;
	int maxiter;
    int codec_id;
	int bin_codec;
	int max_rw;
	short   **hd; 
	short   *syndr;                     

	
//#ifdef KEEP_STATISTIC
	double *sign_counter;
	double *min_abs_llr;
	double *prev_soft;
//#endif

	// Min-Sum Decoder
	MS_DATA *ms_soft;
	short	*ms_BnNS;
	MS_DEC_STATE *ms_dcs;
	MS_DEC_STATE *ms_tmps;
	MS_DATA *ms_buffer;
	MS_DATA *ms_rbuffer;
	MS_DATA *ms_rsoft;

	// Integer Min-Sum Decoder
	IMS_DATA *ims_soft;
	short	 *ims_BnNS;
	IMS_DATA *ims_y;
	IMS_DEC_STATE *ims_dcs;
	IMS_DEC_STATE *ims_tmps;
	IMS_DATA *ims_buffer;
	IMS_DATA *ims_rbuffer;
	IMS_DATA *ims_rsoft;

	// Layered Min-Sum Decoder
	MS_DATA *lms_soft;
	short	*lms_BnNS;
	MS_DEC_STATE *lms_dcs;
	MS_DEC_STATE *lms_tmps;
	MS_DATA *lms_buffer;
	MS_DATA *lms_rbuffer;
	MS_DATA *lms_rsoft;

	// Low complexity-high efficienty Decoder
	double  *lche_data0;
	double	*lche_tmp;
	double  *lche_soft_out;
	double  **lche_state;

}DEC_STATE;


// functions prototypes
DEC_STATE* decod_open( int decoder_id, int q_bits, int mh, int nh, int M );
int decod_init( void* st );
void decod_close( DEC_STATE* st );
int min_sum_decod_qc_lm( DEC_STATE* st, int soft[], int decword[], int maxiter, double alpha );    
int imin_sum_decod_qc_lm( DEC_STATE* st, int soft[], int decword[], int maxiter, double alpha, double thr, int qbits, int dbits );    
int lmin_sum_decod_qc_lm( DEC_STATE* st, int soft[], int decword[], int maxiter, double alpha, double beta );    
int lche_decod( DEC_STATE* st, int soft[], int decword[], int maxiter );    




#endif	//_DECODERS_H_
