#include "stdafx.h"

#include <string.h>
#include <math.h>
#include <assert.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>


#include "decoders.h"

char const * const DEC_FULL_NAME[] = 
{
    "Belief Propagation",
    "Sum-Product",
    "Advanced Sum-Product",
    "Min-Sum",
    "Integer Min-Sum",
    "Integer Advanced Sum-Product",
    "FHT Sum-Product",
    "TDMP Advanced Sum-Product",
	"Layered Min-Sum",
	"Low complexity-high efficiency"
};



#define INP_FPP  16
#define ONE_INP (1 << INP_FPP)

#define PROB_FPP 16//16 //28  // less than 30 !!!
#define ONE_PROB (1 << PROB_FPP)

#define MAP_FPP  16//14//24
#define ONE_MAP (1 << MAP_FPP)

#define HAD_FPP 16//14
#define ONE_HAD (1 << HAD_FPP)

#define TMP_FPP 20
#define ONE_TMP (1 << TMP_FPP)

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#define QMAX  256
#define RWMAX 64

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


#define div_power2( x, n )    ((x) >> (n))
#define div_power2r( x, n )   (((x) + (1 << (n-1))) >> (n))

#define SOFT_FPP 12
#define P_FPP    32//32


#define MAX_UI16 0x7fff
#define ONE_SOFT (1 << SOFT_FPP)


#define MAX_SOFT (ONE_SOFT - 1) 
#define IASP_DEC_MAX_VAL (ONE_SOFT - 1)


#define INPUT_LIMIT 20.0

#define SP_DEC_MIN_VAL 0.000001
#define SP_DEC_MAX_VAL (1.0 - SP_DEC_MIN_VAL)


#define maxi( a, b )  (a) < (b) ? (b) : (a)
#define mini( a, b )  (b) < (a) ? (b) : (a)


static double mind(double a, double b) {if (a<b) return a; else return b;}
static double maxd(double a, double b) {if (a<b) return b; else return a;}
static double absd (double x)
{
	if (x<0) return -x; else return x;
}


void update_statistics( double *prev_soft, double *curr_soft, double *sign_counter, double *min_abs_llr, double thr, int n )
{
	int i;
	for( i = 0; i < n; i++ )
	{
		double val = curr_soft[i] - thr;
		double absval = val < 0.0 ? -val : val;
		min_abs_llr[i] = absval < min_abs_llr[i] ? absval : min_abs_llr[i];
		sign_counter[i] += ((prev_soft[i] - thr) * (curr_soft[i] - thr)) < 0;
		prev_soft[i] = curr_soft[i];
	}
}


double** Alloc2d_double( int b, int c )
{
	double **p;
	int i;

	p = (double**)calloc( b, sizeof(double*) );
	assert(p);
	p[0] = (double*)calloc( b*c, sizeof(double) );
	for( i = 1; i < b; i++ )
	{
		p[i] = p[i-1] + c;
		assert( p[i] );
	}
	return p;
}

short** Alloc2d_short( int b, int c )
{
	short **p;
	int i;

	p = (short**)calloc( b, sizeof(short*) );
	assert(p);
	p[0] = (short*)calloc( b*c, sizeof(short) );
	for( i = 1; i < b; i++ )
	{
		p[i] = p[i - 1] + c;
		assert( p[i] );
	}
	return p;
}

int** Alloc2d_int( int b, int c )
{
	int **p;
	int i;

	p = (int**)calloc( b, sizeof(int*) );
	assert(p);
	p[0] = (int*)calloc( b*c, sizeof(int) );
	for( i = 1; i < b; i++ )
	{
		p[i] = p[i-1] + c;
		assert( p[i] );
	}
	return p;
}


void free2d_int( int **p )
{
	free( p[0] );
	free( p );
}

void free2d_double( double **p )
{
	free( p[0] );
	free( p );
}

void free2d_short( short **p )
{
	free( p[0] );
	free( p );
}


void rotate_sign( short x[], short y[], int shift, int M )
{
	while( shift < 0 )
		shift += M;

	while( shift >= M )
		shift -= M;

	memcpy( y, x + shift, (M-shift)*sizeof(y[0]));
	memcpy( y+(M-shift), x, shift*sizeof(y[0]));
}

void rotate_data( double x[], double y[], int shift, int M )
{
	while( shift < 0 )
		shift += M;

	while( shift >= M )
		shift -= M;

	memcpy( y, x + shift, (M-shift)*sizeof(y[0]));
	memcpy( y+(M-shift), x, shift*sizeof(y[0]));
}

void rotate( void *x, void *y, int shift, int size, int M )
{
	int shift1;
	int shift2;
/*
	int shift1 = size * shift;
	int shift2 = size * (M-shift);
*/
	while( shift < 0 )
		shift += M;

	while( shift >= M )
		shift -= M;

	shift1 = size * shift;
	shift2 = size * (M-shift);

	memcpy( (char*)y, (char*)x + shift1, shift2 );
	memcpy( (char*)y + shift2, (char*)x, shift1 );
}

DEC_STATE* decod_open( int codec_id, int q_bits, int mh, int nh, int M )
{
    DEC_STATE* st;
    
	int N = nh * M;
	int R = mh * M;
	int BBsize = (mh > nh - mh) ? mh : nh - mh;


    st = (DEC_STATE*)calloc( 1, sizeof(DEC_STATE) );
    if( !st ) 
        return NULL;

	st->q_bits  = q_bits;
	st->nh      = nh;
	st->rh      = mh;
	st->m       = M;
    st->n       = N;
    st->codec_id = codec_id;

	switch( codec_id )
	{
	case MS_DEC:
	case IMS_DEC:
	case LMS_DEC:
	case LCHE_DEC:
		st->bin_codec = 1;
		break;
	}

	st->q = st->bin_codec ? 1 : 1 << q_bits;


	st->hd = Alloc2d_short( mh, nh );
	if( st->hd==NULL)
		return NULL;

#ifdef KEEP_STATISTIC
	st->sign_counter = (double*)calloc(N, sizeof(st->sign_counter[0]));
	st->min_abs_llr = (double*)calloc(N, sizeof(st->min_abs_llr[0]));
	st->prev_soft = (double*)calloc(N, sizeof(st->prev_soft[0]));
#endif


	st->syndr = (short*)calloc(R, sizeof(st->syndr[0]) );
	if( !st->syndr )
		return NULL;
    

	switch( codec_id )
	{
	case MS_DEC:
		st->ms_soft = (MS_DATA*)calloc(N, sizeof(st->ms_soft[0]) );
		if( !st->ms_soft )
			return NULL;

		st->ms_BnNS = (short*)calloc(mh*N, sizeof( st->ms_BnNS[0] ) );
		if( !st->ms_BnNS )
			return NULL;

		st->ms_dcs = (MS_DEC_STATE*)calloc(R, sizeof( st->ms_dcs[0] ) );
		if( !st->ms_dcs )
			return NULL;

		st->ms_tmps = (MS_DEC_STATE*)calloc(M, sizeof( st->ms_tmps[0] ) );
		if( !st->ms_tmps )
			return NULL;

		st->ms_buffer = (MS_DATA*)calloc(M, sizeof(st->ms_buffer[0]) );
		if( !st->ms_buffer )
			return NULL;

		st->ms_rbuffer = (MS_DATA*)calloc(M, sizeof(st->ms_rbuffer[0]) );
		if( !st->ms_rbuffer )
			return NULL;

		st->ms_rsoft = (MS_DATA*)calloc(M, sizeof(st->ms_rsoft[0]) );
		if( !st->ms_rsoft )
			return NULL;
		break;

	case IMS_DEC:
		st->ims_y = (IMS_DATA*)calloc(N, sizeof(st->ims_y[0]) );
		if( !st->ims_y )
			return NULL;

		st->ims_soft = (IMS_DATA*)calloc(N, sizeof(st->ims_soft[0]) );
		if( !st->ims_soft )
			return NULL;

		st->ims_BnNS = (short*)calloc(mh*N, sizeof( st->ims_BnNS[0] ) );
		if( !st->ims_BnNS )
			return NULL;

		st->ims_dcs = (IMS_DEC_STATE*)calloc(R, sizeof( st->ims_dcs[0] ) );
		if( !st->ims_dcs )
			return NULL;

		st->ims_tmps = (IMS_DEC_STATE*)calloc(M, sizeof( st->ims_tmps[0] ) );
		if( !st->ims_tmps )
			return NULL;

		st->ims_buffer = (IMS_DATA*)calloc(M, sizeof(st->ims_buffer[0]) );
		if( !st->ims_buffer )
			return NULL;

		st->ims_rbuffer = (IMS_DATA*)calloc(M, sizeof(st->ims_rbuffer[0]) );
		if( !st->ims_rbuffer )
			return NULL;

		st->ims_rsoft = (IMS_DATA*)calloc(M, sizeof(st->ims_rsoft[0]) );
		if( !st->ims_rsoft )
			return NULL;
		break;

	case LMS_DEC:
		st->lms_soft = (MS_DATA*)calloc(N, sizeof(st->lms_soft[0]) );
		if( !st->lms_soft )
			return NULL;

		st->lms_BnNS = (short*)calloc(mh*N, sizeof( st->lms_BnNS[0] ) );
		if( !st->lms_BnNS )
			return NULL;

		st->lms_dcs = (MS_DEC_STATE*)calloc(R, sizeof( st->lms_dcs[0] ) );
		if( !st->lms_dcs )
			return NULL;

		st->lms_tmps = (MS_DEC_STATE*)calloc(M, sizeof( st->lms_tmps[0] ) );
		if( !st->lms_tmps )
			return NULL;

		st->lms_buffer = (MS_DATA*)calloc(M, sizeof(st->lms_buffer[0]) );
		if( !st->lms_buffer )
			return NULL;

		st->lms_rbuffer = (MS_DATA*)calloc(M, sizeof(st->lms_rbuffer[0]) );
		if( !st->lms_rbuffer )
			return NULL;

		st->lms_rsoft = (MS_DATA*)calloc(M, sizeof(st->lms_rsoft[0]) );
		if( !st->lms_rsoft )
			return NULL;
		break;

	case LCHE_DEC:
		st->lche_data0 = (double*)calloc(N, sizeof( st->lche_data0[0] ) );
		if(!st->lche_data0 )
			return NULL;

		st->lche_tmp = (double*)calloc(nh, sizeof( st->lche_tmp[0] ) );
		if(!st->lche_tmp )
			return NULL;

		st->lche_soft_out = (double*)calloc(N, sizeof(st->lche_soft_out[0]) );
		if( !st->lche_soft_out )
			return NULL;

		st->lche_state = Alloc2d_double( mh*M, nh*M );
		if( st->lche_state==NULL)
			break;
		break;

	default: return NULL;
	}

    return st;
    
}

void check_syndrome( short **matr, int rh, int nh, MS_DATA *soft, MS_DATA *rsoft, int m_ldpc, short *synd )
{
	int k, j, n;

	for( j = 0; j < rh; j++ )
	{
		int stateOffset = j*m_ldpc;

		for( k = 0; k < nh; k++ )
		{
			int circ = matr[j][k];

			if( circ != SKIP )
			{
				rotate( &soft[k*m_ldpc], rsoft, circ, sizeof(soft[0]), m_ldpc ); 

				for( n = 0; n < m_ldpc; n++ )
					synd[stateOffset+n] ^= rsoft[n] < 0;
			}
		}
	}
}


void decod_close( DEC_STATE* st )
{
	int nh = st->nh;
	int mh = st->rh;
	int N = st->n;
	int M = st->m;
	int q_bits = st->q_bits;
	int q = st->q;
	int codec_id = st->codec_id;
	int R = mh * M;



	if(st->hd)		{ free2d_short( st->hd );          st->hd       = NULL; }
#ifdef KEEP_STATISTIC
	if(st->sign_counter)	{ free(st->sign_counter);                     st->sign_counter        = NULL; }
	if(st->min_abs_llr)		{ free(st->min_abs_llr);                     st->min_abs_llr        = NULL; }
	if(st->prev_soft)		{ free(st->prev_soft);                     st->prev_soft        = NULL; }
#endif



	if(st->syndr)	{ free( st->syndr ); st->syndr = NULL; }

	switch( codec_id )
	{
	case MS_DEC:
		if( st->ms_soft )    { free(st->ms_soft);    st->ms_soft    = NULL; }
		if( st->ms_BnNS )    { free(st->ms_BnNS);    st->ms_BnNS    = NULL; }
		if( st->ms_dcs )     { free(st->ms_dcs);     st->ms_dcs     = NULL; }
		if( st->ms_tmps )    { free(st->ms_tmps);    st->ms_tmps    = NULL; }
		if( st->ms_buffer )  { free(st->ms_buffer);  st->ms_buffer  = NULL; }
		if( st->ms_rbuffer ) { free(st->ms_rbuffer); st->ms_rbuffer = NULL; }
		if( st->ms_rsoft )   { free(st->ms_rsoft);   st->ms_rsoft   = NULL; }
		break;

	case IMS_DEC:
		if( st->ims_y )       { free(st->ims_y);       st->ims_y       = NULL; }
		if( st->ims_soft )    { free(st->ims_soft);    st->ims_soft    = NULL; }
		if( st->ims_BnNS )    { free(st->ims_BnNS);    st->ims_BnNS    = NULL; }
		if( st->ims_dcs )     { free(st->ims_dcs);     st->ims_dcs     = NULL; }
		if( st->ims_tmps )    { free(st->ims_tmps);    st->ims_tmps    = NULL; }
		if( st->ims_buffer )  { free(st->ims_buffer);  st->ims_buffer  = NULL; }
		if( st->ims_rbuffer ) { free(st->ims_rbuffer); st->ims_rbuffer = NULL; }
		if( st->ims_rsoft )   { free(st->ims_rsoft);   st->ims_rsoft   = NULL; }
		break;

	case LMS_DEC:
		if( st->lms_soft )    { free(st->lms_soft);    st->lms_soft    = NULL; }
		if( st->lms_BnNS )    { free(st->lms_BnNS);    st->lms_BnNS    = NULL; }
		if( st->lms_dcs )     { free(st->lms_dcs);     st->lms_dcs     = NULL; }
		if( st->lms_tmps )    { free(st->lms_tmps);    st->lms_tmps    = NULL; }
		if( st->lms_buffer )  { free(st->lms_buffer);  st->lms_buffer  = NULL; }
		if( st->lms_rbuffer ) { free(st->lms_rbuffer); st->lms_rbuffer = NULL; }
		if( st->lms_rsoft )   { free(st->lms_rsoft);   st->lms_rsoft   = NULL; }
		break;

	case LCHE_DEC:
		if( st->lche_data0 )     { free(st->lche_data0);            st->lche_data0    = NULL; }
		if( st->lche_tmp )     { free(st->lche_tmp);		        st->lche_tmp      = NULL; }
		if( st->lche_soft_out )  { free(st->lche_soft_out);         st->lche_soft_out = NULL; }
		if( st->lche_state )     { free2d_double( st->lche_state ); st->lche_state    = NULL; }
		break;
	default:;
	}

    free( st );
}



#define ROW_WEIGHT_MAX 1024


static inline MS_DATA logexp_int( MS_DATA x )
{
	static double A[] = 
	{	
		1.41e+00, 7.72e-01, 4.54e-01, 2.72e-01, 1.65e-01, 9.97e-02, 6.04e-02, 3.66e-02, 
		2.22e-02, 1.35e-02, 8.17e-03, 4.96e-03, 3.01e-03, 1.82e-03, 1.11e-03, 6.71e-04,
		4.07e-04, 2.47e-04, 1.50e-04, 9.08e-05, 5.51e-05, 3.34e-05, 2.03e-05, 1.23e-05, 
		7.45e-06, 4.52e-06, 2.74e-06, 1.66e-06, 1.01e-06, 6.12e-07, 3.71e-07, 2.25e-07
	};

	static double B[] = 
	{
		3.47, 2.77, 2.37, 2.08, 1.86, 1.69, 1.54, 1.41,
		1.29, 1.19, 1.11, 1.03, 0.95, 0.89, 0.83, 0.77,
		0.72, 0.67, 0.63, 0.59, 0.55, 0.52, 0.48, 0.45,
		0.43, 0.40, 0.37, 0.35, 0.33, 0.31, 0.29, 0.27
	};

	static double C[] =
	{
		6.93, 6.24, 5.83, 5.55, 5.32, 5.14, 4.99, 4.85,
		4.73, 4.63, 4.53, 4.45, 4.37, 4.29, 4.22, 4.16,
		4.10, 4.04, 3.99, 3.94, 3.89, 3.84, 3.80, 3.75,
		3.71, 3.67, 3.64, 3.60, 3.56, 3.53, 3.50, 3.47
	};


    if( x <= 0 ) x = 1.0/4096.0;
	if( x > 16.0 ) x = 16.0;

	if( x >= 2.0 )
		return  -A[(int)(2*x + 0.5) - 1]; 
	else
		if( x > 1.0/16.0 )
			return -B[(int)(16*x + 0.5) - 1];
		else
			if( x > 1.0/512.0 )
				return -C[(int)(512*x + 0.5) - 1];
			else  
			{
                double s = 0;
				while( x < 1.0/512.0 )
				{
					x *= 32;
					s -= 3.46;
				}

				return s-C[(int)(512*x + 0.5) - 1];
			}
}
void map_bin_llr(MS_DATA y[], int n )
{
	int i;
	static int hard[1024];
	static MS_DATA ay[1024];
	static MS_DATA alogpy[1024];
	static MS_DATA a[1024];
	static MS_DATA A[1024];
	int synd;
	MS_DATA sum;

	for( i = 0; i < n; i++ )
	{
		hard[i] = y[i] < 0;
	}	
	
	synd = 0;
	for( i = 0; i < n; i++ )
		synd ^=hard[i];

	for( i = 0; i < n; i++ )
		hard[i] ^= synd;

	for( i = 0; i < n; i++ )
		ay[i] = y[i] < 0.0 ? -y[i] : y[i];
	
	for( i = 0; i < n; i++ )
		alogpy[i] = logexp_int( ay[i] );

	sum = 0;
	for( i = 0; i < n; i++ )
		sum += alogpy[i];

	for( i = 0; i < n; i++ )
	{
		A[i] = alogpy[i] - sum;
	}

	for( i = 0; i < n; i++ )
		a[i] = logexp_int( A[i] );

	for( i = 0; i < n; i++ )
		y[i] = (2 * hard[i] - 1) * a[i];
}


int lche_decod( DEC_STATE* st, int soft[], int decword[], int maxsteps )
{ 
	int i, j, k;
	int steps;
	short **hd = st->hd;
	double **Z = st->lche_state;
	double *soft_out = st->lche_soft_out;
	short *syndr    = st->syndr;
	double *data0 = st->lche_data0;
	double *u = st->lche_data0;
	double *y = st->lche_tmp;
	int synd;

	int rh = st->rh;
	int nh = st->nh;
	int m  = st->m;
	int r = rh * m;
	int n = nh * m;


	for( i = 0; i < r; i++ )
		for(j = 0; j < n; j++ )
			Z[i][j] = 0.0;

	// just to compute syndrome before iterations
	for( i = 0; i < n; i++ )
		soft_out[i] = soft[i];

#ifdef KEEP_STATISTIC
	for( i = 0; i < n; i++ ) st->prev_soft[i] = soft_out[i];
	for( i = 0; i < n; i++ ) st->sign_counter[i] = 0;
	for( i = 0; i < n; i++ ) st->min_abs_llr[i]  = 10000;
#endif

	memset( syndr, 0, r * sizeof(syndr[0]) );
	check_syndrome( hd, rh, nh, soft_out, data0, m, syndr );

	synd = 0;
	for( i = 0; i < r; i++ )	synd |= syndr[i];
	if( synd == 0 )
	{
		for( i = 0; i < n; i++ )
			decword[i] = soft_out[i] < 0;

		return 1; // 0; 
	}

	steps = 0; // number of iterations
	while( steps < maxsteps )
	{
		//	START ITERATIONS

		for( i = 0; i < rh; i++ )	//loop over checks
		{
			int cnt;

			for( k = 0; k < m; k++ )
			{
				double *a = Z[i * m + k];

				cnt = 0;
				for( j = 0; j < nh; j++ )
				{
					int circ = hd[i][j];

					if( circ != -1 )
					{
						int idx = j * m + ((k + circ) % m);

						u[cnt] = y[cnt] = soft_out[idx] - a[idx];

						cnt++;
					}
				}

				map_bin_llr( u, cnt );

				cnt = 0;
				for( j = 0; j < nh; j++ )
				{
					int circ = hd[i][j];

					if( circ != -1 )
					{
						int idx = j * m + ((k + circ) % m);

						soft_out[idx] = u[cnt] + y[cnt];
						a[idx] = u[cnt];
						cnt++;
					}
				}
			}
		}

#ifdef KEEP_STATISTIC
		update_statistics( st->prev_soft, soft_out, st->sign_counter, st->min_abs_llr, 0.5, n );
#endif
		//		synd = check_syndrome_thr( syndr, r, hd, rh, nh, m, soft_out, data0, 0.5 );

		steps = steps+1;

		memset( syndr, 0, r * sizeof(syndr[0]) );
		check_syndrome( hd, rh, nh, soft_out, data0, m, syndr );
		synd = 0;
		for( i = 0; i < r; i++ )	synd |= syndr[i];

		if( synd == 0 ) 
			break;
	}

	for( i = 0; i < n; i++ )
		decword[i] = soft_out[i] < 0.0;

	if( synd == 1 )
		steps = -steps;  // errors detected but not corrected

	return steps;
}







#if 1
#define BPS	15 // bits per soft
#define BPV 15 // bits per inner variables
#define MAX_VAL ((1L << BPV) - 1)
#else
#define BPS	3 // bits per soft
#define BPV 5 // bits per inner variables
#define MAX_VAL ((1 << BPV) - 1)
#endif

#define limit_val( x, max_val )	(x) > max_val ? max_val : ((x) < -max_val ? -max_val : (x))




#ifndef MS_MUL_CORRECTION
int min_sum_decod_qc_lm( DEC_STATE* st, int y[], int decword[], int maxsteps, double alpha )
{
	int   i, j, k, n;
	int iter;
	int parity;
	short **matr       = st->hd;
	short *synd        = st->syndr; 
	MS_DATA *soft      = st->ms_soft;
	short *BnNS        = st->ms_BnNS;  
	MS_DATA *buffer    = st->ms_buffer;
	MS_DATA *rbuffer   = st->ms_rbuffer;
	MS_DATA *rsoft     = st->ms_rsoft;
	MS_DEC_STATE *dcs  = st->ms_dcs;		
	MS_DEC_STATE *tmps = st->ms_tmps;		


	int m  = st->m;
	int rh = st->rh;
	int nh = st->nh;
	
	int c_ldpc = m;
	int r_ldpc = m * rh;
	int n_ldpc = m * nh;


	for( i = 0; i < r_ldpc; i++ )
	{
		dcs[i].min1 = 0;
		dcs[i].min2 = 0;
		dcs[i].pos  = 0;
		dcs[i].sign = 0;
	}
/*	
	for( i = 0; i < n_ldpc; i++ )
	{
		double val = abs(y[i]) * 4;
		int sign = y[i] < 0;
		val = (int)(val + 0.5);
		y[i] = sign ? -val : val;
	}
*/

    {
        double coef;
        double en = 0.0;
  
        for( i = 0; i < n_ldpc; i++ )
            en += y[i] * y[i];
        
        coef = sqrt( n_ldpc / en );
        for( i = 0; i < n_ldpc; i++ )
            y[i] = y[i] * coef;
    }   
    
    memset( BnNS, 0, rh*n_ldpc*sizeof(BnNS[0]) );




	// check input codeword
	for( j = 0; j < rh; j++ )
	{
		for( i = 0; i < nh; i++ ) 
		{
			int pos_r = j * m;
			int pos_n = i * m;
			int circ = matr[j][i];

			if( circ != - 1 )
			{
				rotate( &y[pos_n], rsoft, circ, sizeof(y[0]), m );

				for( k = 0; k < m; k++ )
					synd[pos_r + k] ^= rsoft[k] < 0;
			}
		}
	}

	parity = 0;
	for( i = 0; i < r_ldpc; i++ )
		parity |= synd[i];


	for( iter = 0; iter < maxsteps; iter++ )
	{

		int memOffset;


		// INIT_STAGE
		memset( synd, 0, r_ldpc*sizeof(synd[0]) );
		memset( soft, 0, n_ldpc*sizeof(soft[0]) );
		memOffset = 0;

		// STATE 1: compute sum
		for( j = 0; j < rh; j++ ) 
		{
			int stateOffset = j*c_ldpc;

			for( k = 0; k < nh; k++ )
			{
				int circ = matr[j][k];

				if( circ != SKIP )
				{
					for( n = 0; n < c_ldpc; n++ )
					{
						MS_DATA tmp = dcs[stateOffset+n].pos == k ? dcs[stateOffset+n].min2 : dcs[stateOffset+n].min1;

                        tmp -= alpha;
                        if( tmp < 0 )
                            tmp = 0;
                        
						buffer[n] = (BnNS[memOffset+n] ^ dcs[stateOffset+n].sign) ? -tmp : tmp;
					}

					rotate( buffer, rbuffer, c_ldpc - circ, sizeof(buffer[0]), c_ldpc ); 

					for( n = 0; n < c_ldpc; n++ )
					{
						MS_DATA t = soft[k * c_ldpc + n] + rbuffer[n];
						soft[k * c_ldpc + n] = t;//limit_val( t, MAX_VAL );
					}

				}

				memOffset += c_ldpc;
			}

		} 

		// STATE 2: 
		for( k = 0; k < n_ldpc; k++ )
		{
			soft[k] = y[k] + soft[k];
			decword[k] = soft[k] < 0;
		}


		// STATE 3: update statistic
		memOffset = 0;
		for( j = 0; j < rh; j++ )
		{
			int stateOffset = j*c_ldpc;

			for( k = 0; k < c_ldpc; k++ )
			{
				tmps[k].min1 = MAX_VAL;
				tmps[k].min2 = MAX_VAL;
				tmps[k].pos  = 0;
				tmps[k].sign = 0;
			}

			for( k = 0; k < nh; k++ )
			{
				int circ = matr[j][k];

				if( circ != SKIP )
				{
					rotate( &soft[k*c_ldpc], rsoft, circ, sizeof(soft[0]), c_ldpc ); 

					for( n = 0; n < c_ldpc; n++ )
						synd[stateOffset+n] ^= rsoft[n] < 0;

					for( n = 0; n < c_ldpc; n++ )
                    {
                        MS_DATA tmp = dcs[stateOffset+n].pos == k ? dcs[stateOffset+n].min2 : dcs[stateOffset+n].min1;
                        
                        tmp -= alpha;
                        if( tmp < 0 )
                            tmp = 0;
                        
						buffer[n] = tmp;
                    }

					for( n = 0; n < c_ldpc; n++ )
					{
						short sign = BnNS[memOffset+n] ^ dcs[stateOffset+n].sign;
						MS_DATA val  = buffer[n];// * alpha;
						MS_DATA t    = (BnNS[memOffset+n] ^ dcs[stateOffset+n].sign) ? -val : val;

						t  = rsoft[n] - t;

						sign = t < 0;
						BnNS[memOffset+n] = sign;
						tmps[n].sign ^= sign;

//						val = abs( t );	// incorrect for double t
						val = t < 0.0 ? -t : t;
						val = (val > MAX_VAL) ? MAX_VAL : val;


						if( val < tmps[n].min1 )      
						{
							tmps[n].pos = k;
							tmps[n].min2 = tmps[n].min1; 
							tmps[n].min1 = val;
						}
						else
						{
							if( val < tmps[n].min2 )
								tmps[n].min2 = val;
						}

					}  
				}


				memOffset += c_ldpc;
			}

			for( k = 0; k < c_ldpc; k++ )
				dcs[stateOffset+k] = tmps[k];

		}

		for( parity = 0, i = 0; i < r_ldpc; i++ )
			parity |= synd[i];


		if( parity == 0 )
			break; 

	} // stage

	return parity ? -iter : iter+1; 
}

#else
int min_sum_decod_qc_lm( DEC_STATE* st, int y[], int decword[], int maxsteps, double alpha )
{
	int   i, j, k, n;
	int iter;
	int parity;
	short **matr       = st->hd;
	short *synd        = st->syndr; 
	MS_DATA *soft      = st->ms_soft;
	short *BnNS        = st->ms_BnNS;  
	MS_DATA *buffer    = st->ms_buffer;
	MS_DATA *rbuffer   = st->ms_rbuffer;
	MS_DATA *rsoft     = st->ms_rsoft;
	MS_DEC_STATE *dcs  = st->ms_dcs;		
	MS_DEC_STATE *tmps = st->ms_tmps;		


	int m  = st->m;
	int rh = st->rh;
	int nh = st->nh;
	
	int c_ldpc = m;
	int r_ldpc = m * rh;
	int n_ldpc = m * nh;


	for( i = 0; i < r_ldpc; i++ )
	{
		dcs[i].min1 = 0;
		dcs[i].min2 = 0;
		dcs[i].pos  = 0;
		dcs[i].sign = 0;
	}
/*	
	for( i = 0; i < n_ldpc; i++ )
	{
		double val = abs(y[i]) * 4;
		int sign = y[i] < 0;
		val = (int)(val + 0.5);
		y[i] = sign ? -val : val;
	}
*/

	memset( BnNS, 0, rh*n_ldpc*sizeof(BnNS[0]) );
	
	// check input codeword
	for( j = 0; j < rh; j++ )
	{
		for( i = 0; i < nh; i++ ) 
		{
			int pos_r = j * m;
			int pos_n = i * m;
			int circ = matr[j][i];

			if( circ != - 1 )
			{
				rotate( &y[pos_n], rsoft, circ, sizeof(y[0]), m );

				for( k = 0; k < m; k++ )
					synd[pos_r + k] ^= rsoft[k] < 0;
			}
		}
	}

	parity = 0;
	for( i = 0; i < r_ldpc; i++ )
		parity |= synd[i];


	for( iter = 0; iter < maxsteps; iter++ )
	{

		int memOffset;


		// INIT_STAGE
		memset( synd, 0, r_ldpc*sizeof(synd[0]) );
		memset( soft, 0, n_ldpc*sizeof(soft[0]) );
		memOffset = 0;

		// STATE 1: compute sum
		for( j = 0; j < rh; j++ ) 
		{
			int stateOffset = j*c_ldpc;

			for( k = 0; k < nh; k++ )
			{
				int circ = matr[j][k];

				if( circ != SKIP )
				{
					for( n = 0; n < c_ldpc; n++ )
					{
						MS_DATA tmp = dcs[stateOffset+n].pos == k ? dcs[stateOffset+n].min2 : dcs[stateOffset+n].min1;

						buffer[n] = (BnNS[memOffset+n] ^ dcs[stateOffset+n].sign) ? -tmp : tmp;
					}

					rotate( buffer, rbuffer, c_ldpc - circ, sizeof(buffer[0]), c_ldpc ); 

					for( n = 0; n < c_ldpc; n++ )
					{
						MS_DATA t = soft[k * c_ldpc + n] + rbuffer[n];
						soft[k * c_ldpc + n] = t;//limit_val( t, MAX_VAL );
					}

				}

				memOffset += c_ldpc;
			}

		} 

		// STATE 2: 
		for( k = 0; k < n_ldpc; k++ )
		{
			soft[k] = y[k] + soft[k] * alpha; //alpha-normalization
			decword[k] = soft[k] < 0;
		}


		// STATE 3: update statistic
		memOffset = 0;
		for( j = 0; j < rh; j++ )
		{
			int stateOffset = j*c_ldpc;

			for( k = 0; k < c_ldpc; k++ )
			{
				tmps[k].min1 = MAX_VAL;
				tmps[k].min2 = MAX_VAL;
				tmps[k].pos  = 0;
				tmps[k].sign = 0;
			}

			for( k = 0; k < nh; k++ )
			{
				int circ = matr[j][k];

				if( circ != SKIP )
				{
					rotate( &soft[k*c_ldpc], rsoft, circ, sizeof(soft[0]), c_ldpc ); 

					for( n = 0; n < c_ldpc; n++ )
						synd[stateOffset+n] ^= rsoft[n] < 0;

					for( n = 0; n < c_ldpc; n++ )
						buffer[n] = dcs[stateOffset+n].pos == k ? dcs[stateOffset+n].min2 : dcs[stateOffset+n].min1;

					for( n = 0; n < c_ldpc; n++ )
					{
						short sign = BnNS[memOffset+n] ^ dcs[stateOffset+n].sign;
						MS_DATA val  = buffer[n] * alpha;
						MS_DATA t    = (BnNS[memOffset+n] ^ dcs[stateOffset+n].sign) ? -val : val;

						t  = rsoft[n] - t;

						sign = t < 0;
						BnNS[memOffset+n] = sign;
						tmps[n].sign ^= sign;

//						val = abs( t );	// incorrect for double t
						val = t < 0.0 ? -t : t;
						val = (val > MAX_VAL) ? MAX_VAL : val;


						if( val < tmps[n].min1 )      
						{
							tmps[n].pos = k;
							tmps[n].min2 = tmps[n].min1; 
							tmps[n].min1 = val;
						}
						else
						{
							if( val < tmps[n].min2 )
								tmps[n].min2 = val;
						}

					}  
				}


				memOffset += c_ldpc;
			}

			for( k = 0; k < c_ldpc; k++ )
				dcs[stateOffset+k] = tmps[k];

		}

		for( parity = 0, i = 0; i < r_ldpc; i++ )
			parity |= synd[i];


		if( parity == 0 )
			break; 

	} // stage

	return parity ? -iter : iter+1; 
}
#endif

#ifndef MS_MUL_CORRECTION
int min_sum_decod_qc_lm( DEC_STATE* st, int y[], int decword[], int maxsteps, double alpha )
{
	int   i, j, k, n;
	int iter;
	int parity;
	short **matr       = st->hd;
	short *synd        = st->syndr; 
	MS_DATA *soft      = st->ms_soft;
	short *BnNS        = st->ms_BnNS;  
	MS_DATA *buffer    = st->ms_buffer;
	MS_DATA *rbuffer   = st->ms_rbuffer;
	MS_DATA *rsoft     = st->ms_rsoft;
	MS_DEC_STATE *dcs  = st->ms_dcs;		
	MS_DEC_STATE *tmps = st->ms_tmps;		


	int m  = st->m;
	int rh = st->rh;
	int nh = st->nh;
	
	int c_ldpc = m;
	int r_ldpc = m * rh;
	int n_ldpc = m * nh;


	for( i = 0; i < r_ldpc; i++ )
	{
		dcs[i].min1 = 0;
		dcs[i].min2 = 0;
		dcs[i].pos  = 0;
		dcs[i].sign = 0;
	}
/*	
	for( i = 0; i < n_ldpc; i++ )
	{
		double val = abs(y[i]) * 4;
		int sign = y[i] < 0;
		val = (int)(val + 0.5);
		y[i] = sign ? -val : val;
	}
*/

    {
        double coef;
        double en = 0.0;
  
        for( i = 0; i < n_ldpc; i++ )
            en += y[i] * y[i];
        
        coef = sqrt( n_ldpc / en );
        for( i = 0; i < n_ldpc; i++ )
            y[i] = y[i] * coef;
    }   
    
    memset( BnNS, 0, rh*n_ldpc*sizeof(BnNS[0]) );




	// check input codeword
	for( j = 0; j < rh; j++ )
	{
		for( i = 0; i < nh; i++ ) 
		{
			int pos_r = j * m;
			int pos_n = i * m;
			int circ = matr[j][i];

			if( circ != - 1 )
			{
				rotate( &y[pos_n], rsoft, circ, sizeof(y[0]), m );

				for( k = 0; k < m; k++ )
					synd[pos_r + k] ^= rsoft[k] < 0;
			}
		}
	}

	parity = 0;
	for( i = 0; i < r_ldpc; i++ )
		parity |= synd[i];


	for( iter = 0; iter < maxsteps; iter++ )
	{

		int memOffset;


		// INIT_STAGE
		memset( synd, 0, r_ldpc*sizeof(synd[0]) );
		memset( soft, 0, n_ldpc*sizeof(soft[0]) );
		memOffset = 0;

		// STATE 1: compute sum
		for( j = 0; j < rh; j++ ) 
		{
			int stateOffset = j*c_ldpc;

			for( k = 0; k < nh; k++ )
			{
				int circ = matr[j][k];

				if( circ != SKIP )
				{
					for( n = 0; n < c_ldpc; n++ )
					{
						MS_DATA tmp = dcs[stateOffset+n].pos == k ? dcs[stateOffset+n].min2 : dcs[stateOffset+n].min1;

                        tmp -= alpha;
                        if( tmp < 0 )
                            tmp = 0;
                        
						buffer[n] = (BnNS[memOffset+n] ^ dcs[stateOffset+n].sign) ? -tmp : tmp;
					}

					rotate( buffer, rbuffer, c_ldpc - circ, sizeof(buffer[0]), c_ldpc ); 

					for( n = 0; n < c_ldpc; n++ )
					{
						MS_DATA t = soft[k * c_ldpc + n] + rbuffer[n];
						soft[k * c_ldpc + n] = t;//limit_val( t, MAX_VAL );
					}

				}

				memOffset += c_ldpc;
			}

		} 

		// STATE 2: 
		for( k = 0; k < n_ldpc; k++ )
		{
			soft[k] = y[k] + soft[k];
			decword[k] = soft[k] < 0;
		}


		// STATE 3: update statistic
		memOffset = 0;
		for( j = 0; j < rh; j++ )
		{
			int stateOffset = j*c_ldpc;

			for( k = 0; k < c_ldpc; k++ )
			{
				tmps[k].min1 = MAX_VAL;
				tmps[k].min2 = MAX_VAL;
				tmps[k].pos  = 0;
				tmps[k].sign = 0;
			}

			for( k = 0; k < nh; k++ )
			{
				int circ = matr[j][k];

				if( circ != SKIP )
				{
					rotate( &soft[k*c_ldpc], rsoft, circ, sizeof(soft[0]), c_ldpc ); 

					for( n = 0; n < c_ldpc; n++ )
						synd[stateOffset+n] ^= rsoft[n] < 0;

					for( n = 0; n < c_ldpc; n++ )
                    {
                        MS_DATA tmp = dcs[stateOffset+n].pos == k ? dcs[stateOffset+n].min2 : dcs[stateOffset+n].min1;
                        
                        tmp -= alpha;
                        if( tmp < 0 )
                            tmp = 0;
                        
						buffer[n] = tmp;
                    }

					for( n = 0; n < c_ldpc; n++ )
					{
						short sign = BnNS[memOffset+n] ^ dcs[stateOffset+n].sign;
						MS_DATA val  = buffer[n];// * alpha;
						MS_DATA t    = (BnNS[memOffset+n] ^ dcs[stateOffset+n].sign) ? -val : val;

						t  = rsoft[n] - t;

						sign = t < 0;
						BnNS[memOffset+n] = sign;
						tmps[n].sign ^= sign;

//						val = abs( t );	// incorrect for double t
						val = t < 0.0 ? -t : t;
						val = (val > MAX_VAL) ? MAX_VAL : val;


						if( val < tmps[n].min1 )      
						{
							tmps[n].pos = k;
							tmps[n].min2 = tmps[n].min1; 
							tmps[n].min1 = val;
						}
						else
						{
							if( val < tmps[n].min2 )
								tmps[n].min2 = val;
						}

					}  
				}


				memOffset += c_ldpc;
			}

			for( k = 0; k < c_ldpc; k++ )
				dcs[stateOffset+k] = tmps[k];

		}

		for( parity = 0, i = 0; i < r_ldpc; i++ )
			parity |= synd[i];


		if( parity == 0 )
			break; 

	} // stage

	return parity ? -iter : iter+1; 
}

#else

void process_check_node( MS_DEC_STATE *curr, short curr_v2c_sign, MS_DATA abs_curr_v2c, int index )
{
	curr->sign ^= curr_v2c_sign;

	if( abs_curr_v2c < curr->min1 )      
	{
		curr->pos = index;
		curr->min2 = curr->min1; 
		curr->min1 = abs_curr_v2c;
	}
	else
	{
		if( abs_curr_v2c < curr->min2 )
			curr->min2 = abs_curr_v2c;
	}
}

#define MY_VERSION
//#define ALPHA 0.8
//#define BETA  0.35 

#define ALPHA 1.0
#define BETA  0.4

MS_DATA	get_c2v_val( short sign, int pos, MS_DEC_STATE *stat, double alpha, double beta )
{
	MS_DATA c2v_abs  = stat->pos == pos ? stat->min2 : stat->min1;

	c2v_abs  = c2v_abs * alpha - beta; 
	if( c2v_abs < 0 )
		c2v_abs = 0;

	return sign ? -c2v_abs : c2v_abs;
}

double beta_control( MS_DEC_STATE *state )
{
	double beta = 0.4;
#if 0

//	if( state->min2 > state->min1 * 1.1 )	beta = 0.4;
//	if( state->min2 > state->min1 * 1.2 )	beta = 0.7;


	//if( state->min2 - state->min1 > 0.1 )	beta = 0.6;
	if( state->min2 - state->min1 > 0.2 )	beta = 0.70;
	if( state->min2 - state->min1 > 0.5 )	beta = 0.75;
	
#endif
	return beta;
}

int lmin_sum_decod_qc_lm( DEC_STATE* st, int y[], int decword[], int maxsteps, double alpha, double beta )
{
	int   i, j, k, n;
	int iter;
	int parity;
	short **matr       = st->hd;
	short *synd        = st->syndr; 
	MS_DATA *soft      = st->lms_soft;
	short *signs       = st->lms_BnNS;  
	MS_DATA *buffer    = st->lms_buffer;
	MS_DATA *rbuffer   = st->lms_rbuffer;
	MS_DATA *rsoft     = st->lms_rsoft;
	MS_DEC_STATE *prev = st->lms_dcs;		
	MS_DEC_STATE *curr = st->lms_tmps;		


	int m  = st->m;
	int rh = st->rh;
	int nh = st->nh;
	
	int m_ldpc = m;
	int r_ldpc = m * rh;
	int n_ldpc = m * nh;

	for( i = 0; i < n_ldpc; i++ )
		soft[i] = y[i];

	for( i = 0; i < r_ldpc; i++ )
	{
		prev[i].min1 = 0;
		prev[i].min2 = 0;
		prev[i].pos  = 0;
		prev[i].sign = 0;
	}


#ifdef KEEP_STATISTIC
	for( i = 0; i < n_ldpc; i++ ) st->prev_soft[i] = soft[i];
	for( i = 0; i < n_ldpc; i++ ) st->sign_counter[i] = 0;
	for( i = 0; i < n_ldpc; i++ ) st->min_abs_llr[i]  = 10000;
#endif

#if 01

	memset( signs, 0, rh*n_ldpc*sizeof(signs[0]) );

	// check input codeword
	memset( synd, 0, r_ldpc*sizeof(synd[0]) );
	check_syndrome( matr, rh, nh, soft, rsoft, m_ldpc, synd );  

	parity = 0;
	for( i = 0; i < r_ldpc; i++ )	parity |= synd[i];

	for( iter = 0; iter < maxsteps; iter++ )
	{
		if( parity == 0 )
			break; 

		// INIT_STAGE
		memset( synd, 0, r_ldpc*sizeof(synd[0]) );

		// update statistic
		for( j = 0; j < rh; j++ )
		{
			MS_DATA *prev_c2v_abs = buffer;
			int stateOffset = j*m_ldpc;

			for( k = 0; k < m_ldpc; k++ )
			{
				curr[k].min1 = MAX_VAL;
				curr[k].min2 = MAX_VAL;
				curr[k].pos  = 0;
				curr[k].sign = 0;
			}


#ifdef MY_VERSION
			for( k = 0; k < nh; k++ )
			{
				int memOffset = j * n_ldpc + k * m_ldpc;
				MS_DATA *curr_soft    = &soft[k * m_ldpc];
				int circ = matr[j][k];

				if( circ != SKIP )
				{
					rotate( &soft[k*m_ldpc], rsoft, circ, sizeof(soft[0]), m_ldpc ); 

					for( n = 0; n < m_ldpc; n++ )
						prev_c2v_abs[n] = prev[stateOffset+n].pos == k ? prev[stateOffset+n].min2 : prev[stateOffset+n].min1;

					for( n = 0; n < m_ldpc; n++ )
					{
						short	prev_c2v_sgn  = signs[memOffset+n] ^ prev[stateOffset+n].sign;
						MS_DATA	prev_c2v_val  = prev_c2v_sgn ? -prev_c2v_abs[n] : prev_c2v_abs[n];
						MS_DATA curr_v2c_val  = rsoft[n] - prev_c2v_val;
#if 0
						short	curr_v2c_sgn = curr_v2c_val < 0;
						MS_DATA curr_v2c_abs = (curr_v2c_val < 0.0 ? -curr_v2c_val : curr_v2c_val) *alpha; //scaling 
#else
						const MS_DATA beta = 0.4;// 0.4;
						short	curr_v2c_sgn = curr_v2c_val < 0;
						MS_DATA curr_v2c_abs = (curr_v2c_val < 0.0 ? -curr_v2c_val : curr_v2c_val); //shifting 
						curr_v2c_abs -= beta;
						//curr_v2c_sgn = curr_v2c_abs < 0 ? 0 : curr_v2c_sgn;
						curr_v2c_abs = curr_v2c_abs < 0 ? 0 : curr_v2c_abs;
#endif
						/*curr_v2c*/curr_soft[n] = curr_v2c_val;
						signs[memOffset+n] = curr_v2c_sgn;

						process_check_node( &curr[n], curr_v2c_sgn, curr_v2c_abs, k );

					}  
				}
			}

			for( k = 0; k < m_ldpc; k++ )
				prev[stateOffset+k] = curr[k];

			for( k = 0; k < nh; k++ )
			{
				MS_DATA *curr_c2v_val = buffer;
				MS_DATA *curr_soft    = &soft[k * m_ldpc];
				int memOffset = j * n_ldpc + k * m_ldpc;
				int circ = matr[j][k];

				if( circ != SKIP )
				{
					for( n = 0; n < m_ldpc; n++ )
					{
						MS_DATA curr_c2v_abs = prev[stateOffset+n].pos == k ? prev[stateOffset+n].min2 : prev[stateOffset+n].min1;

						curr_c2v_val[n] = (signs[memOffset+n] ^ prev[stateOffset+n].sign) ? -curr_c2v_abs : curr_c2v_abs;
					}

					for( n = 0; n < m_ldpc; n++ )
						buffer[n] = /*curr_v2c*/curr_soft[n] + curr_c2v_val[n];

					rotate( buffer, rbuffer, m_ldpc - circ, sizeof(curr_c2v_val[0]), m_ldpc ); 

					for( n = 0; n < m_ldpc; n++ )
						curr_soft[n] = rbuffer[n];

				}
			}
#else
			for( k = 0; k < nh; k++ )
			{
				int memOffset = j * n_ldpc + k * m_ldpc;
				MS_DATA *curr_soft    = &soft[k * m_ldpc];
				int circ = matr[j][k];

				if( circ != SKIP )
				{
					rotate( &soft[k*m_ldpc], rsoft, circ, sizeof(soft[0]), m_ldpc ); 

					for( n = 0; n < m_ldpc; n++ )
					{
						short	prev_c2v_sgn  = signs[memOffset+n] ^ prev[stateOffset+n].sign;
#if 0
						MS_DATA	prev_c2v_val  = get_c2v_val( prev_c2v_sgn, k, &prev[stateOffset+n], 1.0/*alpha*/, beta );
#else
						MS_DATA beta          = beta_control( &prev[stateOffset+n] );
						MS_DATA	prev_c2v_val  = get_c2v_val( prev_c2v_sgn, k, &prev[stateOffset+n], 1.0/*alpha*/, beta );
#endif					
						MS_DATA curr_v2c_val  = rsoft[n] - prev_c2v_val;
						short	curr_v2c_sgn = curr_v2c_val < 0;
						MS_DATA curr_v2c_abs = (curr_v2c_val < 0.0 ? -curr_v2c_val : curr_v2c_val); 

						curr_soft[n]       = curr_v2c_val;
						signs[memOffset+n] = curr_v2c_sgn;

						process_check_node( &curr[n], curr_v2c_sgn, curr_v2c_abs, k );
					}  
				}
			}

			for( k = 0; k < m_ldpc; k++ )
				prev[stateOffset+k] = curr[k];

			for( k = 0; k < nh; k++ )
			{
				MS_DATA *curr_c2v_val = buffer;
				MS_DATA *curr_soft    = &soft[k * m_ldpc];
				int memOffset = j * n_ldpc + k * m_ldpc;
				int circ = matr[j][k];

				if( circ != SKIP )
				{
					for( n = 0; n < m_ldpc; n++ )
					{
						short curr_c2v_sgn  = signs[memOffset+n] ^ prev[stateOffset+n].sign;
#if 0
						curr_c2v_val[n]  = get_c2v_val( curr_c2v_sgn, k, &prev[stateOffset+n], 1.0/*alpha*/, beta );
#else
						MS_DATA beta     = beta_control( &prev[stateOffset+n] );
						curr_c2v_val[n]  = get_c2v_val( curr_c2v_sgn, k, &prev[stateOffset+n], 1.0/*alpha*/, beta );
#endif					
					}

					for( n = 0; n < m_ldpc; n++ )
						buffer[n] = curr_soft[n] + curr_c2v_val[n];

					rotate( buffer, rbuffer, m_ldpc - circ, sizeof(curr_c2v_val[0]), m_ldpc ); 

					for( n = 0; n < m_ldpc; n++ )
						curr_soft[n] = rbuffer[n];

				}
			}
#endif
		} 

#ifdef KEEP_STATISTIC
		update_statistics( st->prev_soft, soft, st->sign_counter, st->min_abs_llr, 0.0, n_ldpc );
#endif

		// check syndrome
		check_syndrome( matr, rh, nh, soft, rsoft, m_ldpc, synd );  

		for( parity = 0, i = 0; i < r_ldpc; i++ )
			parity |= synd[i];


		if( parity == 0 )
			break; 

	} // stage
#else
	memset( signs, 0, rh*n_ldpc*sizeof(signs[0]) );

	// check input codeword
	memset( synd, 0, r_ldpc*sizeof(synd[0]) );
	check_syndrome( matr, rh, nh, soft, rsoft, m_ldpc, synd );  

	parity = 0;
	for( i = 0; i < r_ldpc; i++ )
		parity |= synd[i];

	for( iter = 0; iter < maxsteps; iter++ )
	{
		if( parity == 0 )
			break; 

		// INIT_STAGE
		memset( synd, 0, r_ldpc*sizeof(synd[0]) );
	
		// update statistic
		for( j = 0; j < rh; j++ )
		{
			int stateOffset = j*m_ldpc;

			for( k = 0; k < m_ldpc; k++ )
			{
				curr[k].min1 = MAX_VAL;
				curr[k].min2 = MAX_VAL;
				curr[k].pos  = 0;
				curr[k].sign = 0;
			}

			for( k = 0; k < nh; k++ )
			{
				int memOffset = j * n_ldpc + k * m_ldpc;
				int circ = matr[j][k];

				if( circ != SKIP )
				{
					rotate( &soft[k*m_ldpc], rsoft, circ, sizeof(soft[0]), m_ldpc ); 

					for( n = 0; n < m_ldpc; n++ )
						buffer[n] = prev[stateOffset+n].pos == k ? prev[stateOffset+n].min2 : prev[stateOffset+n].min1;

					for( n = 0; n < m_ldpc; n++ )
					{
						short	prev_c2v_sgn = signs[memOffset+n] ^ prev[stateOffset+n].sign;
						MS_DATA	prev_c2v_val = prev_c2v_sgn ? -buffer[n] : buffer[n];
						MS_DATA curr_v2c_val = rsoft[n] - prev_c2v_val;
#if 0
						short	curr_v2c_sgn = curr_v2c_val < 0;
						MS_DATA curr_v2c_abs = (curr_v2c_val < 0.0 ? -curr_v2c_val : curr_v2c_val) * alpha; //scaling 
#else
						const MS_DATA beta = 0.4;
						short	curr_v2c_sgn = curr_v2c_val < 0;
						MS_DATA curr_v2c_abs  = (curr_v2c_val < 0.0 ? -curr_v2c_val : curr_v2c_val); //shifting 
						curr_v2c_abs -= beta;
						if( curr_v2c_abs < 0 )
						{
							curr_v2c_abs = 0;
							curr_v2c_sgn = 0;
						}
#endif
						signs[memOffset+n] = curr_v2c_sgn;

						process_check_node( &curr[n], curr_v2c_sgn, curr_v2c_abs, k );
					}  
				}

				memOffset += m_ldpc;
			}

			for( k = 0; k < m_ldpc; k++ )
				prev[stateOffset+k] = curr[k];
		}

		for( k = 0; k < n_ldpc; k++ )
			soft[k] = y[k];
		
		// STATE: compute soft output
		for( j = 0; j < rh; j++ ) 
		{
			int stateOffset = j*m_ldpc;
			for( k = 0; k < nh; k++ )
			{
				int memOffset = j * n_ldpc + k * m_ldpc;
				int circ = matr[j][k];

				if( circ != SKIP )
				{
					for( n = 0; n < m_ldpc; n++ )
					{
						MS_DATA tmp = prev[stateOffset+n].pos == k ? prev[stateOffset+n].min2 : prev[stateOffset+n].min1;

						buffer[n] = (signs[memOffset+n] ^ prev[stateOffset+n].sign) ? -tmp : tmp;
					}

					rotate( buffer, rbuffer, m_ldpc - circ, sizeof(buffer[0]), m_ldpc ); 

					for( n = 0; n < m_ldpc; n++ )
					{
						soft[k * m_ldpc + n] = soft[k * m_ldpc + n] + rbuffer[n];
					}

				}
			}
		} 

		// check syndrome
		check_syndrome( matr, rh, nh, soft, rsoft, m_ldpc, synd );  

		for( parity = 0, i = 0; i < r_ldpc; i++ )
			parity |= synd[i];


		if( parity == 0 )
			break; 

	} // stage

#endif
	
	for( k = 0; k < n_ldpc; k++ )
		decword[k] = soft[k] < 0;

	return parity ? -iter : iter+1; 
}
#endif



int imin_sum_decod_qc_lm( DEC_STATE* st, int y[], int decword[], int maxsteps, double alpha, double thr, int qbits, int dbits )
{
	int   i, j, k, n;
	int iter;
	int parity;
	short **matr        = st->hd;
	short *synd         = st->syndr; 
	IMS_DATA *soft      = st->ims_soft;
	short *BnNS         = st->ims_BnNS;  
	IMS_DATA *buffer    = st->ims_buffer;
	IMS_DATA *rbuffer   = st->ims_rbuffer;
	IMS_DATA *rsoft     = st->ims_rsoft;
	IMS_DEC_STATE *dcs  = st->ims_dcs;		
	IMS_DEC_STATE *tmps = st->ims_tmps;		
	IMS_DATA *iy        = st->ims_y;
	IMS_DATA max_data   = (IMS_DATA)((1L << (dbits-1)) - 1);
	IMS_DATA max_quant  = (IMS_DATA)((1L << (qbits-1)) - 1);


	int m  = st->m;
	int rh = st->rh;
	int nh = st->nh;
	
	int c_ldpc = m;
	int r_ldpc = m * rh;
	int n_ldpc = m * nh;

#ifdef MS_MUL_CORRECTION
	int ialpha = (int)(alpha * (1L << MS_ALPHA_FPP));
#else
	IMS_DATA delta = (IMS_DATA)(max_data * alpha);
#endif

	for( i = 0; i < r_ldpc; i++ )
	{
		dcs[i].min1 = 0;
		dcs[i].min2 = 0;
		dcs[i].pos  = 0;
		dcs[i].sign = 0;
	}


	{
		double coef;
		double en = 0;
	
		for( i = 0; i < n_ldpc; i++ )
			en += y[i] * y[i];

		coef = sqrt(n_ldpc / en);

		for( i = 0; i < n_ldpc; i++ )
		{
			int ival;
			double val = y[i];
			int   sign = 0;

			if( val < 0 )
			{
				val = -val;
				sign = 1;
			}

			val *= coef;
			if( val > thr ) 
				val = thr;
            ival  = (short)floor( val * max_quant / thr + 0.5 );

			iy[i] = sign ? -ival : ival;
		}
	}

	memset( BnNS, 0, rh*n_ldpc*sizeof(BnNS[0]) );
/*
	// check input codeword
	for( j = 0; j < rh; j++ )
	{
		for( i = 0; i < nh; i++ ) 
		{
			int pos_r = j * m;
			int pos_n = i * m;
			int circ = matr[j][i];

			if( circ != - 1 )
			{
				rotate( &y[pos_n], rsoft, circ, sizeof(y[0]), m );

				for( k = 0; k < m; k++ )
					synd[pos_r + k] ^= rsoft[k] < 0;
			}
		}
	}

	parity = 0;
	for( i = 0; i < r_ldpc; i++ )
		parity |= synd[i];
*/

	for( iter = 0; iter < maxsteps; iter++ )
	{

		int memOffset;


		// INIT_STAGE
		memset( synd, 0, r_ldpc*sizeof(synd[0]) );
		memset( soft, 0, n_ldpc*sizeof(soft[0]) );
		memOffset = 0;

		// STATE 1: compute sum
		for( j = 0; j < rh; j++ ) 
		{
			int stateOffset = j*c_ldpc;

			for( k = 0; k < nh; k++ )
			{
				int circ = matr[j][k];

				if( circ != SKIP )
				{
					for( n = 0; n < c_ldpc; n++ )
					{
						IMS_DATA tmp = dcs[stateOffset+n].pos == k ? dcs[stateOffset+n].min2 : dcs[stateOffset+n].min1;
#ifdef MS_MUL_CORRECTION
						tmp = (IMS_DATA)((tmp * ialpha) >> MS_ALPHA_FPP);
#else
						tmp -= delta;
						if( tmp < 0 ) tmp = 0;
#endif

						buffer[n] = (BnNS[memOffset+n] ^ dcs[stateOffset+n].sign) ? -tmp : tmp;
					}

					rotate( buffer, rbuffer, c_ldpc - circ, sizeof(buffer[0]), c_ldpc ); 

					for( n = 0; n < c_ldpc; n++ )
					{
						IMS_DATA t = soft[k * c_ldpc + n] + rbuffer[n];
						soft[k * c_ldpc + n] = limit_val( t, max_data );
					}

				}

				memOffset += c_ldpc;
			}

		} 

		// STATE 2: 
		for( k = 0; k < n_ldpc; k++ )
		{
#if 0
			soft[k] = iy[k] + (IMS_DATA)((soft[k] * ialpha) >> MS_ALPHA_FPP); //alpha-normalization
#else
			soft[k] = iy[k] + soft[k]; 
#endif
			soft[k] = limit_val( soft[k], max_data );
			decword[k]   = soft[k] < 0;
		}



		// STATE 3: update statistic
		memOffset = 0;
		for( j = 0; j < rh; j++ )
		{
			int stateOffset = j*c_ldpc;

			for( k = 0; k < c_ldpc; k++ )
			{
				tmps[k].min1 = max_data;
				tmps[k].min2 = max_data;
				tmps[k].pos  = 0;
				tmps[k].sign = 0;
			}

			for( k = 0; k < nh; k++ )
			{
				int circ = matr[j][k];

				if( circ != SKIP )
				{
					rotate( &soft[k*c_ldpc], rsoft, circ, sizeof(soft[0]), c_ldpc ); 

					for( n = 0; n < c_ldpc; n++ )
						synd[stateOffset+n] ^= rsoft[n] < 0;

					for( n = 0; n < c_ldpc; n++ )
						buffer[n] = dcs[stateOffset+n].pos == k ? dcs[stateOffset+n].min2 : dcs[stateOffset+n].min1;

					for( n = 0; n < c_ldpc; n++ )
					{
						short sign   = BnNS[memOffset+n] ^ dcs[stateOffset+n].sign;
#ifdef MS_MUL_CORRECTION
						IMS_DATA val = (IMS_DATA)((buffer[n] * ialpha) >> MS_ALPHA_FPP);
#else
						IMS_DATA tmp = buffer[n] - delta;
						IMS_DATA val = tmp < 0 ? 0 : tmp;
#endif
						IMS_DATA t   = sign ? -val : val;
                        IMS_DATA v2c_msg = rsoft[n] - t;
						short v2c_sign = v2c_msg < 0;
						
						BnNS[memOffset+n] = v2c_sign;
						tmps[n].sign     ^= v2c_sign;

						val = v2c_msg < 0 ? -v2c_msg : v2c_msg;
						val = (val > max_data) ? max_data : val;


						if( val < tmps[n].min1 )      
						{
							tmps[n].pos = k;
							tmps[n].min2 = tmps[n].min1; 
							tmps[n].min1 = val;
						}
						else
						{
							if( val < tmps[n].min2 )
								tmps[n].min2 = val;
						}

					}  
				}


				memOffset += c_ldpc;
			}

			for( k = 0; k < c_ldpc; k++ )
				dcs[stateOffset+k] = tmps[k];

		}

		for( parity = 0, i = 0; i < r_ldpc; i++ )
			parity |= synd[i];


		if( parity == 0 )
			break; 

	} // stage

	return parity ? -iter : iter+1; 
}



