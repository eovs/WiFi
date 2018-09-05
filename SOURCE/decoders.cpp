#include <string.h>
#include <math.h>
//#include <assert.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include "matrix.h"
#include "matrix_ext.h"
#include "decoders.h"

//#define CHINA_VERSION

char const * const DEC_FULL_NAME[] = 
{
	"Layered Min-Sum",
	"Integer Layered Min-Sum",
	"Layered Sum-Prod"
	"Integer Layered Sum-Prod"
};

#define LOG_FPP 12
#define ONE_LOG	(1 << LOG_FPP)

#define ILOG(x) ( (int)((x)*ONE_LOG + 0.5) )


#define IL_SOFT_ONE (1 << IL_SOFT_FPP)



#define ILCHE_SOFT_FPP	8 // bits per soft
#define ONE_ILCHE_SOFT (1 << ILCHE_SOFT_FPP)
#define ILCHE_INNER_FPP 8 // bits per inner variables
#define ONE_ILCHE_INNER (1 << ILCHE_INNER_FPP)



#define div_power2( x, n )    ((x) >> (n))
#define div_power2r( x, n )   (((x) + (1 << (n-1))) >> (n))

#define MY_VERSION
//#define ALPHA 0.8
//#define BETA  0.35 

#define ALPHA 1.0
#define BETA  0.4



#define INPUT_LIMIT 20.0

#define SP_DEC_MIN_VAL 0.000001
#define SP_DEC_MAX_VAL (1.0 - SP_DEC_MIN_VAL)


#define maxi( a, b )  (a) < (b) ? (b) : (a)
#define mini( a, b )  (b) < (a) ? (b) : (a)

#define limit_val( x, max_val )	(x) > max_val ? max_val : ((x) < -max_val ? -max_val : (x))

static double mind(double a, double b) {if (a<b) return a; else return b;}
static double maxd(double a, double b) {if (a<b) return b; else return a;}
static double absd (double x)
{
	if (x<0) return -x; else return x;
}


static void update_statistics( double *prev_soft, double *curr_soft, double *sign_counter, double *min_abs_llr, double thr, int n )
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


static double** Alloc2d_double( int b, int c )
{
	double **p;
	int i;

//	p = (double**)calloc( b, sizeof(double*) );
	p = (double**)malloc( b * sizeof(double*) );
	//assert(p);
//	p[0] = (double*)calloc( b*c, sizeof(double) );
	p[0] = (double*)malloc( b*c*sizeof(double) );
	for( i = 1; i < b; i++ )
	{
		p[i] = p[i-1] + c;
		//assert( p[i] );
	}
	return p;
}


static int** Alloc2d_int( int b, int c )
{
	int **p;
	int i;

//	p = (int**)calloc( b, sizeof(int*) );
	p = (int**)malloc( b * sizeof(int*) );
	//assert(p);
//	p[0] = (int*)calloc( b*c, sizeof(int) );
	p[0] = (int*)malloc( b*c*sizeof(int) );
	for( i = 1; i < b; i++ )
	{
		p[i] = p[i-1] + c;
		//assert( p[i] );
	}
	return p;
}


static void free2d_int( int **p )
{
	free( p[0] );
	free( p );
}

static void free2d_double( double **p )
{
	free( p[0] );
	free( p );
}
/*
void rotate_sign( int x[], int y[], int shift, int M )
{
	while( shift < 0 )
		shift += M;

	while( shift >= M )
		shift -= M;

	memcpy( y, x + shift, (M-shift)*sizeof(y[0]));
	memcpy( y+(M-shift), x, shift*sizeof(y[0]));
}
*/
static void rotate_data( double x[], double y[], int shift, int M )
{
	while( shift < 0 )
		shift += M;

	while( shift >= M )
		shift -= M;

	memcpy( y, x + shift, (M-shift)*sizeof(y[0]));
	memcpy( y+(M-shift), x, shift*sizeof(y[0]));
}

static void rotate( void *x, void *y, int shift, int size, int M )
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

DEC_STATE* decod_open( int codec_id, int mh, int nh, int M )
{
    DEC_STATE* st;
    
	int N = nh * M;
	int R = mh * M;
	int BBsize = (mh > nh - mh) ? mh : nh - mh;


//    st = (DEC_STATE*)calloc( 1, sizeof(DEC_STATE) );
    st = (DEC_STATE*)malloc( sizeof(DEC_STATE) );
	memset( st, 0, sizeof(DEC_STATE) );
    if( !st ) 
        return NULL;

	st->nh      = nh;
	st->rh      = mh;
	st->m       = M;
    st->n       = N;
    st->codec_id = codec_id;


	st->hd = Alloc2d_int( mh, nh );
	if( st->hd==NULL)
		return NULL;

#ifdef KEEP_STATISTIC
//	st->sign_counter = (double*)calloc(N, sizeof(st->sign_counter[0]));
//	st->min_abs_llr = (double*)calloc(N, sizeof(st->min_abs_llr[0]));
//	st->prev_soft = (double*)calloc(N, sizeof(st->prev_soft[0]));
#endif


//	st->syndr = (int*)calloc(R, sizeof(st->syndr[0]) );
	st->syndr = (int*)malloc(R * sizeof(st->syndr[0]) );
	if( !st->syndr )
		return NULL;
    

	switch( codec_id )
	{
	case LMS_DEC:
//		st->lms_soft = (MS_DATA*)calloc(N, sizeof(st->lms_soft[0]) );
		st->lms_soft = (MS_DATA*)malloc(N * sizeof(st->lms_soft[0]) );
		if( !st->lms_soft )
			return NULL;

//		st->lms_BnNS = (int*)calloc(mh*N, sizeof( st->lms_BnNS[0] ) );
		st->lms_BnNS = (int*)malloc(mh*N* sizeof( st->lms_BnNS[0] ) );
		if( !st->lms_BnNS )
			return NULL;

//		st->lms_dcs = (MS_DEC_STATE*)calloc(R, sizeof( st->lms_dcs[0] ) );
		st->lms_dcs = (MS_DEC_STATE*)malloc(R* sizeof( st->lms_dcs[0] ) );
		if( !st->lms_dcs )
			return NULL;

//		st->lms_tmps = (MS_DEC_STATE*)calloc(M, sizeof( st->lms_tmps[0] ) );
		st->lms_tmps = (MS_DEC_STATE*)malloc(M* sizeof( st->lms_tmps[0] ) );
		if( !st->lms_tmps )
			return NULL;

//		st->lms_buffer = (MS_DATA*)calloc(M, sizeof(st->lms_buffer[0]) );
		st->lms_buffer = (MS_DATA*)malloc(M* sizeof(st->lms_buffer[0]) );
		if( !st->lms_buffer )
			return NULL;

//		st->lms_rbuffer = (MS_DATA*)calloc(M, sizeof(st->lms_rbuffer[0]) );
		st->lms_rbuffer = (MS_DATA*)malloc(M* sizeof(st->lms_rbuffer[0]) );
		if( !st->lms_rbuffer )
			return NULL;

//		st->lms_rsoft = (MS_DATA*)calloc(M, sizeof(st->lms_rsoft[0]) );
		st->lms_rsoft = (MS_DATA*)malloc(M* sizeof(st->lms_rsoft[0]) );
		if( !st->lms_rsoft )
			return NULL;
		break;

	case ILMS_DEC:
//		st->ilms_y = (IMS_DATA*)calloc(N, sizeof(st->ilms_y[0]) );
		st->ilms_y = (IMS_DATA*)malloc(N* sizeof(st->ilms_y[0]) );
		if( !st->ilms_y )
			return NULL;

//		st->ilms_decword = (IMS_DATA*)calloc(N, sizeof(st->ilms_decword[0]) );
		st->ilms_decword = (IMS_DATA*)malloc(N* sizeof(st->ilms_decword[0]) );
		if( !st->ilms_decword )
			return NULL;

//		st->ilms_soft = (IMS_DATA*)calloc(N, sizeof(st->ilms_soft[0]) );
		st->ilms_soft = (IMS_DATA*)malloc(N* sizeof(st->ilms_soft[0]) );
		if( !st->ilms_soft )
			return NULL;

//		st->ilms_BnNS = (int*)calloc(mh*N, sizeof( st->ilms_BnNS[0] ) );
		st->ilms_BnNS = (int*)malloc(mh*N* sizeof( st->ilms_BnNS[0] ) );
		if( !st->ilms_BnNS )
			return NULL;

//		st->ilms_dcs = (IMS_DEC_STATE*)calloc(R, sizeof( st->ilms_dcs[0] ) );
		st->ilms_dcs = (IMS_DEC_STATE*)malloc(R* sizeof( st->ilms_dcs[0] ) );
		if( !st->ilms_dcs )
			return NULL;

//		st->ilms_tmps = (IMS_DEC_STATE*)calloc(M, sizeof( st->ilms_tmps[0] ) );
		st->ilms_tmps = (IMS_DEC_STATE*)malloc(M* sizeof( st->ilms_tmps[0] ) );
		if( !st->ilms_tmps )
			return NULL;

//		st->ilms_buffer = (IMS_DATA*)calloc(M, sizeof(st->ilms_buffer[0]) );
		st->ilms_buffer = (IMS_DATA*)malloc(M* sizeof(st->ilms_buffer[0]) );
		if( !st->ilms_buffer )
			return NULL;

//		st->ilms_rbuffer = (IMS_DATA*)calloc(M, sizeof(st->ilms_rbuffer[0]) );
		st->ilms_rbuffer = (IMS_DATA*)malloc(M* sizeof(st->ilms_rbuffer[0]) );
		if( !st->ilms_rbuffer )
			return NULL;

//		st->ilms_rsoft = (IMS_DATA*)calloc(M, sizeof(st->ilms_rsoft[0]) );
		st->ilms_rsoft = (IMS_DATA*)malloc(M* sizeof(st->ilms_rsoft[0]) );
		if( !st->ilms_rsoft )
			return NULL;
		break;


	case LCHE_DEC:
//		st->lche_data0 = (double*)calloc(N, sizeof( st->lche_data0[0] ) );
		st->lche_data0 = (double*)malloc(N* sizeof( st->lche_data0[0] ) );
		if(!st->lche_data0 )
			return NULL;

//		st->lche_tmp = (double*)calloc(nh, sizeof( st->lche_tmp[0] ) );
		st->lche_tmp = (double*)malloc(nh* sizeof( st->lche_tmp[0] ) );
		if(!st->lche_tmp )
			return NULL;

//		st->lche_soft_out = (double*)calloc(N, sizeof(st->lche_soft_out[0]) );
		st->lche_soft_out = (double*)malloc(N* sizeof(st->lche_soft_out[0]) );
		if( !st->lche_soft_out )
			return NULL;

		st->lche_state = Alloc2d_double( mh*M, nh*M );
		if( st->lche_state==NULL)
			break;
		break;

	case ILCHE_DEC:
//		st->ilche_data0 = (int*)calloc(N, sizeof( st->ilche_data0[0] ) );
		st->ilche_data0 = (int*)malloc(N* sizeof( st->ilche_data0[0] ) );
		if(!st->ilche_data0 )
			return NULL;

//		st->ilche_tmp = (int*)calloc(nh, sizeof( st->ilche_tmp[0] ) );
		st->ilche_tmp = (int*)malloc(nh* sizeof( st->ilche_tmp[0] ) );
		if(!st->ilche_tmp )
			return NULL;

//		st->ilche_soft_out = (int*)calloc(N, sizeof(st->ilche_soft_out[0]) );
		st->ilche_soft_out = (int*)malloc(N* sizeof(st->ilche_soft_out[0]) );
		if( !st->ilche_soft_out )
			return NULL;

		st->ilche_state = Alloc2d_int( mh*M, nh*M );
		if( st->ilche_state==NULL)
			break;
		break;

	default: return NULL;
	}

    return st;
    
}

static void check_syndrome( int **matr, int rh, int nh, MS_DATA *soft, MS_DATA *rsoft, int m_ldpc, int *synd )
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

void icheck_syndrome( int **matr, int rh, int nh, IMS_DATA *soft, IMS_DATA *rsoft, int m_ldpc, int *synd )
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
	int codec_id = st->codec_id;
	int R = mh * M;



	if(st->hd)		{ free2d_int( st->hd );          st->hd       = NULL; }
#ifdef KEEP_STATISTIC
	if(st->sign_counter)	{ free(st->sign_counter);                     st->sign_counter        = NULL; }
	if(st->min_abs_llr)		{ free(st->min_abs_llr);                     st->min_abs_llr        = NULL; }
	if(st->prev_soft)		{ free(st->prev_soft);                     st->prev_soft        = NULL; }
#endif



	if(st->syndr)	{ free( st->syndr ); st->syndr = NULL; }

	switch( codec_id )
	{
	case LMS_DEC:
		if( st->lms_soft )    { free(st->lms_soft);    st->lms_soft    = NULL; }
		if( st->lms_BnNS )    { free(st->lms_BnNS);    st->lms_BnNS    = NULL; }
		if( st->lms_dcs )     { free(st->lms_dcs);     st->lms_dcs     = NULL; }
		if( st->lms_tmps )    { free(st->lms_tmps);    st->lms_tmps    = NULL; }
		if( st->lms_buffer )  { free(st->lms_buffer);  st->lms_buffer  = NULL; }
		if( st->lms_rbuffer ) { free(st->lms_rbuffer); st->lms_rbuffer = NULL; }
		if( st->lms_rsoft )   { free(st->lms_rsoft);   st->lms_rsoft   = NULL; }
		break;

	case ILMS_DEC:
		if( st->ilms_y )       { free(st->ilms_y);       st->ilms_y       = NULL; }
		if( st->ilms_decword ) { free(st->ilms_decword); st->ilms_decword = NULL; }
		if( st->ilms_soft )    { free(st->ilms_soft);    st->ilms_soft    = NULL; }
		if( st->ilms_BnNS )    { free(st->ilms_BnNS);    st->ilms_BnNS    = NULL; }
		if( st->ilms_dcs )     { free(st->ilms_dcs);     st->ilms_dcs     = NULL; }
		if( st->ilms_tmps )    { free(st->ilms_tmps);    st->ilms_tmps    = NULL; }
		if( st->ilms_buffer )  { free(st->ilms_buffer);  st->ilms_buffer  = NULL; }
		if( st->ilms_rbuffer ) { free(st->ilms_rbuffer); st->ilms_rbuffer = NULL; }
		if( st->ilms_rsoft )   { free(st->ilms_rsoft);   st->ilms_rsoft   = NULL; }
		break;

	case LCHE_DEC:
		if( st->lche_data0 )     { free(st->lche_data0);            st->lche_data0    = NULL; }
		if( st->lche_tmp )     { free(st->lche_tmp);		        st->lche_tmp      = NULL; }
		if( st->lche_soft_out )  { free(st->lche_soft_out);         st->lche_soft_out = NULL; }
		if( st->lche_state )     { free2d_double( st->lche_state ); st->lche_state    = NULL; }
		break;

	case ILCHE_DEC:
		if( st->ilche_data0 )     { free(st->ilche_data0);            st->ilche_data0    = NULL; }
		if( st->ilche_tmp )       { free(st->ilche_tmp);	          st->ilche_tmp      = NULL; }
		if( st->ilche_soft_out )  { free(st->ilche_soft_out);         st->ilche_soft_out = NULL; }
		if( st->ilche_state )     { free2d_int( st->ilche_state );    st->ilche_state    = NULL; }
		break;
	default:;
	}

    free( st );
}



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

static MS_DATA logexp_int( MS_DATA x )
{
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




static void map_bin_llr(MS_DATA y[], int n )
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
	int **hd = st->hd;
	double **Z = st->lche_state;
	double *soft_out = st->lche_soft_out;
	int *syndr    = st->syndr;
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




static int d2i( double val, int one_fpp )
{
	double absval = val < 0.0 ? -val : val;
	int    sign   = val < 0.0 ? 1 : 0;
	int iabsval = (int)(absval * one_fpp + 0.5);
	return sign ? -iabsval : iabsval;
}


static int iA[] = 
{	
	ILOG(1.41e+00), ILOG(7.72e-01), ILOG(4.54e-01), ILOG(2.72e-01), ILOG(1.65e-01), ILOG(9.97e-02), ILOG(6.04e-02), ILOG(3.66e-02), 
	ILOG(2.22e-02), ILOG(1.35e-02), ILOG(8.17e-03), ILOG(4.96e-03), ILOG(3.01e-03), ILOG(1.82e-03), ILOG(1.11e-03), ILOG(6.71e-04),
	ILOG(4.07e-04), ILOG(2.47e-04), ILOG(1.50e-04), ILOG(9.08e-05), ILOG(5.51e-05), ILOG(3.34e-05), ILOG(2.03e-05), ILOG(1.23e-05), 
	ILOG(7.45e-06), ILOG(4.52e-06), ILOG(2.74e-06), ILOG(1.66e-06), ILOG(1.01e-06), ILOG(6.12e-07), ILOG(3.71e-07), ILOG(2.25e-07)
};

static int iB[] = 
{
	ILOG(3.47), ILOG(2.77), ILOG(2.37), ILOG(2.08), ILOG(1.86), ILOG(1.69), ILOG(1.54), ILOG(1.41),
	ILOG(1.29), ILOG(1.19), ILOG(1.11), ILOG(1.03), ILOG(0.95), ILOG(0.89), ILOG(0.83), ILOG(0.77),
	ILOG(0.72), ILOG(0.67), ILOG(0.63), ILOG(0.59), ILOG(0.55), ILOG(0.52), ILOG(0.48), ILOG(0.45),
	ILOG(0.43), ILOG(0.40), ILOG(0.37), ILOG(0.35), ILOG(0.33), ILOG(0.31), ILOG(0.29), ILOG(0.27)
};

static int iC[] =
{
	ILOG(6.93), ILOG(6.24), ILOG(5.83), ILOG(5.55), ILOG(5.32), ILOG(5.14), ILOG(4.99), ILOG(4.85),
	ILOG(4.73), ILOG(4.63), ILOG(4.53), ILOG(4.45), ILOG(4.37), ILOG(4.29), ILOG(4.22), ILOG(4.16),
	ILOG(4.10), ILOG(4.04), ILOG(3.99), ILOG(3.94), ILOG(3.89), ILOG(3.84), ILOG(3.80), ILOG(3.75),
	ILOG(3.71), ILOG(3.67), ILOG(3.64), ILOG(3.60), ILOG(3.56), ILOG(3.53), ILOG(3.50), ILOG(3.47)
};

static MS_DATA ilogexp_int( int ix )
{
	if( ix <= 0 ) ix = ONE_LOG/4096;
	if( ix > 16.0*ONE_LOG ) ix = ONE_LOG * 16;

	if( ix >= ONE_LOG * 2 )
		return  -iA[ (ix * 2) / ONE_LOG - 1]; 
	else
		if( ix > ONE_LOG / 16 )
			return -iB[ (ix * 16) / ONE_LOG  - 1];
		else
			if( ix > ONE_LOG / 512 )
				return -iC[ (ix * 512) / ONE_LOG - 1];
			else  
			{
                int s = 0;
				while( ix < ONE_LOG/512 )
				{
					ix *= 32;
					s -= (int)(3.46 * ONE_LOG);
				}

				return s - iC[(ix * 512) / ONE_LOG - 1];
			}
}

static void imap_bin_llr( int y[], int n )
{
	int i;
	static int hard[1024];
	static int ay[1024];
	static int alogpy[1024];
	static int a[1024];
	static int A[1024];
	int synd;
	int sum;

	for( i = 0; i < n; i++ )
	{
		hard[i] = y[i] < 0;
	}	
	
	synd = 0;
	for( i = 0; i < n; i++ )
		synd ^=hard[i];

	for( i = 0; i < n; i++ )
		hard[i] ^= synd;

#if LOG_FPP >= ILCHE_INNER_FPP
	for( i = 0; i < n; i++ )
		ay[i] = (y[i] < 0.0 ? -y[i] : y[i]) << (LOG_FPP - ILCHE_INNER_FPP);
#else
	for( i = 0; i < n; i++ )
		ay[i] = div_power2r(y[i] < 0.0 ? -y[i] : y[i]), ILCHE_INNER_FPP-LOG_FPP);
#endif
	for( i = 0; i < n; i++ )
		alogpy[i] = d2i( ilogexp_int( ay[i] ), 1 );

	sum = 0;
	for( i = 0; i < n; i++ )
		sum += alogpy[i];

	for( i = 0; i < n; i++ )
		A[i] = alogpy[i] - sum;

	for( i = 0; i < n; i++ )
		a[i] = d2i( ilogexp_int( A[i] ) /  ONE_LOG, ONE_ILCHE_INNER );

	for( i = 0; i < n; i++ )
		y[i] = (2 * hard[i] - 1) * a[i];
}

int ilche_decod( DEC_STATE* st, int soft[], int decword[], int maxsteps )
{ 
	int i, j, k;
	int steps;
	int **hd = st->hd;
	int **Z = st->ilche_state;
	int *soft_out = st->ilche_soft_out;
	int *syndr    = st->syndr;
	int *data0 = st->ilche_data0;
	int *u = st->ilche_data0;
	int *y = st->ilche_tmp;
	int synd;

	int rh = st->rh;
	int nh = st->nh;
	int m  = st->m;
	int r = rh * m;
	int n = nh * m;

	for( i = 0; i < r; i++ )
		for(j = 0; j < n; j++ )
			Z[i][j] = 0;

	// just to compute syndrome before iterations
	for( i = 0; i < n; i++ )
		soft_out[i] = soft[i] * ONE_ILCHE_SOFT;

#ifdef KEEP_STATISTIC
	for( i = 0; i < n; i++ ) st->prev_soft[i] = soft_out[i];
	for( i = 0; i < n; i++ ) st->sign_counter[i] = 0;
	for( i = 0; i < n; i++ ) st->min_abs_llr[i]  = 10000;
#endif

	memset( syndr, 0, r * sizeof(syndr[0]) );
	icheck_syndrome( hd, rh, nh, soft_out, data0, m, syndr );

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
				int *a = Z[i * m + k];

				cnt = 0;
				for( j = 0; j < nh; j++ )
				{
					int circ = hd[i][j];

					if( circ != -1 )
					{
						int idx = j * m + ((k + circ) % m);
#if ONE_ILCHE_INNER < ONE_ILCHE_SOFT
						y[cnt] = u[cnt] = div_power2r( soft_out[idx], (ONE_ILCHE_SOFT - ONE_ILCHE_INNER) ) - a[idx];
#else
						y[cnt] = u[cnt] = ( soft_out[idx] <<  (ONE_ILCHE_INNER - ONE_ILCHE_SOFT) ) - a[idx];
#endif
						cnt++;
					}
				}

				imap_bin_llr( u, cnt );

				cnt = 0;
				for( j = 0; j < nh; j++ )
				{
					int circ = hd[i][j];

					if( circ != -1 )
					{
						int idx = j * m + ((k + circ) % m);
						{
#if ONE_ILCHE_INNER > ONE_ILCHE_SOFT
							soft_out[idx] = div_power2r( u[cnt] + y[cnt], (ONE_ILCHE_INNER - ONE_ILCHE_SOFT);
#else
							soft_out[idx] = (u[cnt] + y[cnt]) << (ONE_ILCHE_SOFT - ONE_ILCHE_INNER);
#endif
							a[idx]        = u[cnt];
						}
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
		icheck_syndrome( hd, rh, nh, soft_out, data0, m, syndr );
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









static void process_check_node( MS_DEC_STATE *curr, int curr_v2c_sign, MS_DATA abs_curr_v2c, int index )
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

static void iprocess_check_node( IMS_DEC_STATE *curr, int curr_v2c_sign, IMS_DATA abs_curr_v2c, int index )
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


static MS_DATA	get_c2v_val( int sign, int pos, MS_DEC_STATE *stat, double alpha, double beta )
{
	MS_DATA c2v_abs  = stat->pos == pos ? stat->min2 : stat->min1;

	c2v_abs  = c2v_abs * alpha - beta; 
	if( c2v_abs < 0 )
		c2v_abs = 0;

	return sign ? -c2v_abs : c2v_abs;
}

static IMS_DATA	iget_c2v_val( int sign, int pos, IMS_DEC_STATE *stat, double alpha, int beta )
{
	IMS_DATA c2v_abs  = stat->pos == pos ? stat->min2 : stat->min1;

	c2v_abs  = (IMS_DATA)(c2v_abs * alpha - beta); 
	if( c2v_abs < 0 )
		c2v_abs = 0;

	return sign ? -c2v_abs : c2v_abs;
}

static double beta_control( MS_DEC_STATE *state )
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

const int betaaa = (int)(0.5 * IL_SOFT_ONE);
static IMS_DATA ibeta_control( IMS_DEC_STATE *state )
{
#if 0

//	if( state->min2 > state->min1 * 1.1 )	beta = 0.4;
//	if( state->min2 > state->min1 * 1.2 )	beta = 0.7;


	//if( state->min2 - state->min1 > 0.1 )	beta = 0.6;
	if( state->min2 - state->min1 > 0.2 )	beta = 0.70;
	if( state->min2 - state->min1 > 0.5 )	beta = 0.75;
	
#endif
	return betaaa;
}

#ifdef CHINA_VERSION
int lmin_sum_decod_qc_lm( DEC_STATE* st, int y[], int decword[], int maxsteps, double alpha, double beta )
{
	int   i, j, k, n;
	int iter;
	int parity;
	int **matr         = st->hd;
	int *synd          = st->syndr; 
	MS_DATA *soft      = st->lms_soft;
	int *signs         = st->lms_BnNS;  
	MS_DATA *buffer    = st->lms_buffer;
	MS_DATA *rbuffer   = st->lms_rbuffer;
	MS_DATA *rsoft     = st->lms_rsoft;
	MS_DEC_STATE *prev = st->lms_dcs;		
	MS_DEC_STATE *curr = st->lms_tmps;		

	static MS_DATA curr_v2c[100000];

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
				curr[k].min1 = 10000;
				curr[k].min2 = 10000;
				curr[k].pos  = 0;
				curr[k].sign = 0;
			}


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
						int	prev_c2v_sgn  = signs[memOffset+n] ^ prev[stateOffset+n].sign;
						MS_DATA	prev_c2v_val  = prev_c2v_sgn ? -prev_c2v_abs[n] : prev_c2v_abs[n];
						MS_DATA curr_v2c_val  = rsoft[n] - prev_c2v_val;

						int	curr_v2c_sgn = curr_v2c_val < 0;
						MS_DATA curr_v2c_abs = (curr_v2c_val < 0.0 ? -curr_v2c_val : curr_v2c_val); //shifting 
						curr_v2c_abs -= beta;

						//curr_v2c_sgn = curr_v2c_abs < 0 ? 0 : curr_v2c_sgn;
						curr_v2c_abs = curr_v2c_abs < 0 ? 0 : curr_v2c_abs;

#if 1
						curr_v2c[k * m_ldpc + n] = curr_v2c_val;
#else
						curr_soft[n] = curr_v2c_val;
#endif
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

#if 1
					for( n = 0; n < m_ldpc; n++ )
						buffer[n] = curr_v2c[k * m_ldpc + n] + curr_c2v_val[n];

					rotate( buffer, rbuffer, m_ldpc - circ, sizeof(curr_c2v_val[0]), m_ldpc ); 
					{
						double omega = 1.0/16.0;
						double omega1 = 1 + omega;
						for( n = 0; n < m_ldpc; n++ )
							curr_soft[n] = rbuffer[n] * omega1 - curr_soft[n] * omega;
					}
#else
					for( n = 0; n < m_ldpc; n++ )
						buffer[n] = curr_soft[n] + curr_c2v_val[n];

					rotate( buffer, rbuffer, m_ldpc - circ, sizeof(curr_c2v_val[0]), m_ldpc ); 

					for( n = 0; n < m_ldpc; n++ )
						curr_soft[n] = rbuffer[n];
#endif

				}
			}
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

	
	for( k = 0; k < n_ldpc; k++ )
		decword[k] = soft[k] < 0;

	return parity ? -iter : iter+1; 
}
#else
int lmin_sum_decod_qc_lm( DEC_STATE* st, int y[], int decword[], int maxsteps, double alpha, double beta )
{
	int   i, j, k, n;
	int iter;
	int parity;
	int **matr         = st->hd;
	int *synd          = st->syndr; 
	MS_DATA *soft      = st->lms_soft;
	int *signs         = st->lms_BnNS;  
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
				curr[k].min1 = 10000;
				curr[k].min2 = 10000;
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
						int	prev_c2v_sgn  = signs[memOffset+n] ^ prev[stateOffset+n].sign;
						MS_DATA	prev_c2v_val  = prev_c2v_sgn ? -prev_c2v_abs[n] : prev_c2v_abs[n];
						MS_DATA curr_v2c_val  = rsoft[n] - prev_c2v_val;
#if 0
						int	curr_v2c_sgn = curr_v2c_val < 0;
						MS_DATA curr_v2c_abs = (curr_v2c_val < 0.0 ? -curr_v2c_val : curr_v2c_val)*alpha; //scaling 
#else
						int	curr_v2c_sgn = curr_v2c_val < 0;
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
						int	prev_c2v_sgn  = signs[memOffset+n] ^ prev[stateOffset+n].sign;
#if 0
						MS_DATA	prev_c2v_val  = get_c2v_val( prev_c2v_sgn, k, &prev[stateOffset+n], 1.0/*alpha*/, beta );
#else
						MS_DATA beta          = beta_control( &prev[stateOffset+n] );
						MS_DATA	prev_c2v_val  = get_c2v_val( prev_c2v_sgn, k, &prev[stateOffset+n], 1.0/*alpha*/, beta );
#endif					
						MS_DATA curr_v2c_val  = rsoft[n] - prev_c2v_val;
						int	curr_v2c_sgn = curr_v2c_val < 0;
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
						int curr_c2v_sgn  = signs[memOffset+n] ^ prev[stateOffset+n].sign;
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
						int	prev_c2v_sgn = signs[memOffset+n] ^ prev[stateOffset+n].sign;
						MS_DATA	prev_c2v_val = prev_c2v_sgn ? -buffer[n] : buffer[n];
						MS_DATA curr_v2c_val = rsoft[n] - prev_c2v_val;
#if 0
						int	curr_v2c_sgn = curr_v2c_val < 0;
						MS_DATA curr_v2c_abs = (curr_v2c_val < 0.0 ? -curr_v2c_val : curr_v2c_val) * alpha; //scaling 
#else
						const MS_DATA beta = 0.4;
						int	curr_v2c_sgn = curr_v2c_val < 0;
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

void il_min_sum_reset( DEC_STATE *st )
{
	int   i;
	int *y = st->ilms_y;
	int *signs       = st->ilms_BnNS;  
	IMS_DEC_STATE *prev = st->ilms_dcs;		

	int m  = st->m;
	int rh = st->rh;
	int nh = st->nh;
	
	int r_ldpc = m * rh;
	int n_ldpc = m * nh;

	for( i = 0; i < r_ldpc; i++ )
	{
		prev[i].min1 = 0;
		prev[i].min2 = 0;
		prev[i].pos  = 0;
		prev[i].sign = 0;
	}

	memset( signs, 0, n_ldpc*rh*sizeof(signs[0]) );

}

int il_min_sum_iterate( DEC_STATE* st, int inner_data_bits )
{
	int i, j, k, n;
	int parity;
	int **matr         = st->hd;
	int *synd          = st->syndr; 
	IMS_DATA *soft     = st->ilms_soft;
	int *signs       = st->ilms_BnNS;  
	IMS_DATA *buffer    = st->ilms_buffer;
	IMS_DATA *rbuffer   = st->ilms_rbuffer;
	IMS_DATA *rsoft     = st->ilms_rsoft;
	IMS_DEC_STATE *prev = st->ilms_dcs;		
	IMS_DEC_STATE *curr = st->ilms_tmps;		

	int dmax = (1 << (inner_data_bits-1)) - 1;
	IMS_DATA ibeta = (IMS_DATA)(0.5 * IL_SOFT_ONE);// 0.4;
	
	int m  = st->m;
	int rh = st->rh;
	int nh = st->nh;
	
	int m_ldpc = m;
	int r_ldpc = m * rh;
	int n_ldpc = m * nh;

	// INIT_STAGE
	memset( synd, 0, r_ldpc*sizeof(synd[0]) );

	// update statistic
	for( j = 0; j < rh; j++ )
	{
		IMS_DATA *prev_c2v_abs = buffer;
		int stateOffset = j*m_ldpc;

		for( k = 0; k < m_ldpc; k++ )
		{
			curr[k].min1 = 10000;
			curr[k].min2 = 10000;
			curr[k].pos  = 0;
			curr[k].sign = 0;
		}


		for( k = 0; k < nh; k++ )
		{
			int memOffset = j * n_ldpc + k * m_ldpc;
			IMS_DATA *curr_soft    = &soft[k * m_ldpc];
			int circ = matr[j][k];

			if( circ != SKIP )
			{
				rotate( &soft[k*m_ldpc], rsoft, circ, sizeof(soft[0]), m_ldpc ); 

				for( n = 0; n < m_ldpc; n++ )
					prev_c2v_abs[n] = prev[stateOffset+n].pos == k ? prev[stateOffset+n].min2 : prev[stateOffset+n].min1;

				for( n = 0; n < m_ldpc; n++ )
				{
					int	prev_c2v_sgn  = signs[memOffset+n] ^ prev[stateOffset+n].sign;
					IMS_DATA	prev_c2v_val  = prev_c2v_sgn ? -prev_c2v_abs[n] : prev_c2v_abs[n];
					IMS_DATA curr_v2c_val  = rsoft[n] - prev_c2v_val;
					int	curr_v2c_sgn = curr_v2c_val < 0;
					IMS_DATA curr_v2c_abs = (curr_v2c_val < 0.0 ? -curr_v2c_val : curr_v2c_val); //shifting 
					curr_v2c_abs -= ibeta;

					//curr_v2c_sgn = curr_v2c_abs < 0 ? 0 : curr_v2c_sgn;
					curr_v2c_abs = curr_v2c_abs < 0 ? 0 : curr_v2c_abs;
					curr_v2c_abs = curr_v2c_abs > dmax ? dmax : curr_v2c_abs;

					/*curr_v2c*/curr_soft[n] = curr_v2c_val;
					signs[memOffset+n] = curr_v2c_sgn;

					iprocess_check_node( &curr[n], curr_v2c_sgn, curr_v2c_abs, k );

				}  
			}
		}

		for( k = 0; k < m_ldpc; k++ )
			prev[stateOffset+k] = curr[k];

		for( k = 0; k < nh; k++ )
		{
			IMS_DATA *curr_c2v_val = buffer;
			IMS_DATA *curr_soft    = &soft[k * m_ldpc];
			int memOffset = j * n_ldpc + k * m_ldpc;
			int circ = matr[j][k];

			if( circ != SKIP )
			{
				for( n = 0; n < m_ldpc; n++ )
				{
					IMS_DATA curr_c2v_abs = prev[stateOffset+n].pos == k ? prev[stateOffset+n].min2 : prev[stateOffset+n].min1;

					curr_c2v_val[n] = (signs[memOffset+n] ^ prev[stateOffset+n].sign) ? -curr_c2v_abs : curr_c2v_abs;
				}

				for( n = 0; n < m_ldpc; n++ )
					buffer[n] = /*curr_v2c*/curr_soft[n] + curr_c2v_val[n];

				rotate( buffer, rbuffer, m_ldpc - circ, sizeof(curr_c2v_val[0]), m_ldpc ); 

				for( n = 0; n < m_ldpc; n++ )
					curr_soft[n] = rbuffer[n];

			}
		}
	}

	// check syndrome
	icheck_syndrome( matr, rh, nh, soft, rsoft, m_ldpc, synd );  

	for( parity = 0, i = 0; i < r_ldpc; i++ )
		parity |= synd[i];

	if( parity == 0 )
		return 2; 
	else
		return 1;
}


#ifdef CHINA_VERSION
#define Z11111111
int il_min_sum_decod_qc_lm( DEC_STATE* st, int y[], int decword[], int maxsteps, double alpha, double beta, int inner_data_bits )
{
	int   i, j, k, n;
	int iter;
	int parity;
	int **matr         = st->hd;
	int *synd          = st->syndr; 
	IMS_DATA *soft     = st->ilms_soft;
	int *signs       = st->ilms_BnNS;  
	IMS_DATA *buffer    = st->ilms_buffer;
	IMS_DATA *rbuffer   = st->ilms_rbuffer;
	IMS_DATA *rsoft     = st->ilms_rsoft;
	IMS_DEC_STATE *prev = st->ilms_dcs;		
	IMS_DEC_STATE *curr = st->ilms_tmps;		

//	IMS_DATA ibeta = (IMS_DATA)(beta * IL_SOFT_ONE);
	IMS_DATA ibeta = (IMS_DATA)(0.5 * IL_SOFT_ONE);// 0.4;

	int dmax = (1 << (inner_data_bits-1)) - 1;

	int m  = st->m;
	int rh = st->rh;
	int nh = st->nh;
	
	int m_ldpc = m;
	int r_ldpc = m * rh;
	int n_ldpc = m * nh;

	static IMS_DATA curr_v2c[100000];



	for( i = 0; i < n_ldpc; i++ )
		soft[i] = y[i] << IL_SOFT_FPP;

//	il_min_sum_reset( prev, r_ldpc, signs, rh*n_ldpc );
	il_min_sum_reset( st );


	// check input codeword
	memset( synd, 0, r_ldpc*sizeof(synd[0]) );
	icheck_syndrome( matr, rh, nh, soft, rsoft, m_ldpc, synd );  

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
			IMS_DATA *prev_c2v_abs = buffer;
			int stateOffset = j*m_ldpc;

			for( k = 0; k < m_ldpc; k++ )
			{
				curr[k].min1 = 10000;
				curr[k].min2 = 10000;
				curr[k].pos  = 0;
				curr[k].sign = 0;
			}

			for( k = 0; k < nh; k++ )
			{
				int memOffset = j * n_ldpc + k * m_ldpc;
				IMS_DATA *curr_soft    = &soft[k * m_ldpc];
				int circ = matr[j][k];

				if( circ != SKIP )
				{
					rotate( &soft[k*m_ldpc], rsoft, circ, sizeof(soft[0]), m_ldpc ); 

					for( n = 0; n < m_ldpc; n++ )
						prev_c2v_abs[n] = prev[stateOffset+n].pos == k ? prev[stateOffset+n].min2 : prev[stateOffset+n].min1;

					for( n = 0; n < m_ldpc; n++ )
					{
						int	prev_c2v_sgn  = signs[memOffset+n] ^ prev[stateOffset+n].sign;
						IMS_DATA	prev_c2v_val  = prev_c2v_sgn ? -prev_c2v_abs[n] : prev_c2v_abs[n];
						IMS_DATA curr_v2c_val  = rsoft[n] - prev_c2v_val;

						int	curr_v2c_sgn = curr_v2c_val < 0;
						IMS_DATA curr_v2c_abs = (curr_v2c_val < 0.0 ? -curr_v2c_val : curr_v2c_val); //shifting 
						curr_v2c_abs -= ibeta;

						//curr_v2c_sgn = curr_v2c_abs < 0 ? 0 : curr_v2c_sgn;
						curr_v2c_abs = curr_v2c_abs < 0 ? 0 : curr_v2c_abs;
						curr_v2c_abs = curr_v2c_abs > dmax ? dmax : curr_v2c_abs;

#ifdef Z11111111
						curr_v2c[k * m_ldpc + n] = curr_v2c_val;
#else
						curr_soft[n] = curr_v2c_val;
#endif
						signs[memOffset+n] = curr_v2c_sgn;

						iprocess_check_node( &curr[n], curr_v2c_sgn, curr_v2c_abs, k );

					}  
				}
			}

			for( k = 0; k < m_ldpc; k++ )
				prev[stateOffset+k] = curr[k];

			for( k = 0; k < nh; k++ )
			{
				IMS_DATA *curr_c2v_val = buffer;
				IMS_DATA *curr_soft    = &soft[k * m_ldpc];
				int memOffset = j * n_ldpc + k * m_ldpc;
				int circ = matr[j][k];

				if( circ != SKIP )
				{
					for( n = 0; n < m_ldpc; n++ )
					{
						IMS_DATA curr_c2v_abs = prev[stateOffset+n].pos == k ? prev[stateOffset+n].min2 : prev[stateOffset+n].min1;

						curr_c2v_val[n] = (signs[memOffset+n] ^ prev[stateOffset+n].sign) ? -curr_c2v_abs : curr_c2v_abs;
					}

#ifdef Z11111111
					for( n = 0; n < m_ldpc; n++ )
						buffer[n] = curr_v2c[k * m_ldpc + n] + curr_c2v_val[n];

					rotate( buffer, rbuffer, m_ldpc - circ, sizeof(curr_c2v_val[0]), m_ldpc ); 

					{
						double omega = 1.0/16.0;
						double omega1 = 1 + omega;

						for( n = 0; n < m_ldpc; n++ )
						{
							double x = omega1 * rbuffer[n] - omega * curr_soft[n];
							double ax = x < 0.0 ? -x : x;
							int s  = x < 0.0 ? 1 : 0;
							int ix = (int)(ax + 0.5);
							curr_soft[n] = s ? -ix : ix;
						}
					}
#else
					for( n = 0; n < m_ldpc; n++ )
						buffer[n] = curr_soft[n] + curr_c2v_val[n];

					rotate( buffer, rbuffer, m_ldpc - circ, sizeof(curr_c2v_val[0]), m_ldpc ); 

					for( n = 0; n < m_ldpc; n++ )
						curr_soft[n] = rbuffer[n];
#endif
				}
			}

		} 

		// check syndrome
		icheck_syndrome( matr, rh, nh, soft, rsoft, m_ldpc, synd );  

		for( parity = 0, i = 0; i < r_ldpc; i++ )
			parity |= synd[i];


		if( parity == 0 )
			break; 

	} // stage

	
	for( k = 0; k < n_ldpc; k++ )
		decword[k] = soft[k] < 0;

	return parity ? -iter : iter+1; 
}
#else //CHINA_VERSION

int il_min_sum_decod_qc_lm( DEC_STATE* st, int y[], int decword[], int maxsteps, double alpha, double beta, int inner_data_bits )
{
	int   i, j, k, n;
	int iter;
	int parity;
	int **matr         = st->hd;
	int *synd          = st->syndr; 
	IMS_DATA *soft     = st->ilms_soft;
	int *signs       = st->ilms_BnNS;  
	IMS_DATA *buffer    = st->ilms_buffer;
	IMS_DATA *rbuffer   = st->ilms_rbuffer;
	IMS_DATA *rsoft     = st->ilms_rsoft;
	IMS_DEC_STATE *prev = st->ilms_dcs;		
	IMS_DEC_STATE *curr = st->ilms_tmps;		

//	IMS_DATA ibeta = (IMS_DATA)(beta * IL_SOFT_ONE);
	IMS_DATA ibeta = (IMS_DATA)(0.5 * IL_SOFT_ONE);// 0.4;

	int dmax = (1 << (inner_data_bits-1)) - 1;

	int m  = st->m;
	int rh = st->rh;
	int nh = st->nh;
	
	int m_ldpc = m;
	int r_ldpc = m * rh;
	int n_ldpc = m * nh;

	for( i = 0; i < n_ldpc; i++ )
		soft[i] = y[i] << IL_SOFT_FPP;

//	il_min_sum_reset( prev, r_ldpc, signs, rh*n_ldpc );
	il_min_sum_reset( st );

#ifdef KEEP_STATISTIC
	for( i = 0; i < n_ldpc; i++ ) st->prev_soft[i] = soft[i];
	for( i = 0; i < n_ldpc; i++ ) st->sign_counter[i] = 0;
	for( i = 0; i < n_ldpc; i++ ) st->min_abs_llr[i]  = 10000;
#endif

#if 01


	// check input codeword
	memset( synd, 0, r_ldpc*sizeof(synd[0]) );
	icheck_syndrome( matr, rh, nh, soft, rsoft, m_ldpc, synd );  

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
			IMS_DATA *prev_c2v_abs = buffer;
			int stateOffset = j*m_ldpc;

			for( k = 0; k < m_ldpc; k++ )
			{
				curr[k].min1 = 10000;
				curr[k].min2 = 10000;
				curr[k].pos  = 0;
				curr[k].sign = 0;
			}


#ifdef MY_VERSION
			for( k = 0; k < nh; k++ )
			{
				int memOffset = j * n_ldpc + k * m_ldpc;
				IMS_DATA *curr_soft    = &soft[k * m_ldpc];
				int circ = matr[j][k];

				if( circ != SKIP )
				{
					rotate( &soft[k*m_ldpc], rsoft, circ, sizeof(soft[0]), m_ldpc ); 

					for( n = 0; n < m_ldpc; n++ )
						prev_c2v_abs[n] = prev[stateOffset+n].pos == k ? prev[stateOffset+n].min2 : prev[stateOffset+n].min1;

					for( n = 0; n < m_ldpc; n++ )
					{
						int	prev_c2v_sgn  = signs[memOffset+n] ^ prev[stateOffset+n].sign;
						IMS_DATA	prev_c2v_val  = prev_c2v_sgn ? -prev_c2v_abs[n] : prev_c2v_abs[n];
						IMS_DATA curr_v2c_val  = rsoft[n] - prev_c2v_val;
#if 0
						int	curr_v2c_sgn = curr_v2c_val < 0;
						MS_DATA curr_v2c_abs = (curr_v2c_val < 0.0 ? -curr_v2c_val : curr_v2c_val)*alpha; //scaling 
#else
						int	curr_v2c_sgn = curr_v2c_val < 0;
						IMS_DATA curr_v2c_abs = (curr_v2c_val < 0.0 ? -curr_v2c_val : curr_v2c_val); //shifting 
						curr_v2c_abs -= ibeta;

						//curr_v2c_sgn = curr_v2c_abs < 0 ? 0 : curr_v2c_sgn;
						curr_v2c_abs = curr_v2c_abs < 0 ? 0 : curr_v2c_abs;
						curr_v2c_abs = curr_v2c_abs > dmax ? dmax : curr_v2c_abs;

#endif
						/*curr_v2c*/curr_soft[n] = curr_v2c_val;
						signs[memOffset+n] = curr_v2c_sgn;

						iprocess_check_node( &curr[n], curr_v2c_sgn, curr_v2c_abs, k );

					}  
				}
			}

			for( k = 0; k < m_ldpc; k++ )
				prev[stateOffset+k] = curr[k];

			for( k = 0; k < nh; k++ )
			{
				IMS_DATA *curr_c2v_val = buffer;
				IMS_DATA *curr_soft    = &soft[k * m_ldpc];
				int memOffset = j * n_ldpc + k * m_ldpc;
				int circ = matr[j][k];

				if( circ != SKIP )
				{
					for( n = 0; n < m_ldpc; n++ )
					{
						IMS_DATA curr_c2v_abs = prev[stateOffset+n].pos == k ? prev[stateOffset+n].min2 : prev[stateOffset+n].min1;

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
				IMS_DATA *curr_soft    = &soft[k * m_ldpc];
				int circ = matr[j][k];

				if( circ != SKIP )
				{
					rotate( &soft[k*m_ldpc], rsoft, circ, sizeof(soft[0]), m_ldpc ); 

					for( n = 0; n < m_ldpc; n++ )
					{
						int	prev_c2v_sgn  = signs[memOffset+n] ^ prev[stateOffset+n].sign;
#if 0
						MS_DATA	prev_c2v_val  = iget_c2v_val( prev_c2v_sgn, k, &prev[stateOffset+n], 1.0/*alpha*/, beta );
#else
						IMS_DATA beta         = ibeta_control( &prev[stateOffset+n] );
						IMS_DATA	prev_c2v_val  = iget_c2v_val( prev_c2v_sgn, k, &prev[stateOffset+n], 1.0/*alpha*/, beta );
#endif					
						IMS_DATA curr_v2c_val  = rsoft[n] - prev_c2v_val;
						int	curr_v2c_sgn = curr_v2c_val < 0;
						IMS_DATA curr_v2c_abs = (curr_v2c_val < 0.0 ? -curr_v2c_val : curr_v2c_val); 

						curr_soft[n]       = curr_v2c_val;
						signs[memOffset+n] = curr_v2c_sgn;

						iprocess_check_node( &curr[n], curr_v2c_sgn, curr_v2c_abs, k );
					}  
				}
			}

			for( k = 0; k < m_ldpc; k++ )
				prev[stateOffset+k] = curr[k];

			for( k = 0; k < nh; k++ )
			{
				IMS_DATA *curr_c2v_val = buffer;
				IMS_DATA *curr_soft    = &soft[k * m_ldpc];
				int memOffset = j * n_ldpc + k * m_ldpc;
				int circ = matr[j][k];

				if( circ != SKIP )
				{
					for( n = 0; n < m_ldpc; n++ )
					{
						int curr_c2v_sgn  = signs[memOffset+n] ^ prev[stateOffset+n].sign;
#if 0
						curr_c2v_val[n]  = get_c2v_val( curr_c2v_sgn, k, &prev[stateOffset+n], 1.0/*alpha*/, beta );
#else
						IMS_DATA beta     = ibeta_control( &prev[stateOffset+n] );
						curr_c2v_val[n]  = iget_c2v_val( curr_c2v_sgn, k, &prev[stateOffset+n], 1.0/*alpha*/, ibeta );
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
		icheck_syndrome( matr, rh, nh, soft, rsoft, m_ldpc, synd );  

		for( parity = 0, i = 0; i < r_ldpc; i++ )
			parity |= synd[i];


		if( parity == 0 )
			break; 

	} // stage
#else
	memset( signs, 0, rh*n_ldpc*sizeof(signs[0]) );

	// check input codeword
	memset( synd, 0, r_ldpc*sizeof(synd[0]) );
	icheck_syndrome( matr, rh, nh, soft, rsoft, m_ldpc, synd );  

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
				curr[k].min1 = 10000;
				curr[k].min2 = 10000;
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
						int	prev_c2v_sgn = signs[memOffset+n] ^ prev[stateOffset+n].sign;
						IMS_DATA	prev_c2v_val = prev_c2v_sgn ? -buffer[n] : buffer[n];
						IMS_DATA curr_v2c_val = rsoft[n] - prev_c2v_val;
#if 0
						int	curr_v2c_sgn = curr_v2c_val < 0;
						MS_DATA curr_v2c_abs = (curr_v2c_val < 0.0 ? -curr_v2c_val : curr_v2c_val) * alpha; //scaling 
#else
						int	curr_v2c_sgn = curr_v2c_val < 0;
						IMS_DATA curr_v2c_abs  = (curr_v2c_val < 0.0 ? -curr_v2c_val : curr_v2c_val); //shifting 
						curr_v2c_abs -= ibeta;
						if( curr_v2c_abs < 0 )
						{
							curr_v2c_abs = 0;
							curr_v2c_sgn = 0;
						}
#endif
						signs[memOffset+n] = curr_v2c_sgn;

						iprocess_check_node( &curr[n], curr_v2c_sgn, curr_v2c_abs, k );
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
						IMS_DATA tmp = prev[stateOffset+n].pos == k ? prev[stateOffset+n].min2 : prev[stateOffset+n].min1;

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
		icheck_syndrome( matr, rh, nh, soft, rsoft, m_ldpc, synd );  

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
#endif //CHINA_VERSION



#define MAX_STATE 4
DEC_STATE *state[MAX_STATE];
int state_num;
void open_ext_il_minsum( int irate, int M, int num )
{
	int i, k, n;
	CODE_CFG code;
	int mh, nh;
	int *matr_ptr;

	state_num = num+1;

	for( i = 0; i < state_num; i++ )
	{
		if( i == 0 )
			code = select_code( irate, M );
		else
			code = select_code_ext( irate, M, i );

		mh = code.nrow;
		nh = code.ncol;
		
		state[i] = decod_open( 1, mh, nh, M ); 

		matr_ptr = code.matrix;
		for( k = 0; k < mh; k++ )
			for( n = 0; n < nh; n++ )
				state[i]->hd[k][n] = *matr_ptr++;

	}
}


int ext_il_min_sum( int *dec_input, int *dec_output, int n_iter, double alpha, double beta, int inner_data_bits )
{
	int stn;
	int iter;
	int **matr_org = NULL;
	IMS_DATA *soft_org = NULL;
	int codelen_org = 0;
	int rh_org = 0;
	int nh_org = 0;

	for( stn = 0; stn < state_num; stn++ )
	{
		DEC_STATE *dec_state = state[stn];

#if 0
		iter = il_min_sum_decod_qc_lm( dec_state, dec_input, dec_output, n_iter, alpha, beta, inner_data_bits );
#else
		int rh = dec_state->rh;
		int nh = dec_state->nh;
		int m  = dec_state->m;
		int *y = dec_input;
		int **matr = dec_state->hd;
		int *syndr = dec_state->syndr;
		IMS_DATA *soft = dec_state->ilms_soft;
		IMS_DATA *rsoft = dec_state->ilms_rsoft;
		int r_ldpc = rh * m;
		int codelen = nh * m;

		int parity;

		if( stn == 0 )
		{
			rh_org = rh;
			nh_org = nh;
			matr_org = matr;
			soft_org = soft;
			codelen_org = codelen;
		}

		for( int i = 0; i < codelen; i++ )
			soft[i] = y[i] << IL_SOFT_FPP; 

		for( int i = 0; i < codelen - codelen_org; i++ )
			soft[codelen_org + i] = 0;

		il_min_sum_reset( dec_state );
		
		for( iter = 0; iter < n_iter; iter++ )
		{
			il_min_sum_iterate( dec_state, inner_data_bits );

			memset( dec_state->syndr, 0, r_ldpc * sizeof( syndr[0]) );
			icheck_syndrome( matr_org, rh_org, nh_org, soft, rsoft, m, syndr );  
//			icheck_syndrome( matr, rh, nh, soft, rsoft, m, syndr );  
		
			parity = 0;
			for( int i = 0; i < r_ldpc; i++ )
			{
				if( syndr[i] )
				{
					parity = 1;
					break;
				}
			}

			if( parity == 0 )
				break;
		}

		if( parity == 0 )
		{
			if( stn > 0 && stn < state_num )
			{
				for( int i = 0; i < codelen_org; i++ )
					soft_org[i] = soft[i];
				stn = stn;
			}

			break;
		}
#endif
	}

	for( int k = 0; k < codelen_org; k++ )
		dec_output[k] = soft_org[k] < 0;

	if( iter == n_iter )
		iter = -iter;
	else
		iter += 1;


	return iter;
}


void close_ext_il_minsum( void )
{
	int i;

	decod_close( state[0] );

	for( i = 1; i < state_num; i++ )
		decod_close( state[i] );
}

int ext_il_minsum( DEC_STATE* st, int soft[], int decword[], int maxiter, double alpha, double beta,  int inner_data_bits )
{
	return maxiter;
}
