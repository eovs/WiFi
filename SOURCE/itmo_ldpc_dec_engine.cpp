#include <vector>

#include "itmo_ldpc_dec_engine.h"

using namespace std;

#define SKIP       -1
#define IL_SOFT_FPP	1
#define IL_SOFT_ONE (1 << IL_SOFT_FPP)

#define my_limit( x, lower, upper ) ((x) < (lower) ? (lower) : ((x) > (upper) ? (upper) : (x)))

typedef struct  
{
	int min1; 
	int min2; 
	int pos;
	int sign;
} ILMS_NODE;


typedef struct
{
	int nh; 
	int rh;
	int m;
    int codelen;
	int synd_len;
	int maxiter;
	int **hd; 

	int *soft;
	int	*sign;
	ILMS_NODE *dcs;
	ILMS_NODE *tmps;
	int *buffer;
}ILMS_STATE;

static int** Alloc2d_int( int b, int c )
{
	int **p;
	int i;

	p = new int* [b];
	p[0] = new int [b*c];
	for( i = 1; i < b; i++ )
		p[i] = p[i-1] + c;

	return p;
}

static void free2d_int( int **p )
{
	delete [] p[0];
	delete [] p;
}

static ILMS_STATE* decoder_open( int mh, int nh, int M )
{
    ILMS_STATE* st;
    
	int N = nh * M;
	int R = mh * M;

	st = new ILMS_STATE;
	memset( st, 0, sizeof(ILMS_STATE) );
    if( !st ) 
        return NULL;

	st->nh      = nh;
	st->rh      = mh;
	st->m       = M;
    st->codelen = N;
	st->synd_len = R;
 

	st->hd = Alloc2d_int( mh, nh );
	if( st->hd==NULL)
		return NULL;

	st->soft    = new int[N];
	st->sign    = new int[mh*N];
	st->dcs     = new ILMS_NODE[R];
	st->tmps    = new ILMS_NODE[M];
	st->buffer  = new int[M];

    return st;
}



static int icheck_syndrome( ILMS_STATE *state )
{
	int parity = 0;
	int r_ldpc = state->rh * state->m;
	int m_ldpc = state->m;
	int *synd = state->buffer;

	for( int j = 0; j < state->rh; j++ )
	{
		for( int n = 0; n < m_ldpc; n++ )
			synd[n] = 0;

		for( int k = 0; k < state->nh; k++ )
		{
			int circ = state->hd[j][k];

			if( circ != SKIP )
			{
				for( int n = 0; n < m_ldpc; n++ )
					synd[n] ^= state->soft[k*m_ldpc + (circ + n) % m_ldpc] < 0;
			}
		}

		for( int n = 0; n < m_ldpc; n++ )
			parity |= synd[n];
		
		if( parity )
			break;
	}

	return parity;
}

static void iprocess_check_node( ILMS_NODE *curr, int curr_v2c_sign, int abs_curr_v2c, int index )
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

static itmo_ldpc_dec_engine_t::ret_status il_min_sum_iterate( ILMS_STATE* st, int inner_data_bits )
{
	int **matr      = st->hd;
	int *soft       = st->soft;
	int *signs      = st->sign;  
	int *buffer     = st->buffer;
	ILMS_NODE *prev = st->dcs;		
	ILMS_NODE *curr = st->tmps;		

	int dmax = (1 << (inner_data_bits-1)) - 1;
	int ibeta = (int)(0.5 * IL_SOFT_ONE);
	
	int m  = st->m;
	int rh = st->rh;
	int nh = st->nh;
	
	int m_ldpc = m;
	int r_ldpc = m * rh;
	int n_ldpc = m * nh;


	// update statistic
	for( int j = 0; j < rh; j++ )
	{
		int stateOffset = j*m_ldpc;

		for( int k = 0; k < m_ldpc; k++ )
		{
			curr[k].min1 = dmax;
			curr[k].min2 = dmax;
			curr[k].pos  = 0;
			curr[k].sign = 0;
		}


		for( int k = 0; k < nh; k++ )
		{
			int memOffset = j * n_ldpc + k * m_ldpc;
			int *curr_soft    = &soft[k * m_ldpc];
			int circ = matr[j][k];

			if( circ != SKIP )
			{
				for( int n = 0; n < m_ldpc; n++ )
				{
					int	prev_c2v_sgn = signs[memOffset+n] ^ prev[stateOffset+n].sign;
					int prev_c2v_abs = prev[stateOffset+n].pos == k ? prev[stateOffset+n].min2 : prev[stateOffset+n].min1;
					int curr_v2c_val = soft[k*m_ldpc + (circ + n) % m_ldpc] - (prev_c2v_sgn ? -prev_c2v_abs : prev_c2v_abs);
					int curr_v2c_abs = (curr_v2c_val < 0.0 ? -curr_v2c_val : curr_v2c_val) - ibeta;  
					
					buffer[n] = curr_v2c_val;
					signs[memOffset+n] = curr_v2c_val < 0;
					iprocess_check_node( &curr[n], signs[memOffset+n], my_limit(curr_v2c_abs, 0, dmax), k );
				}  
				for( int n = 0; n < m_ldpc; n++ )
					curr_soft[n] = buffer[n];
			}
		}

		for( int k = 0; k < m_ldpc; k++ )
			prev[stateOffset+k] = curr[k];

		for( int k = 0; k < nh; k++ )
		{
			int *curr_soft    = &soft[k * m_ldpc];
			int memOffset = j * n_ldpc + k * m_ldpc;
			int circ = matr[j][k];

			if( circ != SKIP )
			{
				for( int n = 0; n < m_ldpc; n++ )
				{
					int curr_c2v_abs = prev[stateOffset+n].pos == k ? prev[stateOffset+n].min2 : prev[stateOffset+n].min1;
					int curr_c2v_val = (signs[memOffset+n] ^ prev[stateOffset+n].sign) ? -curr_c2v_abs : curr_c2v_abs;
					buffer[(n + circ) % m_ldpc] = curr_soft[n] + curr_c2v_val;
				}

				for( int n = 0; n < m_ldpc; n++ )
					curr_soft[n] = buffer[n];

			}
		}
	}

	// check syndrome
	if( icheck_syndrome( st ) == 0 )
		return itmo_ldpc_dec_engine_t::ret_status::ET; 
	else
		return itmo_ldpc_dec_engine_t::ret_status::OK;
}


itmo_ldpc_dec_engine_t::itmo_ldpc_dec_engine_t()
{
	state = NULL;
}


itmo_ldpc_dec_engine_t::~itmo_ldpc_dec_engine_t()
{
	ILMS_STATE *st = (ILMS_STATE*)state;

	if(st->hd)		{ free2d_int( st->hd );          st->hd       = NULL; }

	delete [] st->soft;
	delete [] st->sign;
	delete [] st->dcs;
	delete [] st->tmps;
	delete [] st->buffer;

    delete st;

	state = NULL;
	is_init = false;
	delete decword;
}

void  itmo_ldpc_dec_engine_t::init(const std::vector<std::vector<int>>&check_matrix, int z)
{
	ILMS_STATE *dec_state = (ILMS_STATE*)state;
	int nrow = (int)check_matrix.size();
	int ncol = (int)check_matrix[0].size();
	
	codewordLen = ncol * z;
	decword = new vector<bool>(codewordLen);

	dec_state = decoder_open( nrow, ncol, z );
	state = (void*)dec_state;
	
	for( int i = 0; i < nrow; i++ )
		for( int j = 0; j < ncol; j++ )
			dec_state->hd[i][j] = check_matrix[i][j];
	is_init = true;
}

void  itmo_ldpc_dec_engine_t::reset()
{
	ILMS_STATE *st = (ILMS_STATE*)state;

	int *signs      = st->sign;  
	ILMS_NODE *prev = st->dcs;		

	int m  = st->m;
	int rh = st->rh;
	int nh = st->nh;
	
	int r_ldpc = m * rh;
	int n_ldpc = m * nh;

	for( int i = 0; i < r_ldpc; i++ )
	{
		prev[i].min1 = 0;
		prev[i].min2 = 0;
		prev[i].pos  = 0;
		prev[i].sign = 0;
	}

	memset( signs, 0, n_ldpc*rh*sizeof(signs[0]) );
}

void  itmo_ldpc_dec_engine_t::push(const std::vector<int>&in)
{
	for( int i = 0; i < in.size(); i++ )
		((ILMS_STATE*)state)->soft[i] = in[i] << IL_SOFT_FPP; 
}

itmo_ldpc_dec_engine_t::ret_status  itmo_ldpc_dec_engine_t::iterate()
{
	if( !is_init )
		return ret_status::ERROR;

	return il_min_sum_iterate( (ILMS_STATE*)state, llr_bits + IL_SOFT_FPP );
}

const std::vector<bool>& itmo_ldpc_dec_engine_t::pull()
{
	for( int i = 0; i < codewordLen; i++ )
		(*decword)[i] = ((ILMS_STATE*)state)->soft[i] < 0;
	return *decword;
}

bool  itmo_ldpc_dec_engine_t::calc_parity_check()
{
	return icheck_syndrome( (ILMS_STATE*) state ) ? true : false;
}