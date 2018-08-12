#include "stdafx.h"

#include <vector>

#include "itmo_ldpc_dec_engine.h"
#include "decoders.h"

using namespace std;
itmo_ldpc_dec_engine_t::itmo_ldpc_dec_engine_t()
{
	dec_state = NULL;
}

void  itmo_ldpc_dec_engine_t::init(const std::vector<std::vector<int>>&check_matrix, int z)
{
	int nrow = (int)check_matrix.size();
	int ncol = (int)check_matrix[0].size();
	dec_state = decod_open( IL_MS_DEC, 1, nrow, ncol, z );
	
	for( int i = 0; i < nrow; i++ )
		for( int j = 0; j < ncol; j++ )
			dec_state->hd[i][j] = check_matrix[i][j];

}

void  itmo_ldpc_dec_engine_t::reset()
{
	il_min_sum_init( dec_state->ilms_dcs, dec_state->rh*dec_state->m, dec_state->ilms_BnNS, dec_state->rh*dec_state->nh*dec_state->m );

}

void  itmo_ldpc_dec_engine_t::push(const std::vector<int>&in)
{
	size_t codelen = in.size();
	for( int i = 0; i < codelen; i++ )
		dec_state->ilms_y[i] = in[i]; 
}

int  itmo_ldpc_dec_engine_t::iterate()
{
	return 0;
}

const std::vector<bool>& itmo_ldpc_dec_engine_t::pull()
{
	return (std::vector<bool>&)dec_state->ilms_decword;
}

bool  itmo_ldpc_dec_engine_t::calc_parity_check()
{
	bool res = true;

	icheck_syndrome( dec_state->hd, dec_state->rh, dec_state->nh, dec_state->ilms_soft, dec_state->ilms_rsoft, dec_state->m, dec_state->syndr );  
	for( int i = 0; i < dec_state->m * dec_state->rh; i++ )
	{
		if( dec_state->syndr[i] )
		{
			res = false;
			break;
		}
	}
	return res;
}