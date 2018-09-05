#include <vector>

#include "itmo_ldpc_dec_engine.h"
#include "decoders.h"

using namespace std;
itmo_ldpc_dec_engine_t::itmo_ldpc_dec_engine_t()
{
	dec_state = NULL;
}


itmo_ldpc_dec_engine_t::~itmo_ldpc_dec_engine_t()
{
	decod_close(dec_state);
	dec_state = NULL;
}

void  itmo_ldpc_dec_engine_t::init(const std::vector<std::vector<int>>&check_matrix, int z)
{
	int nrow = (int)check_matrix.size();
	int ncol = (int)check_matrix[0].size();
	
	codewordLen = ncol * z;
	decword = new vector<bool>(codewordLen);

	dec_state = decod_open( ILMS_DEC, nrow, ncol, z );
	
	for( int i = 0; i < nrow; i++ )
		for( int j = 0; j < ncol; j++ )
			dec_state->hd[i][j] = check_matrix[i][j];

}

void  itmo_ldpc_dec_engine_t::reset()
{
	il_min_sum_reset( dec_state );

}

void  itmo_ldpc_dec_engine_t::push(const std::vector<int>&in)
{
	size_t codelen = in.size();
	for( int i = 0; i < codelen; i++ )
		dec_state->ilms_y[i] = in[i]; 

	for( int i = 0; i < codelen; i++ )
		dec_state->ilms_soft[i] = dec_state->ilms_y[i] << IL_SOFT_FPP; 
}

int  itmo_ldpc_dec_engine_t::iterate()
{
	int inner_data_bits = 7;
	int res = il_min_sum_iterate( dec_state, inner_data_bits );
	return res;
}

const std::vector<bool>& itmo_ldpc_dec_engine_t::pull()
{
	for( int i = 0; i < codewordLen; i++ )
		(*decword)[i] = dec_state->ilms_soft[i] < 0;
	return *decword;
}

bool  itmo_ldpc_dec_engine_t::calc_parity_check()
{
	int parity = 0;
	int r_ldpc = dec_state->m * dec_state->rh;

	memset( dec_state->syndr, 0, r_ldpc*sizeof(dec_state->syndr[0]) );

	icheck_syndrome( dec_state->hd, dec_state->rh, dec_state->nh, dec_state->ilms_soft, dec_state->ilms_rsoft, dec_state->m, dec_state->syndr );  
	for( int i = 0; i < r_ldpc; i++ )
	{
		if( dec_state->syndr[i] )
		{
			parity = 1;
			break;
		}
	}
	return parity ? true : false;
}