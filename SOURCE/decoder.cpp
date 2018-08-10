#include "stdafx.h"

#include <vector>

#include "decoder.h"

using namespace std;
itmo_ldpc_dec_engine_t::itmo_ldpc_dec_engine_t()
{
}

void  itmo_ldpc_dec_engine_t::init(const std::vector<std::vector<int>>&check_matrix, int z)
{
}

void  itmo_ldpc_dec_engine_t::reset()
{
}

void  itmo_ldpc_dec_engine_t::push(const std::vector<int>&in)
{
}

/*ret_status*/void  itmo_ldpc_dec_engine_t::iterate()
{
}

const std::vector<bool>& itmo_ldpc_dec_engine_t::pull()
{
	std::vector<bool> res;

	return res;
}

bool  itmo_ldpc_dec_engine_t::calc_parity_check()
{
	return true;
}