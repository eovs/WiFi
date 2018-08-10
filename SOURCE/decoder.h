
#ifndef _DECODER_H_
#define _DECODER_H_

#include <vector>
class itmo_ldpc_dec_engine_t {
  public:
    enum class ret_status { ERROR = 0, OK = 1, ET = 2 };
    itmo_ldpc_dec_engine_t();
    void                     init(const std::vector<std::vector<int>>&check_matrix, int z);
    void                     reset();
    void                     push(const std::vector<int>&in);
    void/*ret_status */              iterate();
    const std::vector<bool>&pull();
    bool  calc_parity_check();
  private:
};
#endif _DECODER_H_
