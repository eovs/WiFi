#ifndef LDPC_DEC_ENGINE_H
#define LDPC_DEC_ENGINE_H

#include <vector>
#include <memory>

struct box_plus_cfg_t;
class ldpc_dec_engine_t {
  public:
    enum class ret_status { ERROR = 0, OK = 1, ET = 2 };
    ldpc_dec_engine_t();
    void                     init(const std::vector<std::vector<int>> &check_matrix, int z);
    void                     reset();
    void                     push(const std::vector<int> &in);
    ret_status               iterate();
    const std::vector<bool> &pull();
    bool                     calc_parity_check();

    ldpc_dec_engine_t(ldpc_dec_engine_t &&);
    ldpc_dec_engine_t &operator=(ldpc_dec_engine_t &&);
    ~ldpc_dec_engine_t();

    std::shared_ptr<box_plus_cfg_t> cfg;

  private:
    bool is_init = false;
    class impl;
    std::unique_ptr<impl> p_impl_;
};

#endif /*LDPC_DEC_ENGINE_H*/
