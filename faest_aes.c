


#include "faest_aes_128.c"
#undef bf_t
#undef bf_load
#undef bf_store
#undef bf_from_bit
#undef bf_zero
#undef bf_one
#undef bf_add
#undef bf_mul
#undef bf_sum_poly
#undef bf_inv
#undef zk_hash



#include "faest_aes_192.c"
#undef bf_t
#undef bf_load
#undef bf_store
#undef bf_from_bit
#undef bf_zero
#undef bf_one
#undef bf_add
#undef bf_mul
#undef bf_sum_poly
#undef bf_inv
#undef zk_hash





#include "faest_aes_256.c"
#undef bf_t
#undef bf_load
#undef bf_store
#undef bf_from_bit
#undef bf_zero
#undef bf_one
#undef bf_add
#undef bf_mul
#undef bf_sum_poly
#undef bf_inv
#undef zk_hash



void aes_prove(const uint8_t* w, const uint8_t* u, uint8_t** V, const uint8_t* in,
               const uint8_t* out, const uint8_t* chall, uint8_t* a_tilde, uint8_t* b_tilde,
               const faest_paramset_t* params) {
 
  switch (params->faest_param.lambda) {
  case 256:
      aes_prove_256(w, u, V, in, out, chall, a_tilde, b_tilde, params);
    break;
  case 192:
      aes_prove_192(w, u, V, in, out, chall, a_tilde, b_tilde, params);
    break;
  default:
      aes_prove_128(w, u, V, in, out, chall, a_tilde, b_tilde, params);
  } 
}

uint8_t* aes_verify(const uint8_t* d, uint8_t** Q, const uint8_t* chall_2, const uint8_t* chall_3,
                    const uint8_t* a_tilde, const uint8_t* in, const uint8_t* out,
                    const faest_paramset_t* params) {
  switch (params->faest_param.lambda) {
  case 256:
      return aes_verify_256(d, Q, chall_2, chall_3, a_tilde, in, out, params);
  case 192:
      return aes_verify_192(d, Q, chall_2, chall_3, a_tilde, in, out, params);
  default:
      return aes_verify_128(d, Q, chall_2, chall_3, a_tilde, in, out, params);
  }
}
