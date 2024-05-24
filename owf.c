/*
 *  SPDX-License-Identifier: MIT
 */

#if defined(HAVE_CONFIG_H)
#include <config.h>
#endif

#include "owf.h"
#include "aes.h"
#include "utils.h"
#include "random_oracle.h"
#include "fields.h"
#include "time.h"
#include "sandwich.h"



bool owf_128(const uint8_t* key, const uint8_t* input, uint8_t* output) {
  const faest_paramset_t paramset =  faest_get_paramset(FAEST_128S);
  const int n = paramset.faest_param.n;
  const int output_len = (n+7)/8;
  int ret = 0;  

  memset(output, 0, output_len);

  bf64_t k0,k1;
  bf64_t witness[9];
  bf64_t mul_inputs[13];

  memcpy(&k0,key,sizeof(k0));
  memcpy(&k1,key+sizeof(k0),sizeof(k1));

  sandwich_128_param_t param;
  init_sandwich_128(&param);

  sandwich_128(&param,k0,k1,(bf64_t*)output,witness,mul_inputs);

  return ret == 0;
}


bool owf_192(const uint8_t* key, const uint8_t* input, uint8_t* output) {
  const faest_paramset_t paramset =  faest_get_paramset(FAEST_192S);
  const int n = paramset.faest_param.n;
  const int output_len = (n+7)/8;
  int ret = 0;  

  memset(output, 0, output_len);

  bf96_t k0,k1;
  bf96_t witness[9];
  bf96_t mul_inputs[13];

  memcpy(&k0,key,sizeof(k0));
  memcpy(&k1,key+sizeof(k0),sizeof(k1));

  sandwich_192_param_t param;
  init_sandwich_192(&param);

  sandwich_192(&param,k0,k1,(bf96_t*)output,witness,mul_inputs);

  return ret == 0;
}


bool owf_256(const uint8_t* key, const uint8_t* input, uint8_t* output) {
  const faest_paramset_t paramset =  faest_get_paramset(FAEST_256S);
  const int n = paramset.faest_param.n;
  const int output_len = (n+7)/8;
  int ret = 0;  

  memset(output, 0, output_len);

  bf128_t k0,k1;
  bf128_t witness[9];
  bf128_t mul_inputs[13];

  memcpy(&k0,key,sizeof(k0));
  memcpy(&k1,key+sizeof(k0),sizeof(k1));

  sandwich_256_param_t param;
  init_sandwich_256(&param);

  sandwich_256(&param,k0,k1,(bf128_t*)output,witness,mul_inputs);

  return ret == 0;
}


bool owf(const uint8_t* key, const uint8_t* input, uint8_t* output, int lambda) {
  if (lambda == 128) {
    return owf_128(key, input, output);
  } else if (lambda == 192) {
    return owf_192(key, input, output);
  } else if (lambda == 256) {
    return owf_256(key, input, output);
  } else {
    return false;
  }
}

bool faest_128s_owf(const uint8_t* key, const uint8_t* input, uint8_t* output){
  return owf(key, input, output, 128);
}
bool faest_128f_owf(const uint8_t* key, const uint8_t* input, uint8_t* output){
  return owf(key, input, output, 128);
}

bool faest_192s_owf(const uint8_t* key, const uint8_t* input, uint8_t* output){
  return owf(key, input, output, 192);
}
bool faest_192f_owf(const uint8_t* key, const uint8_t* input, uint8_t* output){
  return owf(key, input, output, 192);
}

bool faest_256s_owf(const uint8_t* key, const uint8_t* input, uint8_t* output){
  return owf(key, input, output, 256);
}
bool faest_256f_owf(const uint8_t* key, const uint8_t* input, uint8_t* output){
  return owf(key, input, output, 256);
}