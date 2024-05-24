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

bool owf(const uint8_t* key, const uint8_t* input, uint8_t* output, int lambda) {

  const faest_paramset_t paramset =  faest_get_paramset(lambda == 128 ? FAEST_128S : (lambda == 192 ? FAEST_192S : FAEST_256S));
  const int n = paramset.faest_param.n;
  const int m = paramset.faest_param.m;
  const int w = paramset.faest_param.w;
  const int d = paramset.faest_param.d;
  //const int lambda = paramset.faest_param.lambda;
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
}//TODO

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