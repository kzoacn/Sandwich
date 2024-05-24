/*
 *  SPDX-License-Identifier: MIT
 */

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include "faest.h"
#include "faest_aes.h"
#include "fields.h"
#include "vole.h"
#include "universal_hashing.h"
#include "utils.h"
#include "random_oracle.h"
#include "time.h"
#include <string.h>
#include <stdlib.h>
#include "sandwich.h"


#define BITS 64 

#define bf_t bf128_t
#define bf_load bf128_load
#define bf_store bf128_store
#define bf_from_bit bf128_from_bit
#define bf_zero bf128_zero
#define bf_one bf128_one
#define bf_add bf128_add
#define bf_mul bf128_mul
#define bf_sum_poly bf128_sum_poly
#define bf_inv bf128_inv
#define zk_hash zk_hash_128 


static bf_t* column_to_row_major_and_shrink_V_128(uint8_t** v, unsigned int ell) {
  // V is \hat \ell times \lambda matrix over F_2
  // v has \hat \ell rows, \lambda columns, storing in column-major order, new_v has \ell + \lambda
  // rows and \lambda columns storing in row-major order
  bf_t* new_v = malloc((ell + sizeof(bf_t) * 8) * sizeof(bf_t));
  for (unsigned int row = 0; row != ell + sizeof(bf_t) * 8; ++row) {
    uint8_t new_row[sizeof(bf_t)] = {0};
    for (unsigned int column = 0; column != sizeof(bf_t) * 8; ++column) {
      ptr_set_bit(new_row, ptr_get_bit(v[column], row), column);
    }
    new_v[row] = bf_load(new_row);
  }

  return new_v;
} 

static void aes_prove_128(const uint8_t* w, const uint8_t* u, uint8_t** V, const uint8_t* input,
                          const uint8_t* out, const uint8_t* chall, uint8_t* a_tilde,
                          uint8_t* b_tilde, const faest_paramset_t* params) {
  const unsigned int l    = params->faest_param.l;
  const unsigned int n = params->faest_param.n;
  const unsigned int m = params->faest_param.m;
  const unsigned int D = params->faest_param.d;
  const unsigned int lambda = params->faest_param.lambda;
  const unsigned int lambdaBytes = lambda / 8;


  //double ck1=clock();
  // Step: 1..2
  bf_t* bf_v = column_to_row_major_and_shrink_V_128(V, l);

  // Step: 3..4
  // do nothing
  // Step: 6

  // Step: 7
  const unsigned int length_a = m/D*2 + 1;
  bf_t* A0                 = malloc(sizeof(bf_t) * length_a);
  bf_t* A1                 = malloc(sizeof(bf_t) * length_a);

  //fill out A0 A1

  for(uint32_t i=0;i<length_a;i++){
    A0[i] = bf_zero();
    A1[i] = bf_zero();
  }


  
  bf64_t k0,k1;
  bf64_t witness[9];
  bf64_t mul_inputs[13];
  bf64_t _out[2];
  k0 = ((const bf64_t*)w)[0];
  k1 = ((const bf64_t*)w)[1];

  sandwich_128_param_t param;
  init_sandwich_128(&param);

  sandwich_128(&param,k0,k1,_out,witness,mul_inputs);

  bf_t bf_out[2][BITS];
  bf_t bf_witness[9][BITS];
  bf_t bf_mul_inputs[13][BITS];
  bf_t bf_newk[2][BITS];

  bf_t fake_delta=bf_one();

  for(int i=0;i<2;i++){
    for(int j=0;j<BITS;j++){
      bf_out[i][j]=bf_public_128(PROVER,bf64_get_bit(_out[i],j),fake_delta);
    }
  }

  for(int i=0;i<9;i++){
    for(int j=0;j<BITS;j++){
      bf_witness[i][j]=bf_v[i*BITS+j];
    }
  }
  

  sandwich_bitlevel_128(&param,PROVER,bf_out,bf_witness,fake_delta,bf_mul_inputs,bf_newk);

  //bf_t bf_w = bf64_to_bf128(&param,witness[3]);

  int w_pos[] = {2,3,4,6,7,8};
  int mul_0_pos[] = {0,2,4,7,9,11};
  int mul_1_pos[] = {1,3,5,8,10,12};

  for(int i=0;i<6;i++){

    int p = w_pos[i];
    int x = mul_0_pos[i];
    int y = mul_1_pos[i];

    bf_t bf_w_mac = bf128_convert_combine(&param,bf_witness[p]);

    bf_t bf_mul_0 = bf64_to_bf128(&param,mul_inputs[x]);
    bf_t bf_mul_1 = bf64_to_bf128(&param,mul_inputs[y]);

    bf_t bf_mul_0_mac = bf128_convert_combine(&param,bf_mul_inputs[x]);
    bf_t bf_mul_1_mac = bf128_convert_combine(&param,bf_mul_inputs[y]);

    A0[i]=bf_mul(bf_mul_0_mac,bf_mul_1_mac);
    A1[i]=bf_add(bf_add(bf_mul(bf_mul_0_mac,bf_mul_1),bf_mul(bf_mul_1_mac,bf_mul_0)),bf_w_mac);

  }

  bf_t bf_mul_0 = bf64_to_bf128(&param,witness[5]);
  bf_t bf_mul_1 = bf64_to_bf128(&param,mul_inputs[6]);

  bf_t bf_mul_0_mac = bf128_convert_combine(&param,bf_witness[5]);
  bf_t bf_mul_1_mac = bf128_convert_combine(&param,bf_mul_inputs[6]);

  A0[6] = bf_mul(bf_mul_0_mac,bf_mul_1_mac);
  A1[6] = bf_add(bf_add(bf_mul(bf_mul_0_mac,bf_mul_1),bf_mul(bf_mul_1_mac,bf_mul_0)),bf_zero());

  for(int i=0;i<2;i++){
    bf_t newk = bf128_convert_combine(&param,bf_newk[i]);
    bf_t oldk = bf128_convert_combine(&param,bf_witness[i]);
    A0[7+i] = bf_add(newk,oldk); 
    A1[7+i] = bf_zero();
  }

  // Step: 16..18
  A1[length_a - 1] = bf_load(u + l / 8);
  A0[length_a - 1] = bf_sum_poly(bf_v + l);
  free(bf_v);


  zk_hash(a_tilde, chall, A1, length_a - 1);
  zk_hash(b_tilde, chall, A0, length_a - 1);
  //double ck5=clock();

  free(A0);
  free(A1);
 
  // printf("Total time for circuit prove: %f ms\n",(ck5-ck1)/CLOCKS_PER_SEC*1000);
  // printf("Time for computing e: %f ms\n",(ck2-ck1)/CLOCKS_PER_SEC*1000);
  // printf("Time for computing |e_i|<=1: %f ms\n",(ck3-ck2)/CLOCKS_PER_SEC*1000);
  // printf("Time for computing 1*u=1: %f ms\n",(ck4-ck3)/CLOCKS_PER_SEC*1000);
  // printf("Time for computing universal hash: %f ms\n",(ck5-ck4)/CLOCKS_PER_SEC*1000);

}

static uint8_t* aes_verify_128(const uint8_t* d, uint8_t** Q, const uint8_t* chall_2,
                               const uint8_t* chall_3, const uint8_t* a_tilde, const uint8_t* input,
                               const uint8_t* out, const faest_paramset_t* params) {
  const unsigned int lambda      = params->faest_param.lambda;
  const unsigned int tau         = params->faest_param.tau;
  const unsigned int t0          = params->faest_param.t0;
  const unsigned int k0          = params->faest_param.k0;
  const unsigned int t1          = params->faest_param.t1;
  const unsigned int k1          = params->faest_param.k1;
  const unsigned int l           = params->faest_param.l;
  const unsigned int lambdaBytes = lambda / 8;
  const unsigned int n = params->faest_param.n;
  const unsigned int m = params->faest_param.m;
  const unsigned int D = params->faest_param.d; 


  // Step: 1
  const uint8_t* _delta = chall_3;
  // Step: 2,3
  // do nothing

  // Step: 4..10
  for (uint32_t i = 0, col = 0; i < tau; i++) {
    unsigned int depth = i < t0 ? k0 : k1;
    uint8_t fancy_d[MAX_DEPTH];
    ChalDec(chall_3, i, k0, t0, k1, t1, fancy_d);
    for (uint32_t j = 0; j < depth; j++, ++col) {
      if (fancy_d[j] == 1) {
        xor_u8_array(d, Q[col], Q[col], (l + 7) / 8);
      }
    }
  }

  // Step: 11..12
  bf_t* bf_q = column_to_row_major_and_shrink_V_128(Q, l);

  // Step: 13
  const unsigned int length_b = m/D*2  + 1; 
  bf_t* B_0                = malloc(sizeof(bf_t) * length_b);  


  for(uint32_t i=0;i<length_b;i++){
    B_0[i]=bf_zero();
  }


  //fill out B0
  
  //bf64_t k0,k1;
  //bf64_t witness[9];
  //bf64_t mul_inputs[13];
  //bf64_t _out[2];
  //k0 = ((const bf64_t*)w)[0];
  //k1 = ((const bf64_t*)w)[1];

  sandwich_128_param_t param;
  init_sandwich_128(&param);

  //sandwich_128(&param,k0,k1,_out,witness,mul_inputs);

  bf_t bf_out[2][BITS];
  bf_t bf_witness[9][BITS];
  bf_t bf_mul_inputs[13][BITS];
  bf_t bf_newk[2][BITS];
  bf64_t *_out=(bf64_t*)out;

  bf_t delta = bf_load(_delta);

  for(int i=0;i<2;i++){
    for(int j=0;j<BITS;j++){
      bf_out[i][j]=bf_public_128(VERIFIER,bf64_get_bit(_out[i],j),delta);
    }
  }

  for(int i=0;i<9;i++){
    for(int j=0;j<BITS;j++){
      bf_witness[i][j]=bf_q[i*BITS+j];
    }
  }

  sandwich_bitlevel_128(&param,VERIFIER,bf_out,bf_witness,delta,bf_mul_inputs,bf_newk);


  int w_pos[] = {2,3,4,6,7,8};
  int mul_0_pos[] = {0,2,4,7,9,11};
  int mul_1_pos[] = {1,3,5,8,10,12};

  for(int i=0;i<6;i++){

    int p = w_pos[i];
    int x = mul_0_pos[i];
    int y = mul_1_pos[i];

    bf_t bf_w_mac = bf128_convert_combine(&param,bf_witness[p]);

    bf_t bf_mul_0_mac = bf128_convert_combine(&param,bf_mul_inputs[x]);
    bf_t bf_mul_1_mac = bf128_convert_combine(&param,bf_mul_inputs[y]);

    B_0[i] = bf_add( bf_mul(bf_mul_0_mac,bf_mul_1_mac) , bf_mul(bf_w_mac,delta));
  }


  bf_t bf_mul_0_mac = bf128_convert_combine(&param,bf_witness[5]);
  bf_t bf_mul_1_mac = bf128_convert_combine(&param,bf_mul_inputs[6]);

  B_0[6] = bf_add( bf_mul(bf_mul_0_mac,bf_mul_1_mac) , bf_mul(delta,delta));


  for(int i=0;i<2;i++){
    bf_t newk = bf128_convert_combine(&param,bf_newk[i]);
    bf_t oldk = bf128_convert_combine(&param,bf_witness[i]);
    B_0[7+i] = bf_add(newk,oldk); 
  }

  // Step: 20
  B_0[length_b - 1] = bf_sum_poly(bf_q + l);
  free(bf_q);

  // Step 21
  uint8_t* q_tilde = malloc(lambdaBytes);
  zk_hash(q_tilde, chall_2, B_0, length_b - 1);
  free(B_0);

  bf_t bf_qtilde = bf_load(q_tilde);
  bf_store(q_tilde, bf_add(bf_qtilde, bf_mul(bf_load(a_tilde), delta)));

  return q_tilde;
}
