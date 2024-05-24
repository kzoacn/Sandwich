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
 

#define bf_t bf256_t
#define bf_load bf256_load
#define bf_store bf256_store
#define bf_from_bit bf256_from_bit
#define bf_zero bf256_zero
#define bf_one bf256_one
#define bf_add bf256_add
#define bf_mul bf256_mul
#define bf_sum_poly bf256_sum_poly
#define bf_inv bf256_inv
#define zk_hash zk_hash_256 


static bf_t* column_to_row_major_and_shrink_V_256(uint8_t** v, unsigned int ell) {
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

static void aes_prove_256(const uint8_t* w, const uint8_t* u, uint8_t** V, const uint8_t* input,
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
  bf_t* bf_v = column_to_row_major_and_shrink_V_256(V, l);

  // Step: 3..4
  // do nothing
  // Step: 6

  // Step: 7
  const unsigned int length_a = m/D*2 + 1;
  bf_t* A0                 = malloc(sizeof(bf_t) * length_a);
  bf_t* A1                 = malloc(sizeof(bf_t) * length_a);
  uint8_t* e = malloc(m);
  uint8_t* compact_e = malloc((m-n)/D*(D-1));
  bf_t* bf_e = malloc(sizeof(bf_t) * m);
  uint8_t* y=malloc(n);
  uint8_t *buffer;
  uint8_t *R = malloc(m*lambdaBytes);
  uint8_t *S = malloc(m*lambdaBytes);
  uint8_t *buffer2 = malloc(m*2*lambdaBytes);
  bf_t* tmp = malloc(sizeof(bf_t) * m);

  for(uint32_t i=0;i<length_a;i++){
    A0[i] = bf_zero();
    A1[i] = bf_zero();
  }

  //generate mat H
 
  buffer=generate_H_mat(n,m,input,lambda); 
  
  //unpack out to y

  for(uint32_t i=0;i<n;i++){
    y[i] = out[i/8] >> (i%8) & 1;
  }
  
  //unpack w to compact_e

  for(uint32_t i=0;i<(m-n)/D*(D-1);i++){
    compact_e[i] = w[i/8] >> (i%8) & 1;
  }
  
  uint32_t cur = 0;
  for(uint32_t i=n;i<m;i++){
    if(i%D==D-1){
      e[i]=1;
      bf_e[i]=bf_zero();
      for(uint32_t j=1;j<D;j++){
         e[i]^=e[i-j];
         bf_e[i]=bf_add(bf_e[i],bf_e[i-j]);
      }
    }else{
      e[i]=compact_e[cur];
      bf_e[i]=bf_add(bf_v[cur],bf_zero());
      cur++;
    }
  }
 
  for(uint32_t i=0;i<n;i++){
    e[i]=y[i];
    bf_e[i]=bf_zero();

    for(uint32_t j=n;j<m;j++){
      int index = i*(m-n)+j-n;
      int h=buffer[index/8]>>(index%8)&1;
      e[i]^=(h&e[j]);
      if(h) 
        bf_e[i]=bf_add(bf_e[i],bf_e[j]);
      else
        bf_e[i]=bf_add(bf_e[i],bf_zero());
    }
 
  }  

  //double ck2=clock();
  
    H_c_context_t ctx;
    H_c_init(&ctx, lambda);
    H_c_update(&ctx, chall, lambdaBytes);
    H_c_final(&ctx, buffer2,m*2*lambdaBytes);
    memcpy(R,buffer2,m*lambdaBytes);
    memcpy(S,buffer2+m*lambdaBytes,m*lambdaBytes); 
  



  for(uint32_t i=0;i<m/D;i++){
    bf_t z1 = bf_zero();
    bf_t z2 = bf_zero();
    bf_t z3 = bf_zero();
    bf_t Mz1 = bf_zero();
    bf_t Mz2 = bf_zero();
    bf_t Mz3 = bf_zero();

    for(uint32_t j=0;j<D;j++){ 
      uint32_t index = i*D+j;
      bf_t r = bf_load(R+index*lambdaBytes);
      bf_t s = bf_load(S+index*lambdaBytes);
      bf_t rs = bf_mul(r,s);
      if(e[index]==1){
        z1 = bf_add(z1,r);
        z2 = bf_add(z2,s);
        z3 = bf_add(z3,rs);
      }else{
        z1 = bf_add(z1,bf_zero());
        z2 = bf_add(z2,bf_zero());
        z3 = bf_add(z3,bf_zero());
      }

      Mz1 = bf_add(Mz1,bf_mul(r,bf_e[index]));
      Mz2 = bf_add(Mz2,bf_mul(s,bf_e[index]));
      Mz3 = bf_add(Mz3,bf_mul(rs,bf_e[index]));
    }
    // prove z1*z2=z3


    A0[i]=bf_mul(Mz1,Mz2);
    A1[i]=bf_add(bf_add(bf_mul(Mz1,z2),bf_mul(Mz2,z1)),Mz3);
  }
 
  //double ck3=clock();

  for(uint32_t i=0;i<m/D;i++){
    bf_t s = bf_zero();
    for(uint32_t j=0;j<D;j++){
      uint32_t index=i*D+j;
      s=bf_add(s,bf_e[index]);
    }
    A0[i+m/D]=bf_zero();
    A1[i+m/D]=s;
  }

  //double ck4=clock();

  free(tmp);
  free(e);
  free(compact_e);
  free(bf_e);
  free(y);
  free(buffer);
  free(R);
  free(S);
  free(buffer2);
  
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

static uint8_t* aes_verify_256(const uint8_t* d, uint8_t** Q, const uint8_t* chall_2,
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

  uint8_t* compact_e = malloc((m-n)/D*(D-1));
  bf_t* bf_e = malloc(sizeof(bf_t) * m);
  uint8_t* y=malloc(n);
  uint8_t *buffer;
  uint8_t *R = malloc(m*lambdaBytes);
  uint8_t *S = malloc(m*lambdaBytes);
  uint8_t *buffer2 = malloc(m*2*lambdaBytes);

  

  // Step: 1
  const uint8_t* delta = chall_3;
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
  bf_t* bf_q = column_to_row_major_and_shrink_V_256(Q, l);

  // Step: 13
  const unsigned int length_b = m/D*2  + 1; 
  bf_t* B_0                = malloc(sizeof(bf_t) * length_b);  


  for(uint32_t i=0;i<length_b;i++){
    B_0[i]=bf_zero();
  }


  //generate mat H
  buffer=generate_H_mat(n,m,input,lambda);
  
  //unpack out to y

  for(uint32_t i=0;i<n;i++){
    y[i] = out[i/8] >> (i%8) & 1;
  }
  
  //unpack w to compact_e
  uint32_t cur = 0;
  for(uint32_t i=n;i<m;i++){
    if(i%D==D-1){
      bf_e[i]=bf_load(delta);
      for(uint32_t j=1;j<D;j++){
         bf_e[i]=bf_add(bf_e[i],bf_e[i-j]);
      }
    }else{
      bf_e[i]=bf_add(bf_q[cur],bf_zero());
      cur++;
    }
  }
  for(uint32_t i=0;i<n;i++){
    bf_e[i]=bf_zero();
    if(y[i])
      bf_e[i]=bf_load(delta);
    for(uint32_t j=n;j<m;j++){
      if(getH(i,j,n,m,buffer)) 
        bf_e[i]=bf_add(bf_e[i],bf_e[j]);
    }
  } 
   
    H_c_context_t ctx;
    H_c_init(&ctx, lambda);
    H_c_update(&ctx, chall_2, lambdaBytes);
    H_c_final(&ctx, buffer2,m*2*lambdaBytes);
    memcpy(R,buffer2,m*lambdaBytes);
    memcpy(S,buffer2+m*lambdaBytes,m*lambdaBytes);  
 
  for(uint32_t i=0;i<m/D;i++){
    bf_t Kz1 = bf_zero();
    bf_t Kz2 = bf_zero();
    bf_t Kz3 = bf_zero();

    for(uint32_t j=0;j<D;j++){ 
      uint32_t index = i*D+j;
      bf_t r = bf_load(R+index*lambdaBytes);
      bf_t s = bf_load(S+index*lambdaBytes); 
      Kz1 = bf_add(Kz1,bf_mul(r,bf_e[index]));
      Kz2 = bf_add(Kz2,bf_mul(s,bf_e[index]));

      Kz3 = bf_add(Kz3,bf_mul(bf_mul(r,s),bf_e[index]));
    }
    // prove z1*z2=z3

    B_0[i] = bf_add( bf_mul(Kz1,Kz2) , bf_mul(Kz3,bf_load(delta)));
  }
  for(uint32_t i=0;i<m/D;i++){
    bf_t s = bf_zero();
    for(uint32_t j=0;j<D;j++){
      uint32_t index=i*D+j;
      s=bf_add(s,bf_e[index]);
    }
    B_0[i+m/D]=bf_add( bf_mul(s,bf_load(delta)),bf_mul(bf_load(delta),bf_load(delta)) );
  }

  
  free(compact_e);
  free(bf_e);
  free(y);
  free(buffer);
  free(R);
  free(S);
  free(buffer2);
  
  // Step: 20
  B_0[length_b - 1] = bf_sum_poly(bf_q + l);
  free(bf_q);

  // Step 21
  uint8_t* q_tilde = malloc(lambdaBytes);
  zk_hash(q_tilde, chall_2, B_0, length_b - 1);
  free(B_0);

  bf_t bf_qtilde = bf_load(q_tilde);
  bf_store(q_tilde, bf_add(bf_qtilde, bf_mul(bf_load(a_tilde), bf_load(delta))));

  return q_tilde;
}
