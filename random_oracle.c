/*
 *  SPDX-License-Identifier: MIT
 */

#if defined(HAVE_CONFIG_H)
#include <config.h>
#endif

#include "random_oracle.h"

static const uint8_t domain_sep_H0 = 0;
static const uint8_t domain_sep_H1 = 1;
static const uint8_t domain_sep_H2 = 2;
static const uint8_t domain_sep_H3 = 3;

static const uint8_t domain_sep_H_mat = 4;
static const uint8_t domain_sep_H_e = 5;

static const uint8_t domain_sep_H_c = 6;

// H_0
void H0_init(H0_context_t* ctx, unsigned int security_param) {
  hash_init(ctx, security_param == 128 ? 128 : 256);
}

void H0_update(H0_context_t* ctx, const uint8_t* src, size_t len) {
  hash_update(ctx, src, len);
}

void H0_final(H0_context_t* ctx, uint8_t* seed, size_t seed_len, uint8_t* commitment,
              size_t commitment_len) {
  hash_update(ctx, &domain_sep_H0, sizeof(domain_sep_H0));
  hash_final(ctx);
  hash_squeeze(ctx, seed, seed_len);
  hash_squeeze(ctx, commitment, commitment_len);
  hash_clear(ctx);
}

void H0_x4_init(H0_context_x4_t* ctx, unsigned int security_param) {
  hash_init_x4(ctx, security_param == 128 ? 128 : 256);
}

void H0_x4_update(H0_context_x4_t* ctx, const uint8_t* src0, const uint8_t* src1,
                  const uint8_t* src2, const uint8_t* src3, size_t len) {
  hash_update_x4_4(ctx, src0, src1, src2, src3, len);
}

void H0_x4_final(H0_context_x4_t* ctx, uint8_t* seed0, uint8_t* seed1, uint8_t* seed2,
                 uint8_t* seed3, size_t seed_len, uint8_t* commitment0, uint8_t* commitment1,
                 uint8_t* commitment2, uint8_t* commitment3, size_t commitment_len) {
  hash_update_x4_1(ctx, &domain_sep_H0, sizeof(domain_sep_H0));
  hash_final_x4(ctx);
  hash_squeeze_x4_4(ctx, seed0, seed1, seed2, seed3, seed_len);
  hash_squeeze_x4_4(ctx, commitment0, commitment1, commitment2, commitment3, commitment_len);
  hash_clear_x4(ctx);
}

// H_1
void H1_init(H1_context_t* ctx, unsigned int security_param) {
  hash_init(ctx, security_param == 128 ? 128 : 256);
}

void H1_update(H1_context_t* ctx, const uint8_t* src, size_t len) {
  hash_update(ctx, src, len);
}

void H1_final(H1_context_t* ctx, uint8_t* digest, size_t len) {
  hash_update(ctx, &domain_sep_H1, sizeof(domain_sep_H1));
  hash_final(ctx);
  hash_squeeze(ctx, digest, len);
  hash_clear(ctx);
}

// H_2
void H2_init(H2_context_t* ctx, unsigned int security_param) {
  hash_init(ctx, security_param == 128 ? 128 : 256);
}

void H2_update(H2_context_t* ctx, const uint8_t* src, size_t len) {
  hash_update(ctx, src, len);
}

void H2_final(H2_context_t* ctx, uint8_t* digest, size_t len) {
  hash_update(ctx, &domain_sep_H2, sizeof(domain_sep_H2));
  hash_final(ctx);
  hash_squeeze(ctx, digest, len);
  hash_clear(ctx);
}

// H_3
void H3_init(H3_context_t* ctx, unsigned int security_param) {
  hash_init(ctx, security_param == 128 ? 128 : 256);
}

void H3_update(H3_context_t* ctx, const uint8_t* src, size_t len) {
  hash_update(ctx, src, len);
}

void H3_final(H3_context_t* ctx, uint8_t* digest, size_t len, uint8_t* iv) {
  hash_update(ctx, &domain_sep_H3, sizeof(domain_sep_H3));
  hash_final(ctx);
  hash_squeeze(ctx, digest, len);
  hash_squeeze(ctx, iv, 16);
  hash_clear(ctx);
}



// H_mat
void H_mat_init(H_mat_context_t* ctx, unsigned int security_param) {
  hash_init(ctx, security_param == 128 ? 128 : 256);
}

void H_mat_update(H_mat_context_t* ctx, const uint8_t* src, size_t len) {
  hash_update(ctx, src, len);
}

void H_mat_final(H_mat_context_t* ctx, uint8_t* digest, size_t len) {
  hash_update(ctx, &domain_sep_H_mat, sizeof(domain_sep_H_mat));
  hash_final(ctx);
  hash_squeeze(ctx, digest, len);
  hash_clear(ctx);
}
uint8_t* generate_H_mat(int n,int m,const uint8_t* seed,int lambda){
  int size=(n*(m-n)+7)/8;
  uint8_t* buffer = (uint8_t*)malloc(size);
  H_mat_context_t ctx;
  H_mat_init(&ctx, lambda);
  H_mat_update(&ctx, seed, lambda/8);
  H_mat_final(&ctx, buffer, size);
  return buffer;
}

int getH(int i,int j,int n,int m,const uint8_t *buffer){
    if(i==j){
      return 1;
    }else{
      if(j<n)
        return 0;
      else{
        int index = i*(m-n)+j-n;
        return buffer[index/8]>>(index%8)&1;
      }
    }
}

// H_e
void H_e_init(H_e_context_t* ctx, unsigned int security_param) {
  hash_init(ctx, security_param == 128 ? 128 : 256);
}

void H_e_update(H_e_context_t* ctx, const uint8_t* src, size_t len) {
  hash_update(ctx, src, len);
}

void H_e_final(H_e_context_t* ctx, uint8_t* digest, size_t len) {
  hash_update(ctx, &domain_sep_H_e, sizeof(domain_sep_H_e));
  hash_final(ctx);
  hash_squeeze(ctx, digest, len);
  hash_clear(ctx);
}


void generate_e(uint8_t *buffer,int m,int w,int d,const uint8_t* seed,int lambda){
  H_e_context_t ctx;
  H_e_init(&ctx, lambda);
  H_e_update(&ctx, seed, lambda/8);
  uint64_t *tmp = malloc(sizeof(uint64_t)*w);
  H_e_final(&ctx, (uint8_t*)tmp, sizeof(uint64_t)*w);
  memset(buffer,0,m);
  for(int i=0;i<m;i+=d){
    int offset = tmp[i/d]%d;
    buffer[i+offset] = 1;
  }
  free(tmp);
}


void H_c_init(H2_context_t* ctx, unsigned int security_param) {
  hash_init(ctx, security_param == 128 ? 128 : 256);
}

void H_c_update(H2_context_t* ctx, const uint8_t* src, size_t len) {
  hash_update(ctx, src, len);
}

void H_c_final(H2_context_t* ctx, uint8_t* digest, size_t len) {
  hash_update(ctx, &domain_sep_H_c, sizeof(domain_sep_H_c));
  hash_final(ctx);
  hash_squeeze(ctx, digest, len);
  hash_clear(ctx);
}
