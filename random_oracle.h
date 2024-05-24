/*
 *  SPDX-License-Identifier: MIT
 */

#ifndef FAEST_RANDOM_ORACLE_H
#define FAEST_RANDOM_ORACLE_H

#include "hash_shake.h"
#include "macros.h"
#include "stdlib.h"

FAEST_BEGIN_C_DECL

// implementation of H_0

typedef hash_context H0_context_t;

void H0_init(H0_context_t* H0_ctx, unsigned int security_param);
void H0_update(H0_context_t* H0_ctx, const uint8_t* src, size_t len);
void H0_final(H0_context_t* H0_ctx, uint8_t* seed, size_t seed_len, uint8_t* commitment,
              size_t commitment_len);

typedef hash_context_x4 H0_context_x4_t;

void H0_x4_init(H0_context_x4_t* H0_ctx, unsigned int security_param);
void H0_x4_update(H0_context_x4_t* H0_ctx, const uint8_t* src0, const uint8_t* src1,
                  const uint8_t* src2, const uint8_t* src3, size_t len);
void H0_x4_final(H0_context_x4_t* ctx, uint8_t* seed0, uint8_t* seed1, uint8_t* seed2,
                 uint8_t* seed3, size_t seed_len, uint8_t* commitment0, uint8_t* commitment1,
                 uint8_t* commitment2, uint8_t* commitment3, size_t commitment_len);

// implementation of H_1

typedef hash_context H1_context_t;

void H1_init(H1_context_t* H1_ctx, unsigned int security_param);
void H1_update(H1_context_t* H1_ctx, const uint8_t* src, size_t len);
void H1_final(H1_context_t* H1_ctx, uint8_t* digest, size_t len);

// implementation of H_2

typedef hash_context H2_context_t;

void H2_init(H2_context_t* ctx, unsigned int security_param);
void H2_update(H2_context_t* ctx, const uint8_t* src, size_t len);
void H2_final(H2_context_t* ctx, uint8_t* digest, size_t len);

// implementation for H_3

typedef hash_context H3_context_t;

void H3_init(H3_context_t* ctx, unsigned int security_param);
void H3_update(H3_context_t* ctx, const uint8_t* src, size_t len);
void H3_final(H3_context_t* ctx, uint8_t* digest, size_t len, uint8_t* iv);


// implementation for H_mat

typedef hash_context H_mat_context_t;

void H_mat_init(H_mat_context_t* ctx, unsigned int security_param);
void H_mat_update(H_mat_context_t* ctx, const uint8_t* src, size_t len);
void H_mat_final(H_mat_context_t* ctx, uint8_t* digest, size_t len);

uint8_t* generate_H_mat(int n,int m,const uint8_t* seed,int lambda);

int getH(int i,int j,int n,int m,const uint8_t *buffer);

// implementation for H_e

typedef hash_context H_e_context_t;

void H_e_init(H_e_context_t* ctx, unsigned int security_param);
void H_e_update(H_e_context_t* ctx, const uint8_t* src, size_t len);
void H_e_final(H_e_context_t* ctx, uint8_t* digest, size_t len);

void generate_e(uint8_t *buffer,int m,int w,int d,const uint8_t* seed,int lambda);


typedef hash_context H_c_context_t;

void H_c_init(H_c_context_t* ctx, unsigned int security_param);
void H_c_update(H_c_context_t* ctx, const uint8_t* src, size_t len);
void H_c_final(H_c_context_t* ctx, uint8_t* digest, size_t len);


FAEST_END_C_DECL

#endif
