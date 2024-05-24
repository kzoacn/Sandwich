#ifndef SANDWICHER_SANDWICH_256_H
#define SANDWICHER_SANDWICH_256_H

#include "stdio.h"
#include "stdlib.h"
#include "string.h"
#include "fields.h"
#include "sandwich.h"

#define BITS 128

#define bf_t bf256_t
#define sf_t bf128_t

typedef struct sandwich_256_param_t {

    bf_t alpha;

    bf64_t matMDS[3][3];
    sf_t top_consts[4][2], bot_consts[4][2];
    sf_t top_P[4][BITS], bot_P[4][BITS];
    sf_t inv_P[BITS];
    sf_t inv_const; 

} sandwich_256_param_t;


bf_t bf_public_256(int id, int bit,bf_t delta);

bf_t bf128_to_bf256(sandwich_256_param_t* para,sf_t in);

bf_t bf256_convert_combine(sandwich_256_param_t* para,const bf_t in[BITS]);

void init_sandwich_256(sandwich_256_param_t* para);

void sandwich_256(sandwich_256_param_t* para,sf_t k0, sf_t k1, sf_t out[2],sf_t witness[9],sf_t mul_inputs[13]);

void sandwich_bitlevel_256(sandwich_256_param_t* para,int id, const bf_t out[2][BITS],const bf_t witness[9][BITS],const bf_t delta,
                        bf_t mul_inputs[13][BITS],bf_t newk[2][BITS]);



#undef BITS 

#undef bf_t
#undef sf_t

#endif