#ifndef SANDWICHER_SANDWICH_128_H
#define SANDWICHER_SANDWICH_128_H

#include "stdio.h"
#include "stdlib.h"
#include "string.h"
#include "fields.h"
#include "sandwich.h"

#define BITS 64

#define bf_t bf128_t
#define sf_t bf64_t

typedef struct sandwich_128_param_t {

    bf_t alpha;

    sf_t matMDS[3][3];
    sf_t top_consts[4][2], bot_consts[4][2];
    sf_t top_P[4][BITS], bot_P[4][BITS];
    sf_t inv_P[BITS];
    sf_t inv_const; 

} sandwich_128_param_t;


bf_t bf_public_128(int id, int bit,bf_t delta);

bf_t bf64_to_bf128(sandwich_128_param_t* para,sf_t in);

bf_t bf128_convert_combine(sandwich_128_param_t* para,const bf_t in[BITS]);

void init_sandwich_128(sandwich_128_param_t* para);

void sandwich_128(sandwich_128_param_t* para,sf_t k0, sf_t k1, sf_t out[2],sf_t witness[9],sf_t mul_inputs[13]);

void sandwich_bitlevel_128(sandwich_128_param_t* para,int id, const bf_t out[2][BITS],const bf_t witness[9][BITS],const bf_t delta,
                        bf_t mul_inputs[13][BITS],bf_t newk[2][BITS]);



#undef BITS 

#undef bf_t
#undef sf_t

#endif