#ifndef SANDWICHER_SANDWICH_192_H
#define SANDWICHER_SANDWICH_192_H

#include "stdio.h"
#include "stdlib.h"
#include "string.h"
#include "fields.h"


#define topRound 3
#define bottomRound 3

#define BITS  64

#define PROVER 0
#define VERIFIER 1


#define bf_t bf128_t
#define sf_t bf64_t


#define bf_load bf128_load
#define bf_store bf128_store
#define bf_from_bit bf128_from_bit
#define bf_zero bf128_zero
#define bf_one bf128_one
#define bf_add bf128_add
#define bf_mul bf128_mul
#define bf_sum_poly bf128_sum_poly
#define bf_inv bf128_inv


#define sf_load bf64_load
#define sf_store bf64_store
#define sf_from_bit bf64_from_bit
#define sf_zero bf64_zero
#define sf_one bf64_one
#define sf_add bf64_add
#define sf_mul bf64_mul
#define sf_sum_poly bf64_sum_poly
#define sf_inv bf64_inv
#define sf_count bf64_count
#define sf_from_bf64 bf64_from_bf64
#define sf_modulus bf64_modulus
#define sf_flip bf64_flip
#define sf_get_bit bf64_get_bit
#define sf_and bf64_and
//binary field
//subfield



#define alpha_128_s "11001100000011011110100011110100100000000001110000001011100010100001101000111011000000110000000011010011100001000011011000000000"
#define alpha_192_s "101100110000001101000000010011100010110010111000000111110010110011100001111110111100010011110111011000010011001010000001000010010110000111001110100100001111111100101001000010101000110000000000"
#define alpha_256_s "1101111101110101100001001110000111000001010011100110010011101000100100100010010110110010101001001101101111011001011001011011100011111001101111000001010000010001100101000011101010101011110101010001001110100000110110101001011001100000100011010111110000000000"


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


#undef bf_load 
#undef bf_store 
#undef bf_from_bit 
#undef bf_zero 
#undef bf_one 
#undef bf_add 
#undef bf_mul 
#undef bf_sum_poly 
#undef bf_inv 


#undef sf_load 
#undef sf_store 
#undef sf_from_bit 
#undef sf_zero 
#undef sf_one 
#undef sf_add 
#undef sf_mul 
#undef sf_sum_poly 
#undef sf_inv 
#undef sf_count 
#undef sf_from_bf64 
#undef sf_modulus 
#undef sf_flip 
#undef sf_get_bit 
#undef sf_and 

#endif