#include "stdio.h"
#include "stdlib.h"
#include "string.h"
#include "fields.h"
#include "sandwich.h"

#define BITS  64

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
#define bf_get_bit bf128_get_bit


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


#define alpha para->alpha
#define inv_const para->inv_const
#define top_consts para->top_consts
#define bot_consts para->bot_consts
#define top_P para->top_P
#define bot_P para->bot_P
#define inv_P para->inv_P
#define matMDS para->matMDS


bf_t bf_public_128(int id, int bit,bf_t delta){
    if(bit==0)
        return bf_zero();
    if(id==PROVER)
        return bf_zero();
    else
        return delta;
}

sf_t LinP_eval_128(sf_t in, const sf_t LinP[BITS]){
    sf_t res = 0;
    for(int i=0;i<BITS;i++){ 
        if(sf_count(sf_and(LinP[i],in)) % 2)
            res = sf_flip(res,i);
    }
    return res;
} 

void LinP_eval_bitlevel_128(const bf_t in[BITS], const sf_t LinP[BITS], bf_t out[BITS]){

    for(int i=0;i<BITS;i++){
        out[i] = bf_zero();
        for(int j=0;j<BITS;j++){
            if(sf_get_bit(LinP[i],j))
                out[i] = bf_add(out[i],in[j]); // TODO : constant time
        }
    }
    
}


bf_t bf64_to_bf128(sandwich_128_param_t* para,sf_t in){
    bf_t res = bf_zero();
    bf_t a = bf_one();
    for(int i=0;i<BITS;i++){
        if(sf_get_bit(in,i))
            res = bf_add(res,a);
        a=bf_mul(a,alpha);
    }
    return res;
}

void init_sandwich_128(sandwich_128_param_t* para){

    // TODO : random constants

    alpha = bf_zero();
    for(int i=127;i>=0;i--){
        alpha = bf128_shift_left_1(alpha);
        if(alpha_128_s[i] == '1')
            alpha = bf128_add(alpha,bf128_one());
    }
    para->power_of_alpha[0] = bf_one();
    for(int i=1;i<BITS;i++)
        para->power_of_alpha[i] = bf_mul(para->power_of_alpha[i-1],alpha);
    bf_t two2  = bf128_from_bf64(2);

    para->power_of_two[0] = bf_one();
    for(int i=1;i<BITS*2;i++)
        para->power_of_two[i] = bf_mul(para->power_of_two[i-1],two2);

    matMDS[0][0] = sf_from_bf64(1);
    matMDS[0][1] = sf_from_bf64(2);
    matMDS[0][2] = sf_from_bf64(1);
    
    matMDS[1][0] = sf_from_bf64(2);
    matMDS[1][1] = sf_from_bf64(1);
    matMDS[1][2] = sf_from_bf64(1);
    
    matMDS[2][0] = sf_from_bf64(1);
    matMDS[2][1] = sf_from_bf64(1);
    matMDS[2][2] = sf_from_bf64(2);

    inv_const = sf_zero();
    for(int i=0;i<4;i++){
        for(int j=0;j<2;j++){
            top_consts[i][j] = sf_zero();
            bot_consts[i][j] = sf_zero();
        } 
    }
    
    memset(top_P,0,sizeof(top_P));
    memset(bot_P,0,sizeof(bot_P));
    memset(inv_P,0,sizeof(inv_P));

    srand(123);
    for(int i=0;i<4;i++){ 
        int num = BITS*2 / 6;

        for(int j=0;j<BITS;j++){
            top_P[i][j] = sf_zero();
            bot_P[i][j] = sf_zero();
            inv_P[j] = sf_zero();
        }

        for(int j=0;j<num;j++){
            int x = rand() % BITS;
            int y = rand() % BITS;
            int z = rand() % BITS;
            top_P[i][x] = sf_from_bf64(1);
            bot_P[i][y] = sf_from_bf64(1);
            
            if(i==0)
                inv_P[z] = sf_from_bf64(1);
        }

    }

}



void feistel_128(const sf_t state[3],const sf_t consts[2],const sf_t P[BITS], sf_t out_state[3], sf_t *witness, sf_t mul_in[2]){
	out_state[0] = state[1];
	sf_t in0 = sf_add(state[0], consts[0]);
	sf_t in1 = sf_add(state[1], consts[1]);
    sf_t out = sf_mul(in0 , in1);

    if(witness != NULL)
        *witness = out;

    if(mul_in != NULL){
        mul_in[0] = in0;
        mul_in[1] = in1;
    }

    sf_t tmp = LinP_eval_128(out, P);
	out_state[1] = sf_add(tmp , state[2]);
	out_state[2] = state[0];
}
 



void feistel_bitlevel_128(int id,bf_t state[3][BITS],const sf_t consts[2],const sf_t P[BITS], const bf_t witness[BITS], const bf_t delta, bf_t out_state[3][BITS],bf_t mul_in0[BITS],bf_t mul_in1[BITS]){
    memcpy(out_state[0],state[1],sizeof(out_state[0]));
    bf_t in0[BITS],in1[BITS];

    for(int i=0;i<BITS;i++){
        in0[i] = bf_add(state[0][i],bf_public_128(id,consts[0]>>i & 1,delta));
        in1[i] = bf_add(state[1][i],bf_public_128(id,consts[1]>>i & 1,delta));
    }
    
    memcpy(mul_in0,in0,sizeof(in0));
    memcpy(mul_in1,in1,sizeof(in1));
    bf_t out[BITS],tmp[BITS];
    memcpy(out,witness,sizeof(out));
    LinP_eval_bitlevel_128(out,P,tmp);
    for(int i=0;i<BITS;i++){
        out_state[1][i] = bf_add(tmp[i],state[2][i]);
        out_state[2][i] = state[0][i];
    }    
}




void inverse_128(const sf_t state[3],const sf_t consts,const sf_t P[BITS],sf_t out_state[3], sf_t *witness, sf_t *mul_in){
	sf_t in = sf_add(state[1],consts);
    sf_t out = sf_inv(in);
    
    if(witness != NULL)
        *witness = out;
    
    if(mul_in != NULL)
        *mul_in = in;

    out_state[0] = state[0];
    out_state[1] = LinP_eval_128(out, P);
    out_state[2] = state[2];
}



void inverse_bitlevel_128(int id,bf_t state[3][BITS],const sf_t consts,const sf_t P[BITS], const bf_t witness[BITS], const bf_t delta, bf_t out_state[3][BITS],bf_t mul_in[BITS]){
    bf_t in[BITS];

    for(int i=0;i<BITS;i++)
        in[i] = bf_add(state[1][i],bf_public_128(id,consts>>i & 1,delta));

    memcpy(mul_in,in,sizeof(in));
    
    bf_t out[BITS],tmp[BITS];
    memcpy(out,witness,sizeof(out));
    LinP_eval_bitlevel_128(out,P,tmp);
    for(int i=0;i<BITS;i++){
        out_state[0][i] = state[0][i];
        out_state[1][i] = tmp[i];
        out_state[2][i] = state[2][i];
    }
}


void mds_128(sandwich_128_param_t* para,const sf_t state[3], sf_t out_state[3]){
    for(int i=0;i<3;i++){
        out_state[i] = 0;
        for(int j=0;j<3;j++){
            out_state[i] = sf_add(out_state[i],
                sf_mul(matMDS[i][j] , state[j])
            );
        }
    }
}

void mds_bitlevel_128(sandwich_128_param_t* para,bf_t state[3][BITS],bf_t out_state[3][BITS]){
    for(int i=0;i<3;i++){
        //out_state[i] = bf128_zero();
        memset(&out_state[i],0,sizeof(out_state[i]));
        for(int j=0;j<3;j++){
            if(matMDS[i][j]==1){
                for(int k=0;k<BITS;k++){
                    out_state[i][k] = bf_add(out_state[i][k],state[j][k]);
                }
            }else{// matMDS[i][j]==2

                for(int k=1;k<BITS;k++)
                    out_state[i][k] = bf_add(out_state[i][k],state[j][k-1]);
                for(int k=0;k<20;k++) if(sf_modulus & (1ULL<<k))
                    out_state[i][k] = bf_add(out_state[i][k],state[j][BITS-1]);

            }
        }
    }
}


void sandwich_128(sandwich_128_param_t* para,sf_t k0, sf_t k1, sf_t out[2],sf_t witness[9],sf_t mul_inputs[13]){

    //sf_t *witness = (sf_t *)w;
    int w_count = 0;
    witness[w_count++] = k0;
    witness[w_count++] = k1;

    int m_count = 0;

    sf_t top_state[4][3];
    sf_t bot_state[4][3];
    sf_t inv_state[3];
    sf_t final_state[3];

    top_state[0][0] = k0;
    top_state[0][1] = sf_add(k0,k1);
    top_state[0][2] = k1;
    for(int i=0;i<topRound;i++){
        feistel_128(top_state[i],top_consts[i],top_P[i],top_state[i+1],&witness[w_count++],&mul_inputs[m_count]);
        m_count+=2;
    }
    
    inverse_128(top_state[topRound],inv_const,inv_P,inv_state,&witness[w_count++],&mul_inputs[m_count++]);
    mds_128(para,inv_state,bot_state[0]);

    for(int i=0;i<bottomRound;i++){
        feistel_128(bot_state[i],bot_consts[i],bot_P[i],bot_state[i+1],&witness[w_count++],&mul_inputs[m_count]);
        m_count+=2;
    }
    
    mds_128(para,bot_state[bottomRound],final_state);

    out[0] = sf_add(final_state[0],k0);
    out[1] = sf_add(final_state[2],k1);
}



void sandwich_bitlevel_128(sandwich_128_param_t* para,int id, const bf_t out[2][BITS],const bf_t witness[9][BITS],const bf_t delta,
                        bf_t mul_inputs[13][BITS],bf_t newk[2][BITS]){

    bf_t k0[BITS],k1[BITS];

    int w_count = 0;
    memcpy(k0,witness[w_count++],sizeof(k0));
    memcpy(k1,witness[w_count++],sizeof(k1));

    int m_count = 0;

    
   
    bf_t top_state[4][3][BITS];
    bf_t bot_state[4][3][BITS];
    bf_t inv_state[3][BITS];
    bf_t final_state[3][BITS];

    for(int i=0;i<BITS;i++){
        top_state[0][0][i] = k0[i];
        top_state[0][1][i] = bf_add(k0[i],k1[i]);
        top_state[0][2][i] = k1[i];
    }

    for(int i=0;i<topRound;i++){
        feistel_bitlevel_128(id,top_state[i],top_consts[i],top_P[i],witness[w_count++],delta,top_state[i+1],&mul_inputs[m_count][0],&mul_inputs[m_count+1][0]);
        m_count+=2;
    }

    inverse_bitlevel_128(id,top_state[topRound],inv_const,inv_P,witness[w_count++],delta,inv_state,&mul_inputs[m_count++][0]);



    mds_bitlevel_128(para,inv_state,bot_state[0]);

    for(int i=0;i<bottomRound;i++){
        feistel_bitlevel_128(id,bot_state[i],bot_consts[i],bot_P[i],witness[w_count++],delta,bot_state[i+1],&mul_inputs[m_count][0],&mul_inputs[m_count+1][0]);
        m_count+=2;
    }    


    mds_bitlevel_128(para,bot_state[bottomRound],final_state);

    //out[0] = bf64_add(final_state[0],k0);
    //out[1] = bf64_add(final_state[2],k1);
    for(int i=0;i<BITS;i++){
        newk[0][i] = bf_add(final_state[0][i],out[0][i]);
        newk[1][i] = bf_add(final_state[2][i],out[1][i]);
    }

}


void bf64_to_bf128_bitlevel(sandwich_128_param_t* para,const bf_t in[BITS], bf_t out[BITS*2]){
    bf_t a = bf_one();

    for(int j=0;j<BITS*2;j++)
        out[j]=bf_zero();
    
    bf_t mask[2];
    mask[0]=bf128_zero();
    mask[1]=bf128_all_one();

    for(int i=0;i<BITS;i++){
        a = para->power_of_alpha[i];
        for(int j=0;j<BITS*2;j++){
            bf64_t bit = bf_get_bit(a,j);
            bf_t m = mask[bit];
            //if(bit)
            out[j]=bf_add(out[j],bf128_and(m,in[i]));
        }
        //a=bf_mul(a,alpha);
    }
}

bf_t bf128_convert_combine(sandwich_128_param_t* para,const bf_t in[BITS]){
    bf_t out[BITS*2];
    bf64_to_bf128_bitlevel(para,in,out);
    bf_t two = bf128_from_bf64(2);
    bf_t a = bf_one();
    bf_t res = bf_zero();

    for(int i=BITS*2-1;i>=0;i--){
        a = para->power_of_two[i];
        //res = res<<1 
        const uint64_t mask = bf128_bit_to_uint64_mask(res, 127);
        res                 = bf128_shift_left_1(res);
        res.values[0] ^= (mask & bf128_modulus);

        res = bf_add(res,out[i]);
    }
    return res;
}


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
