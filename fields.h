/*
 *  SPDX-License-Identifier: MIT
 */

#ifndef FAEST_FIELDS_H
#define FAEST_FIELDS_H

#include "macros.h"
#include "endian_compat.h"

#include <stdint.h>
#include <string.h>

FAEST_BEGIN_C_DECL



// GF(2^8) with X^8 + X^4 + X^3 + X^1 + 1
#define bf8_modulus (UINT8_C((1 << 4) | (1 << 3) | (1 << 1) | 1))
// GF(2^64) with X^64 + X^4 + X^3 + X^1 + 1
#define bf64_modulus (UINT64_C((1 << 4) | (1 << 3) | (1 << 1) | 1))
// GF(2^128) with X^128 + X^7 + X^2 + X^1 + 1
#define bf128_modulus (UINT64_C((1 << 7) | (1 << 2) | (1 << 1) | 1))
// GF(2^192) with X^192 + X^7 + X^2 + X^1 + 1
#define bf192_modulus (UINT64_C((1 << 7) | (1 << 2) | (1 << 1) | 1))
// GF(2^256) with X^256 + X^10 + X^5 + X^2 + 1
#define bf256_modulus (UINT64_C((1 << 10) | (1 << 5) | (1 << 2) | 1))

 




typedef uint8_t bf8_t;
typedef uint64_t bf64_t;

typedef struct {
  uint64_t values[2];
} bf128_t;

typedef struct {
  uint64_t values[3];
} bf192_t;

typedef struct {
  uint64_t values[4];
} bf256_t;

// GF(2^8) implementation

ATTR_PURE ATTR_ALWAYS_INLINE static inline bf8_t bf8_load(const uint8_t* src) {
  return *src;
}

ATTR_ALWAYS_INLINE static inline void bf8_store(uint8_t* dst, bf8_t src) {
  *dst = src;
}

ATTR_CONST ATTR_ALWAYS_INLINE static inline bf8_t bf8_zero() {
  return 0;
}

ATTR_CONST ATTR_ALWAYS_INLINE static inline bf8_t bf8_one() {
  return 1;
}

//bf8_t bf8_rand(); 

ATTR_CONST ATTR_ALWAYS_INLINE inline bf8_t bf8_add(bf8_t lhs, bf8_t rhs) {
  return lhs ^ rhs;
}

ATTR_CONST bf8_t bf8_mul(bf8_t lhs, bf8_t rhs);
ATTR_CONST bf8_t bf8_inv(bf8_t lhs);

ATTR_CONST ATTR_ALWAYS_INLINE static inline bf8_t bf8_from_bit(uint8_t bit) {
  return bit & 1;
}

// GF(2^64) implementation

ATTR_PURE ATTR_ALWAYS_INLINE inline bf64_t bf64_load(const uint8_t* src) {
  bf64_t ret;
  memcpy(&ret, src, sizeof(ret));
#if defined(FAEST_IS_BIG_ENDIAN)
  ret = le64toh(ret);
#endif
  return ret;
}

ATTR_ALWAYS_INLINE inline void bf64_store(uint8_t* dst, bf64_t src) {
#if defined(FAEST_IS_BIG_ENDIAN)
  src = htole64(src);
#endif
  memcpy(dst, &src, sizeof(src));
}

ATTR_CONST ATTR_ALWAYS_INLINE static inline bf64_t bf64_zero() {
  return 0;
}

ATTR_CONST ATTR_ALWAYS_INLINE static inline bf64_t bf64_one() {
  return 1;
}

//bf64_t bf64_rand();

ATTR_CONST ATTR_ALWAYS_INLINE static inline bf64_t bf64_add(bf64_t lhs, bf64_t rhs) {
  return lhs ^ rhs;
}



ATTR_CONST ATTR_ALWAYS_INLINE static inline bf64_t bf64_flip(bf64_t x, bf64_t pos) {
  return x ^ ((bf64_t)1 << pos);
}

ATTR_CONST ATTR_ALWAYS_INLINE static inline bf64_t bf64_get_bit(bf64_t x, bf64_t pos) {
  return (x >> pos) & 1;
}
ATTR_CONST ATTR_ALWAYS_INLINE static inline bf64_t bf64_and(bf64_t lhs, bf64_t rhs) {
  return lhs & rhs;
}


ATTR_CONST bf64_t bf64_mul(bf64_t lhs, bf64_t rhs);
ATTR_CONST bf64_t bf64_inv(bf64_t lhs);

ATTR_CONST ATTR_ALWAYS_INLINE static inline bf64_t bf64_from_bf64(bf64_t src) {
  bf64_t ret = src;
  return ret;
}

ATTR_CONST ATTR_ALWAYS_INLINE static inline bf64_t bf64_from_bit(uint8_t bit) {
  return bit & 1;
}

// bf64_count

ATTR_CONST ATTR_ALWAYS_INLINE static inline unsigned int bf64_count(bf64_t x) {
  unsigned int count = 0;
  for (unsigned int i = 0; i < 64; ++i) {
    if(x& (1ULL << i)) 
      count++;
  }
  return count;
} 
// GF(2^128) implementation



ATTR_CONST
static inline bf128_t bf128_and(bf128_t lhs, bf128_t rhs) {
  for (unsigned int i = 0; i != ARRAY_SIZE(lhs.values); ++i) {
    lhs.values[i] &= rhs.values[i];
  }
  return lhs;
}

ATTR_CONST ATTR_ALWAYS_INLINE static inline bf128_t bf128_flip(bf128_t x, bf64_t pos) {
  if (pos < 64) {
    x.values[0] ^= (UINT64_C(1) << pos);
  } else {
    x.values[1] ^= (UINT64_C(1) << (pos - 64));
  }
  return x;
}

ATTR_CONST
static inline bf128_t bf128_shift_left_1(bf128_t value) {
  value.values[1] = (value.values[1] << 1) | (value.values[0] >> 63);
  value.values[0] = value.values[0] << 1;
  return value;
}

ATTR_PURE ATTR_ALWAYS_INLINE static inline bf128_t bf128_load(const uint8_t* src) {
  bf128_t ret; 
  memcpy(&ret, src, sizeof(ret)); 
  return ret;
}

ATTR_ALWAYS_INLINE static inline void bf128_store(uint8_t* dst, bf128_t src) { 
  memcpy(dst, &src, sizeof(src)); 
}

ATTR_CONST ATTR_ALWAYS_INLINE static inline bf128_t bf128_from_bf64(bf64_t src) {
  bf128_t ret;
  ret.values[0] = src;
  ret.values[1] = 0;
  return ret;
}

ATTR_CONST ATTR_ALWAYS_INLINE static inline bf64_t bf128_get_bit(bf128_t x, bf64_t pos) {
  return (x.values[pos / 64] >> (pos % 64)) & 1;
}

ATTR_CONST ATTR_ALWAYS_INLINE static inline bf64_t bf256_get_bit(bf256_t x, bf64_t pos) {
  return (x.values[pos / 64] >> (pos % 64)) & 1;
}

ATTR_CONST ATTR_ALWAYS_INLINE static inline unsigned int bf128_count(bf128_t x) {
  unsigned int count = 0;
  for (unsigned int i = 0; i < 128; ++i) {
    if(bf128_get_bit(x, i)) 
      count++;
  }
  return count;
}
ATTR_CONST ATTR_ALWAYS_INLINE static inline bf128_t bf128_from_bf8(bf8_t src) {
  bf128_t ret;
  ret.values[0] = src;
  ret.values[1] = 0;
  return ret;
}

ATTR_CONST ATTR_ALWAYS_INLINE static inline bf128_t bf128_from_bit(uint8_t bit) {
  return bf128_from_bf8(bit & 1);
}

ATTR_CONST ATTR_ALWAYS_INLINE static inline bf128_t bf128_zero() {
  bf128_t r = {0};
  return r;
}

ATTR_CONST ATTR_ALWAYS_INLINE static inline bf128_t bf128_one() {
  bf128_t r = {{1, 0}};
  return r;
}
 
ATTR_CONST static inline bf128_t bf128_add(bf128_t lhs, bf128_t rhs) {
  for (unsigned int i = 0; i != ARRAY_SIZE(lhs.values); ++i) {
    lhs.values[i] ^= rhs.values[i];
  }
  return lhs;
}

ATTR_CONST bf128_t bf128_mul(bf128_t lhs, bf128_t rhs);
//ATTR_CONST bf128_t bf128_mul_bit(bf128_t lhs, uint8_t rhs);
ATTR_CONST bf128_t bf128_inv(bf128_t lhs);
ATTR_PURE bf128_t bf128_sum_poly(const bf128_t* xs);

// GF(2^192) implemenation

ATTR_PURE ATTR_ALWAYS_INLINE static inline bf192_t bf192_load(const uint8_t* src) {
  bf192_t ret;
#if defined(FAEST_IS_BIG_ENDIAN)
  for (unsigned int i = 0; i != ARRAY_SIZE(ret.values); ++i, src += sizeof(uint64_t)) {
    memcpy(&ret.values[i], src, sizeof(ret.values[i]));
    ret.values[i] = le64toh(ret.values[i]);
  }
#else
  memcpy(&ret, src, sizeof(ret));
#endif
  return ret;
}

ATTR_ALWAYS_INLINE static inline void bf192_store(uint8_t* dst, bf192_t src) {
#if defined(FAEST_IS_BIG_ENDIAN)
  for (unsigned int i = 0; i != ARRAY_SIZE(src.values); ++i, dst += sizeof(uint64_t)) {
    uint64_t tmp = htole64(src.values[i]);
    memcpy(dst, &tmp, sizeof(tmp));
  }
#else
  memcpy(dst, &src, sizeof(src));
#endif
}

ATTR_CONST ATTR_ALWAYS_INLINE static inline bf192_t bf192_from_bf64(bf64_t src) {
  bf192_t ret;
  ret.values[0] = src;
  ret.values[1] = ret.values[2] = 0;
  return ret;
}

ATTR_CONST ATTR_ALWAYS_INLINE static inline bf192_t bf192_from_bf8(bf8_t src) {
  bf192_t ret;
  ret.values[0] = src;
  ret.values[1] = ret.values[2] = 0;
  return ret;
}

ATTR_CONST ATTR_ALWAYS_INLINE static inline bf192_t bf192_from_bit(uint8_t bit) {
  return bf192_from_bf8(bit & 1);
}

ATTR_CONST ATTR_ALWAYS_INLINE static inline bf192_t bf192_zero() {
  bf192_t r = {0};
  return r;
}

ATTR_CONST ATTR_ALWAYS_INLINE static inline bf192_t bf192_one() {
  bf192_t r = {{1, 0, 0}};
  return r;
}


ATTR_CONST static inline bf192_t bf192_add(bf192_t lhs, bf192_t rhs) {
  for (unsigned int i = 0; i != ARRAY_SIZE(lhs.values); ++i) {
    lhs.values[i] ^= rhs.values[i];
  }
  return lhs;
}

ATTR_CONST bf192_t bf192_mul(bf192_t lhs, bf192_t rhs);
//ATTR_CONST bf192_t bf192_mul_bit(bf192_t lhs, uint8_t rhs);
ATTR_CONST bf192_t bf192_inv(bf192_t lhs);
ATTR_PURE bf192_t bf192_sum_poly(const bf192_t* xs);

// GF(2^256) implementation


ATTR_CONST
static inline bf256_t bf256_shift_left_1(bf256_t value) {
  value.values[3] = (value.values[3] << 1) | (value.values[2] >> 63);
  value.values[2] = (value.values[2] << 1) | (value.values[1] >> 63);
  value.values[1] = (value.values[1] << 1) | (value.values[0] >> 63);
  value.values[0] = value.values[0] << 1;
  return value;
}

ATTR_PURE ATTR_ALWAYS_INLINE static inline bf256_t bf256_load(const uint8_t* src) {
  bf256_t ret;
#if defined(FAEST_IS_BIG_ENDIAN)
  for (unsigned int i = 0; i != ARRAY_SIZE(ret.values); ++i, src += sizeof(uint64_t)) {
    memcpy(&ret.values[i], src, sizeof(ret.values[i]));
    ret.values[i] = le64toh(ret.values[i]);
  }
#else
  memcpy(&ret, src, sizeof(ret));
#endif
  return ret;
}

ATTR_ALWAYS_INLINE static inline void bf256_store(uint8_t* dst, bf256_t src) {
#if defined(FAEST_IS_BIG_ENDIAN)
  for (unsigned int i = 0; i != ARRAY_SIZE(src.values); ++i, dst += sizeof(uint64_t)) {
    uint64_t tmp = htole64(src.values[i]);
    memcpy(dst, &tmp, sizeof(tmp));
  }
#else
  memcpy(dst, &src, sizeof(src));
#endif
}

ATTR_CONST ATTR_ALWAYS_INLINE static inline bf256_t bf256_from_bf64(bf64_t src) {
  bf256_t ret;
  ret.values[0] = src;
  ret.values[1] = ret.values[2] = ret.values[3] = 0;
  return ret;
}

ATTR_CONST ATTR_ALWAYS_INLINE static inline bf256_t bf256_from_bf8(bf8_t src) {
  bf256_t ret;
  ret.values[0] = src;
  ret.values[1] = ret.values[2] = ret.values[3] = 0;
  return ret;
}

ATTR_CONST ATTR_ALWAYS_INLINE static inline bf256_t bf256_from_bit(uint8_t bit) {
  return bf256_from_bf8(bit & 1);
}

ATTR_CONST ATTR_ALWAYS_INLINE static inline bf256_t bf256_zero() {
  bf256_t r = {0};
  return r;
}

ATTR_CONST ATTR_ALWAYS_INLINE static inline bf256_t bf256_one() {
  bf256_t r = {{1, 0, 0, 0}};
  return r;
}
 
ATTR_CONST static inline bf256_t bf256_add(bf256_t lhs, bf256_t rhs) {
  for (unsigned int i = 0; i != ARRAY_SIZE(lhs.values); ++i) {
    lhs.values[i] ^= rhs.values[i];
  }
  return lhs;
}

ATTR_CONST bf256_t bf256_mul(bf256_t lhs, bf256_t rhs);
//ATTR_CONST bf256_t bf256_mul_bit(bf256_t lhs, uint8_t rhs);
ATTR_CONST bf256_t bf256_inv(bf256_t lhs);
ATTR_PURE bf256_t bf256_sum_poly(const bf256_t* xs);

FAEST_END_C_DECL

#endif
