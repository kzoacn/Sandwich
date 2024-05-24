/*
 *  SPDX-License-Identifier: MIT
 */

#ifndef OWF_H
#define OWF_H

#include "macros.h"

#include <stdbool.h>
#include <stdint.h>

FAEST_BEGIN_C_DECL

bool owf(const uint8_t* key, const uint8_t* input, uint8_t* output, int lambda); 

bool faest_128s_owf(const uint8_t* key, const uint8_t* input, uint8_t* output); 
bool faest_128f_owf(const uint8_t* key, const uint8_t* input, uint8_t* output);

bool faest_192s_owf(const uint8_t* key, const uint8_t* input, uint8_t* output);
bool faest_192f_owf(const uint8_t* key, const uint8_t* input, uint8_t* output);

bool faest_256s_owf(const uint8_t* key, const uint8_t* input, uint8_t* output);
bool faest_256f_owf(const uint8_t* key, const uint8_t* input, uint8_t* output);

FAEST_END_C_DECL

#endif
