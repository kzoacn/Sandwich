/*
 *  SPDX-License-Identifier: MIT
 */

#include "faest_aes.h"
#include "instances.hpp"
#include "randomness.h"
#include "tvs_aes.hpp"

#include <boost/test/unit_test.hpp>
#include <boost/test/data/test_case.hpp>
#include <cmath>
#include <array>
#include <vector>

BOOST_AUTO_TEST_SUITE(test_aes_prove)

BOOST_DATA_TEST_CASE(aes_prove_verify, all_parameters, param_id) {
  BOOST_TEST_CONTEXT("Parameter set: " << faest_get_param_name(param_id)) {
    
  }
}

BOOST_AUTO_TEST_SUITE_END()
