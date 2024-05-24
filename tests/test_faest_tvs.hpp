/*
 *  SPDX-License-Identifier: MIT
 */

#ifndef TEST_FAEST_TVS_HPP
#define TEST_FAEST_TVS_HPP

#include <array>
#include <cstdint>
#include <string>

namespace faest_tvs {

  extern const std::string message;
 namespace faest_128s_tvs {
extern const std::array<uint8_t, 32> packed_sk;
extern const std::array<uint8_t, 86> packed_pk;
extern const std::array<uint8_t, 16> randomness;
extern const uint64_t signature_checksum;
  } // namespace faest_128s_tvs

  namespace faest_128f_tvs {
extern const std::array<uint8_t, 32> packed_sk;
extern const std::array<uint8_t, 86> packed_pk;
extern const std::array<uint8_t, 16> randomness;
extern const uint64_t signature_checksum;
  } // namespace faest_128f_tvs

  namespace faest_192s_tvs {
extern const std::array<uint8_t, 48> packed_sk;
extern const std::array<uint8_t, 142> packed_pk;
extern const std::array<uint8_t, 24> randomness;
extern const uint64_t signature_checksum;
  } // namespace faest_192s_tvs

  namespace faest_192f_tvs {
extern const std::array<uint8_t, 48> packed_sk;
extern const std::array<uint8_t, 142> packed_pk;
extern const std::array<uint8_t, 24> randomness;
extern const uint64_t signature_checksum;
  } // namespace faest_192f_tvs

  namespace faest_256s_tvs {
extern const std::array<uint8_t, 64> packed_sk;
extern const std::array<uint8_t, 190> packed_pk;
extern const std::array<uint8_t, 32> randomness;
extern const uint64_t signature_checksum;
  } // namespace faest_256s_tvs

  namespace faest_256f_tvs {
extern const std::array<uint8_t, 64> packed_sk;
extern const std::array<uint8_t, 190> packed_pk;
extern const std::array<uint8_t, 32> randomness;
extern const uint64_t signature_checksum;
  } // namespace faest_256f_tvs


} // namespace faest_tvs

#endif // TEST_FAEST_TVS_HPP
