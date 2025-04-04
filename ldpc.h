//
// Created by nadav on 4/3/2025.
//
#include <vector>
// Type aliases for clarity
using Vec = std::vector<int>;
using Matrix = std::vector<Vec>;

#ifndef LDPC_H
#define LDPC_H



std::pair<Vec, bool> decode_bit_flipping(const Matrix& H, const Vec& received_vector, int max_iterations);



#endif //LDPC_H
