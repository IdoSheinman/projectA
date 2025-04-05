#ifndef LDPC_SYSTEM_H
#define LDPC_SYSTEM_H

#include <vector>
#include <utility> // For std::pair
#include <string>  // For error messages (optional)

// --- Define Code Parameters ---
// These constants define the expected dimensions for matrices
// passed to the functions declared below.
// You might make these parameters to the functions instead if you
// need more flexibility, but for this example, we keep them fixed.
const int LDPC_K = 8;  // Message length
const int LDPC_M = 4;  // Number of parity checks
const int LDPC_N = 12; // Codeword length (N = K + M assumed for systematic)

// --- Type Aliases ---
using Vec = std::vector<int>;
using Matrix = std::vector<Vec>;

// --- Function Declarations ---

/**
 * @brief Encodes a binary message using a specified LDPC Generator matrix (G).
 *        (Adapts the flow of Algorithm 5, simplified for standard LDPC).
 * @param message The binary message vector of length K.
 * @param G The KxN Generator matrix (expected systematic form [I_k | P^T] for easy message recovery later).
 * @return The encoded codeword vector of length N.
 * @throws std::invalid_argument if message size != K or G dimensions are incorrect.
 */
Vec algorithm5_encode_simple(const Vec& message, const Matrix& G);

/**
 * @brief Decodes a received (potentially erroneous) codeword using a specified
 *        LDPC Parity Check matrix (H) and the simple bit-flipping algorithm.
 *        Attempts to recover the original message if decoding is successful.
 *        (Adapts the flow of Algorithm 6, simplified for standard LDPC).
 * @param received_codeword The received binary vector of length N.
 * @param H The MxN Parity Check matrix (must be consistent with G used for encoding).
 * @param max_iterations The maximum number of iterations for the bit-flipping decoder.
 * @return A pair:
 *         - first: The decoded *original message* vector of length K (empty if failed).
 *         - second: A boolean flag indicating success (true) or failure (false).
 * @throws std::invalid_argument if received_codeword size != N or H dimensions are incorrect.
 */
std::pair<Vec, bool> algorithm6_decode_simple(const Vec& received_codeword, const Matrix& H, int max_iterations);


// --- Optional: Declare Specific Matrices (Less Flexible Approach) ---
/*
   // If you prefer to have THE specific H and G matrices globally available,
   // declare them as 'extern const' here and DEFINE them exactly ONCE
   // in the corresponding .cpp file. This makes the header less reusable
   // for different codes.

   extern const Matrix H_example; // Must be defined in ldpc_system.cpp
   extern const Matrix G_example; // Must be defined in ldpc_system.cpp
*/

// --- Matrix and Vector Definitions (Larger Example) ---

// H (Parity-Check Matrix, MxN). Constructed as H = [P | I_m]
const Matrix H_example = {
   // P matrix part (4x8)                       | I_m part (4x4)
   {1, 1, 0, 1, 0, 0, 0, 1,  /* | */            1, 0, 0, 0},
   {0, 1, 1, 0, 1, 0, 1, 0,  /* | */            0, 1, 0, 0},
   {1, 0, 0, 1, 1, 1, 0, 0,  /* | */            0, 0, 1, 0},
   {0, 0, 1, 0, 0, 1, 1, 1,  /* | */            0, 0, 0, 1}
};

// G (Generator Matrix, KxN). Constructed as G = [I_k | P^T]
const Matrix G_example = {
   // I_k part (8x8)                               | P^T part (8x4)
   {1, 0, 0, 0, 0, 0, 0, 0,  /* | */             1, 0, 1, 0}, // Row 0
   {0, 1, 0, 0, 0, 0, 0, 0,  /* | */             1, 1, 0, 0}, // Row 1
   {0, 0, 1, 0, 0, 0, 0, 0,  /* | */             0, 1, 0, 1}, // Row 2
   {0, 0, 0, 1, 0, 0, 0, 0,  /* | */             1, 0, 1, 0}, // Row 3
   {0, 0, 0, 0, 1, 0, 0, 0,  /* | */             0, 1, 1, 0}, // Row 4
   {0, 0, 0, 0, 0, 1, 0, 0,  /* | */             0, 0, 1, 1}, // Row 5
   {0, 0, 0, 0, 0, 0, 1, 0,  /* | */             0, 1, 0, 1}, // Row 6
   {0, 0, 0, 0, 0, 0, 0, 1,  /* | */             1, 0, 0, 1}  // Row 7
};


#endif // LDPC_SYSTEM_H
