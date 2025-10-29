/// Note this code was created by Gemini.
/// It should handle ldpc encoding / decoding
#include "ldpc.h"

#include <iostream>
#include <vector>
#include <numeric>
#include <algorithm>
#include <iterator>
#include <stdexcept>
#include <utility> // For std::pair
#include <string> // For std::string
#include <random> // For potential random message generation

// --- Define Code Parameters (Larger Example) ---
const int K = 8;
const int M = 4;
const int N = 12; // N = K + M

// --- Helper Functions (Printing, Math - vec_mat_mult_mod2 corrected) ---
void print_vector(const std::string &name, const Vec &vec) {
    std::cout << name << ": [";
    for (size_t i = 0; i < vec.size(); ++i) {
        std::cout << vec[i] << (i == vec.size() - 1 ? "" : ", ");
    }
    std::cout << "]" << std::endl;
}

void print_matrix(const std::string &name, const Matrix &matrix) {
    if (matrix.empty()) {
        std::cout << name << " (empty)" << std::endl;
        return;
    }
    size_t rows = matrix.size();
    size_t cols = (rows > 0 && !matrix[0].empty()) ? matrix[0].size() : 0;
    std::cout << name << " (" << rows << "x" << cols << "):" << std::endl;
    for (size_t i = 0; i < rows; ++i) {
        std::cout << "  [";
        size_t current_cols = matrix[i].size();
        for (size_t j = 0; j < current_cols; ++j) {
            std::cout << matrix[i][j] << (j == current_cols - 1 ? "" : ", ");
        }
        std::cout << "]" << std::endl;
    }
}

// Multiply vector (row) by matrix: result = vec * matrix (mod 2)
// vec size = matrix rows, result size = matrix cols
Vec vec_mat_mult_mod2(const Vec &vec, const Matrix &matrix) {
    if (matrix.empty()) return {};
    size_t matrix_rows = matrix.size();
    if (matrix_rows == 0) return {};
    size_t matrix_cols = matrix[0].size();
    if (vec.size() != matrix_rows) {
        throw std::invalid_argument("vec_mat_mult_mod2: Vector size (" + std::to_string(vec.size()) +
                                    ") must match matrix rows (" + std::to_string(matrix_rows) + ").");
    }

    Vec result(matrix_cols, 0);
    for (size_t j = 0; j < matrix_cols; ++j) {
        // Iterate over columns of matrix (output vector elements)
        int sum = 0;
        for (size_t i = 0; i < matrix_rows; ++i) {
            // Iterate over rows of matrix (vector elements)
            // ** CORRECTED INDEXING **
            sum += vec[i] * matrix[i][j];
        }
        result[j] = sum % 2;
    }
    return result;
}


// Multiply matrix by vector (column): result = matrix * vec (mod 2)
// matrix cols = vec size, result size = matrix rows
Vec mat_vec_mult_mod2(const Matrix &matrix, const Vec &vec) {
    if (matrix.empty()) return {};
    size_t matrix_rows = matrix.size();
    if (matrix_rows == 0) return {};
    size_t matrix_cols = matrix[0].size();
    if (vec.size() != matrix_cols) {
        throw std::invalid_argument("mat_vec_mult_mod2: Vector size (" + std::to_string(vec.size()) +
                                    ") must match matrix columns (" + std::to_string(matrix_cols) + ").");
    }

    Vec result(matrix_rows, 0);
    for (size_t i = 0; i < matrix_rows; ++i) {
        // Iterate over rows of matrix (output vector elements)
        int sum = 0;
        for (size_t j = 0; j < matrix_cols; ++j) {
            // Iterate over columns of matrix (vector elements)
            sum += matrix[i][j] * vec[j];
        }
        result[i] = sum % 2;
    }
    return result;
}


// --- Bit-Flipping Decoder (Unchanged) ---
std::pair<Vec, bool> decode_bit_flipping(const Matrix &H, const Vec &received_vector, int max_iterations) {
    if (H.empty() || received_vector.size() != H[0].size()) {
        throw std::invalid_argument("Invalid H matrix or received vector size for decoder.");
    }
    size_t n_bits = received_vector.size();
    size_t m_checks = H.size();
    Vec current_vector = received_vector;

    for (int iter = 0; iter < max_iterations; ++iter) {
        Vec syndrome = mat_vec_mult_mod2(H, current_vector);
        if (std::all_of(syndrome.begin(), syndrome.end(), [](int s) { return s == 0; })) {
            return {current_vector, true};
        }
        Vec flip_counts(n_bits, 0);
        int total_unsatisfied = 0;
        for (size_t i = 0; i < m_checks; ++i) {
            if (syndrome[i] == 1) {
                total_unsatisfied++;
                for (size_t j = 0; j < n_bits; ++j) {
                    if (H[i][j] == 1) flip_counts[j]++;
                }
            }
        }
        if (total_unsatisfied == 0) {
            // Should not happen if syndrome not zero, but safety check
            return {current_vector, false}; // Stuck
        }
        int max_count = 0;
        if (!flip_counts.empty()) max_count = *std::max_element(flip_counts.begin(), flip_counts.end());
        if (max_count == 0) {
            return {current_vector, false}; // Failed - no bit helps
        }
        int bit_to_flip = std::distance(flip_counts.begin(),
                                        std::find(flip_counts.begin(), flip_counts.end(), max_count));
        current_vector[bit_to_flip] = 1 - current_vector[bit_to_flip];
    }
    Vec final_syndrome = mat_vec_mult_mod2(H, current_vector);
    bool success = std::all_of(final_syndrome.begin(), final_syndrome.end(), [](int s) { return s == 0; });
    return {current_vector, success};
}


// --- Algorithm 5 Adaptation (Encoding - Unchanged Flow) ---
Vec algorithm5_encode_simple(const Vec &message, const Matrix &G) {
    // std::cout << "\n--- Algorithm 5 (Encoding Adaptation) ---" << std::endl;
    // print_vector("Input message (x, K=" + std::to_string(message.size()) + ")", message);
    if (message.size() != K) {
        throw std::invalid_argument("Algorithm 5: Input message length incorrect.");
    }
    // std::cout << "Step 1: Obtain message oligos (Skipped B-BIC)" << std::endl;
    // std::cout << "Step 2: LDPC Encode (Simplified: c = message * G)" << std::endl;
    Vec codeword = vec_mat_mult_mod2(message, G);
    // std::cout << "Step 3-12: Skipped steps related to B-BIC/RSB/GC" << std::endl;
    // std::cout << "Step 13: Obtain final codeword" << std::endl;
    //    print_vector("Encoded Codeword (N=" + std::to_string(codeword.size()) + ")", codeword);
    return codeword;
}

// --- Algorithm 6 Adaptation (Decoding - Includes message recovery) ---
std::pair<Vec, bool> algorithm6_decode_simple(const Vec &received_codeword, const Matrix &H, int max_iterations) {
    // std::cout << "\n--- Algorithm 6 (Decoding Adaptation) ---" << std::endl;
    //    print_vector("Input received codeword (N=" + std::to_string(received_codeword.size()) + ")", received_codeword);
    if (received_codeword.size() != N) {
        throw std::invalid_argument("Algorithm 6: Received codeword length incorrect.");
    }

    // std::cout << "Initialization (Skipped d''_j)" << std::endl;
    // std::cout << "Steps 1-8: Perform LDPC decoding (Simplified)" << std::endl;
    std::pair<Vec, bool> decoding_result = decode_bit_flipping(H, received_codeword, max_iterations);
    Vec &estimated_codeword = decoding_result.first;
    bool ldpc_success = decoding_result.second;

    Vec final_syndrome = mat_vec_mult_mod2(H, estimated_codeword);
    bool is_valid_codeword = std::all_of(final_syndrome.begin(), final_syndrome.end(), [](int s) { return s == 0; });

    if (!ldpc_success || !is_valid_codeword) {
        // std::cout << "Step 9/10: LDPC Decoding FAILED or result is not a valid codeword. Cannot recover message." << std::endl;
        // print_vector("Final codeword state", estimated_codeword);
        // print_vector("Resulting Syndrome", final_syndrome);
        return {{}, false};
    }

    // std::cout << "Step 9: Obtain estimated message oligos (from successfully decoded codeword)" << std::endl;
    // print_vector("Estimated Corrected Codeword (ĉ_j)", estimated_codeword);

    // --- FINAL STEP: DECODE ORIGINAL MESSAGE ---
    // std::cout << "Step 10: Decode Original Message (Extract first K=" << K << " bits due to systematic G)" << std::endl;
    if (estimated_codeword.size() < K) {
        std::cerr << "Error: Estimated codeword too short (" << estimated_codeword.size() << ") to extract message!" <<
                std::endl;
        return {{}, false};
    }
    // Because G = [I_k | P^T], the first K bits of the codeword are the message bits
    Vec decoded_message(estimated_codeword.begin(), estimated_codeword.begin() + K);

    // print_vector(">>> Decoded Original Message (x) <<<", decoded_message);
    return {decoded_message, true};
}
