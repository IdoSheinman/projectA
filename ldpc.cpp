/*
 * This code is based on chatGPT. But it implements a completely standard LDPC
 * code and there is no novel or original ideas in it that may be traced to a specific
 * source chatGPT was trained on probably.
 */
#include <iostream>
#include <vector>
#include <numeric>   // For std::accumulate (optional, can use loops)
#include <algorithm> // For std::all_of, std::max_element, std::find
#include <iterator>  // For std::distance
#include <stdexcept> // For std::invalid_argument (optional error check)
#include <utility>   // For std::pair

// --- Define Code Parameters ---
const int M = 3; // Number of parity checks (rows in H)
const int N = 6; // Codeword length (columns in H, columns in G)
const int K = 3; // Message length (rows in G) (N - M for this systematic code)

// --- Matrix and Vector Definitions ---

// H (Parity-Check Matrix, MxN)
const Matrix H_example = {
    {1, 1, 0, 1, 0, 0},
    {0, 1, 1, 0, 1, 0},
    {1, 0, 1, 0, 0, 1}
};

// G (Generator Matrix, KxN)
const Matrix G_example = {
    {1, 0, 0, 1, 0, 1},
    {0, 1, 0, 1, 1, 0},
    {0, 0, 1, 0, 1, 1}
};

// --- Helper Functions ---

// Function to print a vector
void print_vector(const std::string& name, const Vec& vec) {
    std::cout << name << ": [";
    for (size_t i = 0; i < vec.size(); ++i) {
        std::cout << vec[i] << (i == vec.size() - 1 ? "" : ", ");
    }
    std::cout << "]" << std::endl;
}

// Function to print a matrix
void print_matrix(const std::string& name, const Matrix& matrix) {
    if (matrix.empty()) {
        std::cout << name << " (empty)" << std::endl;
        return;
    }
    size_t rows = matrix.size();
    size_t cols = matrix[0].size();
    std::cout << name << " (" << rows << "x" << cols << "):" << std::endl;
    for (size_t i = 0; i < rows; ++i) {
        std::cout << "  [";
        for (size_t j = 0; j < cols; ++j) {
            std::cout << matrix[i][j] << (j == cols - 1 ? "" : ", ");
        }
        std::cout << "]" << std::endl;
    }
}

// Multiply vector (row) by matrix: result = vec * matrix (mod 2)
// vec size = matrix rows, result size = matrix cols
Vec vec_mat_mult_mod2(const Vec& vec, const Matrix& matrix) {
    if (matrix.empty()) return {};
    size_t matrix_rows = matrix.size();
    size_t matrix_cols = matrix[0].size();
    if (vec.size() != matrix_rows) {
        throw std::invalid_argument("Vector size must match matrix rows for vec*mat multiplication.");
    }

    Vec result(matrix_cols, 0);
    for (size_t j = 0; j < matrix_cols; ++j) {
        int sum = 0;
        for (size_t i = 0; i < matrix_rows; ++i) {
            sum += vec[i] * matrix[i][j];
        }
        result[j] = sum % 2;
    }
    return result;
}

// Multiply matrix by vector (column): result = matrix * vec (mod 2)
// matrix cols = vec size, result size = matrix rows
Vec mat_vec_mult_mod2(const Matrix& matrix, const Vec& vec) {
    if (matrix.empty()) return {};
    size_t matrix_rows = matrix.size();
    size_t matrix_cols = matrix[0].size();
     if (vec.size() != matrix_cols) {
        throw std::invalid_argument("Vector size must match matrix columns for mat*vec multiplication.");
    }

    Vec result(matrix_rows, 0);
    for (size_t i = 0; i < matrix_rows; ++i) {
        int sum = 0;
        for (size_t j = 0; j < matrix_cols; ++j) {
            sum += matrix[i][j] * vec[j];
        }
        result[i] = sum % 2;
    }
    return result;
}

// --- LDPC Encoding ---
Vec encode(const Vec& message, const Matrix& G) {
    if (message.size() != G.size()) {
         throw std::invalid_argument("Message length must match Generator matrix rows.");
    }
    std::cout << "\n--- Encoding ---" << std::endl;
    print_vector("Input Message (k)", message);
    Vec codeword = vec_mat_mult_mod2(message, G);
    print_vector("Encoded Codeword (n)", codeword);
    return codeword;
}

// --- LDPC Verification (Checks if H*c^T == 0) ---
// Returns pair: {is_valid, syndrome_vector}
std::pair<bool, Vec> verify(const Matrix& H, const Vec& codeword) {
     if (codeword.size() != H[0].size()) {
        throw std::invalid_argument("Codeword length must match Parity Check matrix columns.");
     }
     Vec syndrome = mat_vec_mult_mod2(H, codeword);
     bool is_valid = std::all_of(syndrome.begin(), syndrome.end(), [](int s){ return s == 0; });
     return {is_valid, syndrome};
}


// --- LDPC Bit-Flipping Decoder ---
// Returns pair: {decoded_vector_attempt, success_flag}
std::pair<Vec, bool> decode_bit_flipping(const Matrix& H, const Vec& received_vector, int max_iterations) {
    if (H.empty() || received_vector.size() != H[0].size()) {
         throw std::invalid_argument("Invalid H matrix or received vector size.");
    }
    size_t n_bits = received_vector.size();
    size_t m_checks = H.size();

    Vec current_vector = received_vector; // Work on a copy

    std::cout << "\n--- Starting Bit-Flipping Decoder ---" << std::endl;
    print_vector("Initial Received Vector", current_vector);

    for (int iter = 0; iter < max_iterations; ++iter) {
        std::cout << "\n--- Iteration " << iter + 1 << " ---" << std::endl;

        // 1. Calculate Syndrome
        Vec syndrome = mat_vec_mult_mod2(H, current_vector);
        print_vector("Syndrome (s = H*y^T)", syndrome);

        // 2. Check for Success (Syndrome is all zero?)
        if (std::all_of(syndrome.begin(), syndrome.end(), [](int s){ return s == 0; })) {
            std::cout << "Syndrome is zero. Decoding successful!" << std::endl;
            return {current_vector, true}; // Success
        }

        // 3 & 4. Count Votes for Unsatisfied Checks
        Vec flip_counts(n_bits, 0); // Counts how many unsatisfied checks each bit participates in
        std::cout << "Unsatisfied checks (indices): [";
        bool first_unsat = true;
        for (size_t i = 0; i < m_checks; ++i) { // Iterate through checks (rows of H)
            if (syndrome[i] == 1) { // If check 'i' is unsatisfied
                 if (!first_unsat) std::cout << ", ";
                 std::cout << i;
                 first_unsat = false;
                // Add votes to bits involved in this check
                for (size_t j = 0; j < n_bits; ++j) { // Iterate through bits (columns of H)
                    if (H[i][j] == 1) {
                        flip_counts[j]++;
                    }
                }
            }
        }
        std::cout << "]" << std::endl;
        print_vector("Flip Counts (votes per bit)", flip_counts);

        // 5. Find Bit(s) with Maximum Votes
        int max_count = 0;
        if (!flip_counts.empty()) {
             max_count = *std::max_element(flip_counts.begin(), flip_counts.end());
        }


        if (max_count == 0) {
            // Should not happen if syndrome is non-zero and H is valid,
            // but good to check. Indicates algorithm is stuck.
            std::cout << "No bit participates significantly in unsatisfied checks (max_count=0). Cannot proceed." << std::endl;
            return {current_vector, false}; // Failed
        }

        // Find the *first* bit index with the maximum count
        int bit_to_flip = -1;
        auto it = std::find(flip_counts.begin(), flip_counts.end(), max_count);
        if (it != flip_counts.end()) {
            bit_to_flip = std::distance(flip_counts.begin(), it);
        }

        // Print all candidates (optional but informative)
        std::cout << "Highest count is " << max_count << " for bit(s): [";
        bool first_max = true;
        for(size_t j = 0; j < n_bits; ++j) {
            if (flip_counts[j] == max_count) {
                if (!first_max) std::cout << ", ";
                std::cout << j;
                first_max = false;
            }
        }
        std::cout << "]" << std::endl;

        if (bit_to_flip == -1) { // Should technically not happen if max_count > 0
             std::cerr << "Error: Could not find index for max_count." << std::endl;
             return {current_vector, false};
        }


        // 6. Flip the Selected Bit
        std::cout << "Flipping bit at index: " << bit_to_flip << std::endl;
        current_vector[bit_to_flip] = 1 - current_vector[bit_to_flip]; // Flip 0 to 1 or 1 to 0
        print_vector("Current vector after flip", current_vector);

    } // End of iteration loop

    // 7. Termination (Max iterations reached without success)
    std::cout << "\nMaximum iterations (" << max_iterations << ") reached. Decoding failed." << std::endl;
    return {current_vector, false}; // Failed
}


// --- Main Program ---
int main() {
    Vec message = {1, 0, 1};
    Vec codeword;
    Vec codeword_with_error;
    Vec decoded_codeword;
    std::pair<bool, Vec> verification_result;

    // --- Demonstrate Encoding ---
    print_matrix("H Matrix", H_example);
    print_matrix("G Matrix", G_example);
    codeword = encode(message, G_example);

    // --- Verify the Original Codeword ---
    std::cout << "\n--- Verifying Original Codeword ---" << std::endl;
    print_vector("Original Codeword", codeword);
    verification_result = verify(H_example, codeword);
    print_vector("Syndrome (H*c^T)", verification_result.second);
    std::cout << "Is original codeword valid? " << (verification_result.first ? "Yes" : "No") << std::endl;


    // --- Introduce an Error ---
    std::cout << "\n--- Introducing Error ---" << std::endl;
    codeword_with_error = codeword; // Copy
    int error_pos = 2;
    if (error_pos >= 0 && error_pos < N) {
        codeword_with_error[error_pos] = 1 - codeword_with_error[error_pos]; // Flip bit
    }
    print_vector("Codeword with error", codeword_with_error);

    // --- Verify the Erroneous Codeword ---
     std::cout << "\n--- Verifying Erroneous Codeword ---" << std::endl;
    verification_result = verify(H_example, codeword_with_error);
    print_vector("Syndrome (H*c^T)", verification_result.second);
    std::cout << "Is erroneous codeword valid? " << (verification_result.first ? "Yes" : "No") << std::endl;


    // --- Attempt to Decode/Correct the Error ---
    int max_iter = 5;
    std::pair<Vec, bool> decoding_result = decode_bit_flipping(H_example, codeword_with_error, max_iter);
    decoded_codeword = decoding_result.first;
    bool success = decoding_result.second;


    // --- Check Decoding Result ---
    std::cout << "\n--- Decoding Result ---" << std::endl;
    print_vector("Final Decoded Vector", decoded_codeword);
    if (success) {
        std::cout << "Decoding reported SUCCESS." << std::endl;
        // Optional: Compare with original codeword
        if (decoded_codeword == codeword) {
            std::cout << "Corrected vector matches the original codeword." << std::endl;
        } else {
             // This could happen if BF converges to a different valid codeword
             std::cout << "Warning: Decoded vector is valid but differs from original." << std::endl;
        }
    } else {
        std::cout << "Decoding reported FAILURE." << std::endl;
         // Verify the final state again
        verification_result = verify(H_example, decoded_codeword);
        print_vector("Syndrome of final vector", verification_result.second);
    }

    return 0; // Indicate successful execution
}
