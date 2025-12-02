//
// Created by idosh on 21/11/2025.
//

#include "parallel_algo.h"
#include <iostream>
#include <errno.h>
#include <signal.h>
#include <sstream>
#include <fstream>
#include <algorithm>
#include <list>
#include <math.h>
#include <utility>
#include <vector>
#include <random>
#include "ldpc.h"
using namespace std;

#define NUC_A 0
#define NUC_T 1
#define NUC_C 2
#define NUC_G 3

// GLOBAL VAR
vector<bool> data_vec_gl;
vector<pair<int, vector<bool> > > bbic_enc_oligo_vec_gl;
vector<bool> restored_vec_gl;
int size_m_gl;
int size_k_gl;
int size_beta_gl;
int parallel_number_gl;
int dummy_bit_counter_gl;
int used_all_data_gl;
int kb_idxs[8];
void knuth_balance(int idx);
int real_kb_idxs(int stored_idx);

void printNucVectorBool(vector<bool> vec) {
    auto it = vec.begin();
    for (; (it != vec.end()) && (it + 1 != vec.end()); it = it + 2) {
        if ((!*it) && (!*(it + 1))) {
            cout << "A" << " ";
        }
        if ((!*it) && (*(it + 1))) {
            cout << "T" << " ";
        }
        if ((*it) && (!*(it + 1))) {
            cout << "C" << " ";
        }
        if (*it && *(it + 1)) {
            cout << "G" << " ";
        }
    }
    cout << endl;
}


bool areEqual(const vector<bool> &vec1, const vector<bool> &vec2, bool verbose = false) {
    if (verbose) {
        cout << "Differences between the vectors:" << endl;
    }
    bool hasDifferences = false;
    for (size_t i = 0; i < vec1.size(); ++i) {
        if (vec1[i] != vec2[i]) {
            if (verbose) {
                cout << "Index " << i << ": vec1 = " << vec1[i] << ", vec2 = " << vec2[i] << endl;
            }
            hasDifferences = true;
        }
    }

    if (verbose && !hasDifferences) {
        cout << "The vectors are identical." << endl;
    }
    return !hasDifferences;
}

bool add_dummy_bits()
/**
 * This gives us dummy bits to read if there is no more data, in a way that doesnt create long running A(/C/T/G) only sequences.
 */
{
    used_all_data_gl = true;
    bool dummy_arr[] = {false, false, false, true, true, false, true, true};
    bool ret = dummy_arr[dummy_bit_counter_gl % 8];
    dummy_bit_counter_gl++;
    return ret;
}

int add_empty_items_bbic(std::vector<bool> &result, int m)
/**
 * Fills result with zeros up to its size.
 */
{
    for (int i = 2 * m + 1; i <= result.size(); i += 2 * m) {
        result.insert(result.begin() + i, false);
    }
    return 0;
}


int restart_param(int m, int k, int p)
/**
 * Initialize the input parameters.
 * :param m: longest allowed RL value.
 * :param k: As appears in the paper.
 */
{
    // BBIC
    size_m_gl = m; //max_sequence_nuc_gl
    size_k_gl = k; //size_oligo_nuc_gl
    size_beta_gl = (size_k_gl + size_m_gl - 1) / size_m_gl - 1;
    parallel_number_gl = p;
    dummy_bit_counter_gl = 0;

    // Knuth
    // We want the knuth balance indexes we write to the oligo
    // to be CG/AT balanced (1 CG and 1 AT)
    kb_idxs[0] = 2; // AC
    kb_idxs[1] = 3; // AG
    kb_idxs[2] = 6; // TC
    kb_idxs[3] = 7; // TG

    kb_idxs[4] = 8; // CA
    kb_idxs[5] = 9; // CT
    kb_idxs[6] = 12; // GA
    kb_idxs[7] = 13; // GT

    return 0;
}

void fill_random_data_vec_gl(int data_size)
/**
 * This function creates random data for us to encode and decode to verify our success.
 */
{
    // Constants for 0 to 1 MB size in bits
    const size_t MAX_BITS = data_size; //128 * 8 * 8; // 1 KB in bits

    // Seed random number generator
    random_device rd; // Random device for seed
    mt19937 gen(rd()); // Mersenne Twister RNG
    uniform_int_distribution<size_t> size_dist(0, MAX_BITS); // Size between 0 and 1MB
    uniform_int_distribution<int> bit_dist(0, 1); // Values 0 or 1

    // Generate a random size
    size_t test_size = size_dist(gen) * 8;

    // Resize the vector
    data_vec_gl.resize(test_size);

    // Fill with random values
    for (size_t i = 0; i < test_size; ++i) {
        data_vec_gl[i] = static_cast<bool>(bit_dist(gen));
    }
}

int copy_data_chunk_without_Q (int ol_i, int idx_s, int idx_e)
/**
 * this writes the raw data into the oligo, while not writing anything
 * for the Q indexes.
 * :param ol_i: the index of the oligo to fill.
 */
{
    auto it = bbic_enc_oligo_vec_gl.begin();
    for (; it != bbic_enc_oligo_vec_gl.end(); it++) {
        if (it->first == ol_i) {
            break;
        }
    }

    if (it == bbic_enc_oligo_vec_gl.end()) {
        return 1;
    }

    for (int l = idx_s; l <= idx_e; l++) {
        if (l < data_vec_gl.size()) {
            it->second[l - idx_s + 4] = data_vec_gl[l]; // +4 for KB idx
        } else {
            it->second[l - idx_s + 4] = add_dummy_bits(); // +4 for KB idx
        }
    }

    return (add_empty_items_bbic((it->second), size_m_gl));

}

bool is_condition_1_satisfied(vector<bool> oli, int position)
/**
 * Check for condition 1 as appears in the paper in page 768.
 * :param oli: The index of the oligo to check.
 * :param position: The index of the position to check (where the Q value is inserted.
 */
{
    int first_upper_pos = position - (2 * size_m_gl + 1);
    int first_lower_pos = position - (2 * size_m_gl);

    bool curr_upper_bit;
    bool curr_lower_bit;
    int current_streak = 1;
    for (int i = 2; i < 2 * (2 * size_m_gl); i += 2) {
        if (first_lower_pos + i >= oli.size()) {
            return false;
        }
        if (i == 2 * size_m_gl) {
            continue;
        }
        curr_upper_bit = oli[first_upper_pos + i];
        curr_lower_bit = oli[first_lower_pos + i];

            if (curr_upper_bit == oli[first_upper_pos + i] &&
                curr_lower_bit == oli[first_lower_pos + i]) {
                    current_streak += 1;
                }
            else {
                    current_streak = 1;
                }

        if (current_streak == size_m_gl) {
            return true;
        }
    }
    return false;
}

bool is_condition_2_satisfied(vector<bool> oli, int position)
/**
 * Check for condition 2 as appears in the paper in page 768.
 * :param oli: The index of the oligo to check.
 * :param position: The index of the position to check (where the Q value is inserted.
 */
{
    return (oli[position - 1] == oli[position - 3]);
}


int add_bbic_bits (int ol_i, int idx_s)
/**
 * This function write into the bbic positions their correct values ( either stops a sequence of repeating ATCG or data )
 * :param ol_i: The  oligo index
 */
{
    int d_i = 0;
    auto it = bbic_enc_oligo_vec_gl.begin();
    for (; it != bbic_enc_oligo_vec_gl.end(); it++) {
        if (it->first == ol_i) {
            break;
        }
    }

    if (it == bbic_enc_oligo_vec_gl.end()) {
        // This shouldnt happen, but for hw purposes not sure return 1 is ideal.
        return -1;
    }

    for (int l = 2 * size_m_gl + 1; l < it->second.size(); l += 2 * size_m_gl) {
        // If need to stop sequence.
        if (is_condition_1_satisfied(it->second, l) && is_condition_2_satisfied(it->second, l)) {
            // TODO: Shouldnt this be l-1
            it->second[l] = !(it->second[l - 2]);
        }
        // Else can fill data
        else {
            if (idx_s + d_i < data_vec_gl.size()) {
                it->second[l] = data_vec_gl[idx_s + d_i];
                if (idx_s + d_i == data_vec_gl.size() - 1) {
                    used_all_data_gl = true;
                }
                d_i++;
            }
            // We want to fill data but there is no more data to fill
            else {
                //cout << "ADD DUMMY BITS" << endl;
                it->second[l] = add_dummy_bits(); //need to see what val at the end
            }
        }
    }
    return d_i;
}


int encode_bbic() {
    int ptr_s = 0;
    int ptr_e = 0;
    int oligo_num = 0;
    used_all_data_gl = false;
    int rest =signed(data_vec_gl.size() - ptr_s);
    if (rest <= 0) {
        used_all_data_gl = true;
    }
    int oligo_data_size = 2 * size_k_gl - size_beta_gl;
    while (!used_all_data_gl) {
        // copy data without bbic bits
        for (int ol_n = oligo_num; ol_n < oligo_num + parallel_number_gl; ol_n++) {
            ptr_e = ptr_s + oligo_data_size - 1;
            vector<bool> oligo(oligo_data_size + 4); // +4 for KB index
            // Add the pair (i, oligo) to result_vec
            bbic_enc_oligo_vec_gl.emplace_back(ol_n, oligo);
            if (copy_data_chunk_without_Q(ol_n, ptr_s, ptr_e))  return 1;
            ptr_s = ptr_e + 1;
        }
        // end that parallel_number_gl oligos fill with data
        // and ptr_s = ptr_e + 1;

        // add bbic bits
        int added = 0;
        for (int ol_n = oligo_num; ol_n < oligo_num + parallel_number_gl; ol_n++) {
            knuth_balance(ol_n);
            added = add_bbic_bits(ol_n, ptr_s); //add q bits

            ptr_s = ptr_s + added;
        }

        rest = signed(data_vec_gl.size() - ptr_s);
        if (rest <= 0) {
            used_all_data_gl = true;
        }
        oligo_num = oligo_num + parallel_number_gl;
    }
    return 0;

}

int decode_bbic()
{
    restored_vec_gl = vector<bool>();
    int oligo_num = 0;
    vector<bool> Q_data = vector<bool>();

    while (oligo_num < bbic_enc_oligo_vec_gl.size()) {
        int first_q = 2 * size_m_gl + 1; // e.g., 7

        // 1. EXTRACT Q DATA
        for (int ol_n = oligo_num; ol_n < oligo_num + parallel_number_gl; ol_n++) {
            if (ol_n >= bbic_enc_oligo_vec_gl.size()) break;

            // Start checking from first_q.
            // Since we pre-allocated, the vector size is correct.
            for (int l = first_q; l < bbic_enc_oligo_vec_gl[ol_n].second.size(); l += 2 * size_m_gl) {

                if (!(is_condition_1_satisfied(bbic_enc_oligo_vec_gl[ol_n].second, l) && is_condition_2_satisfied(
                bbic_enc_oligo_vec_gl[ol_n].second, l))) {
                    Q_data.insert(Q_data.end(), 1, bbic_enc_oligo_vec_gl[ol_n].second[l]);
                }
            }
        }

        // 2. DECODE DATA
        for (int ol_n = oligo_num; ol_n < oligo_num + parallel_number_gl; ol_n++) {
            // Read Header
            int kbt_idx =
                            (bbic_enc_oligo_vec_gl[ol_n].second[0] << 3) | // Index 0 is MSB
                            (bbic_enc_oligo_vec_gl[ol_n].second[1] << 2) |
                            (bbic_enc_oligo_vec_gl[ol_n].second[2] << 1) |
                            (bbic_enc_oligo_vec_gl[ol_n].second[3] << 0); // Index 3 is LSB

            // NOTE: In knuth_balance you wrote:
            // [0] = val & 1
            // [1] = (val & 2) >> 1
            // So [0] is LSB.
            // In your previous decode code you had (vec[0] << 3). This was reversed.
            // I corrected the bitwise logic above to match your knuth_balance write order.

            // Unbalance Payload (skip header 0-3)
            for (int j = 4; j <= real_kb_idxs(kbt_idx) + 4; j += 2) {
                 if (j < bbic_enc_oligo_vec_gl[ol_n].second.size())
                    bbic_enc_oligo_vec_gl[ol_n].second[j] = !bbic_enc_oligo_vec_gl[ol_n].second[j];
            }

            // Extract Data (skip header 0-3)
            for (int l = 4; l < bbic_enc_oligo_vec_gl[ol_n].second.size(); l++) {
                // Skip Q positions
                if (l < first_q || ((l - first_q) % (2 * size_m_gl) != 0)) {
                    restored_vec_gl.insert(restored_vec_gl.end(), 1, bbic_enc_oligo_vec_gl[ol_n].second[l]);
                }
            }
        }

        restored_vec_gl.insert(restored_vec_gl.end(), Q_data.begin(), Q_data.end());
        Q_data = vector<bool>();
        oligo_num = oligo_num + parallel_number_gl;
    }
    return 0;
}

// Knuth Balance
int real_kb_idxs(int stored_idx)
/**
 * After writing CG/AT balanced values, we want to actually use
 * optimal knuth_balance indexes and not use non-optimal indexes just because they are balanced.
 * This converts the optimal for writing to optimal for using indexes.
 */
{
    if (stored_idx == kb_idxs[0]) {
        return 0;
    }
    if (stored_idx == kb_idxs[1]) {
        return 4;
    }
    if (stored_idx == kb_idxs[2]) {
        return 6;
    }
    if (stored_idx == kb_idxs[3]) {
        return 8;
    }
    if (stored_idx == kb_idxs[4]) {
        return 10;
    }
    if (stored_idx == kb_idxs[5]) {
        return 12;
    }
    if (stored_idx == kb_idxs[6]) {
        return 14;
    }
    if (stored_idx == kb_idxs[7]) {
        return 14;
    }
    return -1;
}

int find_best_idx(int idx)
/**
 * This function finds the best index, for a given idx,
 * to knuth balance it out of the possible values.
 * :param idx: The index of the oligo to find the best index for.
 */ {
    vector<bool> olig = bbic_enc_oligo_vec_gl[idx].second;
    int scores[8] = {0};
    for (int i = 0; i < 8; i++) {
        vector<bool> tmp_olig;
        copy(olig.begin(), olig.end(), back_inserter(tmp_olig));
        for (int j = 4; j <= real_kb_idxs(kb_idxs[i]) + 4; j += 2) {
            tmp_olig[j].flip();
        }
        int counter = 0;
        for (int j = 4; j < tmp_olig.size(); j += 2) {
            if (tmp_olig[j] == 0) {
                counter += 1;
            } else {
                counter -= 1;
            }
        }
        scores[i] = abs(counter);
    }
    int index = 0;

    for (int i = 0; i < 8; i++) {
        if (scores[i] < scores[index])
            index = i;
    }
    return kb_idxs[index];
}

void knuth_balance(int idx)
/**
 * Flips bits according to the knuth balance technique and inserts the knuth balance index in the beggining of the
 * given oligo.
 * :param idx: The index of the oligo to balance.
 */ {
    vector<bool> Q_data = vector<bool>();
    vector<bool>& olig = bbic_enc_oligo_vec_gl[idx].second;
    int zero_cnt = 0;
    int one_cnt = 0;
    for (int i = 0; i < olig.size(); i += 2) {
        if (olig[i] == 0) {
            zero_cnt += 1;
        }
        if (olig[i] == 1) {
            one_cnt += 1;
        }
    }
    int idx_to_balance = 0;
    idx_to_balance = find_best_idx(idx);
    for (int i = 4; i <= real_kb_idxs(idx_to_balance) + 4; i += 2) {
        bbic_enc_oligo_vec_gl[idx].second[i] = !bbic_enc_oligo_vec_gl[idx].second[i];
    }
    // 4 index bits, depend on oligo size (this should work with 20)
    olig[0] = ((idx_to_balance & 8) >> 3) & 1; // Bit 3 (MSB)
    olig[1] = ((idx_to_balance & 4) >> 2) & 1; // Bit 2
    olig[2] = ((idx_to_balance & 2) >> 1) & 1; // Bit 1
    olig[3] = (idx_to_balance & 1);            // Bit 0 (LSB)
}

int verify_oligo(vector<bool> oligo)
/**
 * Verifies the given oligo is valid (RL and CG constraints are satisified)
 * :param oligo: The oligo ( a vector of 0/1s)
 */
{
    int atcg_counts[4] = {0, 0, 0, 0};

    for (int x=0; x<oligo.size(); x+=2) {
        int nuc = oligo[x] << 1 | oligo[x + 1];
        atcg_counts[nuc]++;
    }
    int at = atcg_counts[0] + atcg_counts[1];
    int cg =  atcg_counts[2] + atcg_counts[3];
    //  GC - content constraint
    if (at > 7 || cg > 7) { // The oligo length is 12. the at/cg contents needs to be [0.4,0.6] so between 5 to 7
        return 1;
    }

    int successive = 1;
    int prev = (oligo[0]<<1)|oligo[1];
    int current;
    for (int i=2; i < oligo.size(); i+=2) {
        current = (oligo[i]<<1) | oligo[i + 1];
        if (prev == current) {
            successive++;
            if (successive == 4) {
                return 1;
            }
        }
        else {
            successive = 1;
        }
        prev = current;
    }
    return 0;
}

// LDPC
int IGNORED_POSITIONS[] = {6, 10, 14, 18};

// Algo 5-6
int balance_parity_oligo(Vec& parity_oligo, int x)
/**
 * We balance the parity oligo not with the knuth balance technique but rather
 * by inserting errors. THIS MAY LEAD TO DECODING errors if we are unlucky.
 * (Specifically if two oligos used together for decoding have an error in the same index we may have an error).
 * :param parity_oligo: The oligo to balance.
*/
{
    int current_balance = 0;
    for (int i=0;i<parity_oligo.size();i+=2) {
        bool skip=false;
        for (int j=0;j<4;j++) {
            if (i == IGNORED_POSITIONS[j]) {
                skip=true;
            }
        }
        if (skip) {
            continue;
        }

        if (parity_oligo[i] == 0) {
            current_balance++;
        } else {
            current_balance--;
        }
    }

    for (int i=0;i<4;i++) {
        if (current_balance > 0) { // More 0s
            parity_oligo[IGNORED_POSITIONS[i]] = 1;
            current_balance -= 1;
        } else {
            parity_oligo[IGNORED_POSITIONS[i]] = 0;
            current_balance += 1;
        }
    }
    // Handle edge case of 8 AT/CG before bbic bits
    if (current_balance < -3 || current_balance > 3) {
        parity_oligo[IGNORED_POSITIONS[x]-2] = !parity_oligo[IGNORED_POSITIONS[x]-2];
    }
    return 0;
}

vector<Vec> parity = vector<Vec>();
// currently for every 8 data oligos we use 4 parity oligos
int encode_ldpc()
/**
 * Takes the encoded BBIC oligos and creates parity bits for them using LDPC algorithm.
 * The parity oligos satisfy the RL and CG constraints for the price of possible failures to decode.
 */
{
    Vec par_1;
    Vec par_2;
    Vec par_3;
    Vec par_4;
    for (int i = 0; i < bbic_enc_oligo_vec_gl.size(); i += 8) {
        int p = 8;
        if (bbic_enc_oligo_vec_gl.size() - i < 8) {
            p = bbic_enc_oligo_vec_gl.size() - i;
        }
        par_1 = {
            0, 0, 0, 0, 0, 0, 0, 0,
            0, 0, 0, 0, 0, 0, 0, 0,
            0, 0, 0, 0, 0, 0, 0, 0
        }; // K=8
        par_2 = {
            0, 0, 0, 0, 0, 0, 0, 0,
            0, 0, 0, 0, 0, 0, 0, 0,
            0, 0, 0, 0, 0, 0, 0, 0
        }; // K=8
        par_3 = {
            0, 0, 0, 0, 0, 0, 0, 0,
            0, 0, 0, 0, 0, 0, 0, 0,
            0, 0, 0, 0, 0, 0, 0, 0
        }; // K=8
        par_4 = {
            0, 0, 0, 0, 0, 0, 0, 0,
            0, 0, 0, 0, 0, 0, 0, 0,
            0, 0, 0, 0, 0, 0, 0, 0
        }; // K=8

        for (int j = 0; j < 2 * size_k_gl + 4; j++) {
            bool ignored = false;
            Vec word = {0, 0, 0, 0, 0, 0, 0, 0}; // K=8
            for (int k = 0; k < p; k++) {
                word[k] = bbic_enc_oligo_vec_gl[i + k].second[j];
            }
            Vec encoded_data = algorithm5_encode_simple(word, G_example);
            for (int k = 0; k < size(IGNORED_POSITIONS); k++) {
                if (j == IGNORED_POSITIONS[k]|| j == (IGNORED_POSITIONS[k]+1)) {
                    ignored = true;
                }
            }
            if (ignored) {
                par_1[j] = !par_1[j - 2];
                par_2[j] = !par_2[j - 2];
                par_3[j] = !par_3[j - 2];
                par_4[j] = !par_4[j - 2];
            } else {
                par_1[j] = encoded_data[8];
                par_2[j] = encoded_data[9];
                par_3[j] = encoded_data[10];
                par_4[j] = encoded_data[11];
            }
        }

        balance_parity_oligo(par_1,0);
        balance_parity_oligo(par_2,1);
        balance_parity_oligo(par_3,2);
        balance_parity_oligo(par_4,3);

        parity.emplace_back(par_1);
        parity.emplace_back(par_2);
        parity.emplace_back(par_3);
        parity.emplace_back(par_4);
    }

    return 0;
}


int decode_ldpc()
/**
 * Decodes the given LDPC encoded oligos. If a group of oligos and parity oligos fails to decode
 * (because of errors inserted in the oligos / parity bits) we ignore the
 * parity oligos (as our LDPC code (/ matrix) is systematic).
 */
{
    for (int i = 0; i < bbic_enc_oligo_vec_gl.size(); i += 8) {
        bool tmp[8][24] = {false};

        int p = 8;
        if (bbic_enc_oligo_vec_gl.size() - i < 8) {
            p = bbic_enc_oligo_vec_gl.size() - i;
        }

        for (int j = 0; j < 24; j++) {
            Vec word = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}; // K=8
            for (int k = 0; k < p; k++) {
                word[k] = bbic_enc_oligo_vec_gl[i + k].second[j];
            }
            for (int k = 0; k < 4; k++) {
                word[8 + k] = parity[k + 4 * i / 8][j];
            }
            for (int k = 0; k < 8; k++) {
                bool is_ignored = false;
                for (int r = 0; r < size(IGNORED_POSITIONS); r++) {
                    if (j == IGNORED_POSITIONS[r] || j == (IGNORED_POSITIONS[r] + 1)) {
                        is_ignored = true;
                    }
                }
                if (!is_ignored) {
                    Vec t = algorithm6_decode_simple(word, H_example, 5).first;

                    if (t.size() == 0) {
                        // Best effort decoding, if the problem is in the ldpc parity bits this is still fine.
                        cout << "Failed to decode" << endl;
                        tmp[k][j] = word[k];
                        continue;
                    }
                    tmp[k][j] = t[k];
                } else {
                    tmp[k][j] = word[k];
                }
            }
        }

        for (int j = 0; j < p; j++) {
            vector tmp_vec(24, false);
            for (int k = 0; k < 24; k++) {
                bbic_enc_oligo_vec_gl[i+j].second[k] = tmp[j][k];
            }
        }
    }
    return 0;
}



int main() {
    for (int t=0; t<100; t++) {
        data_vec_gl.clear();
        bbic_enc_oligo_vec_gl.clear();
        restored_vec_gl.clear();

        //data_vec_gl = {1,1,0,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1};
        fill_random_data_vec_gl(6400); //create random input data for test (must be even length)
        //restart algo param according to paper
        unsigned long data_size = data_vec_gl.size();
        restart_param(3, 10, 8);

        //check param
        // size_m >= because of the knuth balance index
        if (2 * size_m_gl >= size_k_gl || size_m_gl <= 2) {
            cout << "K and M is invalid" << endl;
            exit(1);
        }

        //print original data
        cout << endl << endl << endl << endl;
        cout << "ORIGINAL DATA" << endl;
        cout << endl;
        printNucVectorBool(data_vec_gl);
        fflush(stdout); // Force the output immediately



        encode_bbic();

        //print encode bbic data
        cout << endl << endl << endl << endl;
        cout << "ENCODE BBIC DATA" << endl;
        cout << endl;

        for (int i = 0; i < bbic_enc_oligo_vec_gl.size(); i++) {
            // printf("%d, ", i);

            printNucVectorBool(bbic_enc_oligo_vec_gl[i].second);
        }

        // knuth_balance();
        //print encode bbic data
        cout << endl << endl << endl << endl;
        cout << "KNUTH BALANCE DATA" << endl;
        cout << endl;
        for (int i = 0; i < bbic_enc_oligo_vec_gl.size(); i++) {
            // printf("%d, ", i);
            printNucVectorBool(bbic_enc_oligo_vec_gl[i].second);
        }

        // Verify oligos are balanced
        for (int i = 0; i < bbic_enc_oligo_vec_gl.size(); i++) {
            // printf("%d, ", i);

            printNucVectorBool(bbic_enc_oligo_vec_gl[i].second);
            if (verify_oligo(bbic_enc_oligo_vec_gl[i].second) != 0) {
                cout << "BAD OLIGO" << endl;
                return 1; // BAD OLIGO
            };
        }

        encode_ldpc();

        // Check if we have enough oligos to run this loop safely
        if (bbic_enc_oligo_vec_gl.size() > 32) {
            for (size_t x = 0; x < bbic_enc_oligo_vec_gl.size() - 32; x += 16) {

                // Safety check for the specific offsets used inside
                if (x + 7 < bbic_enc_oligo_vec_gl.size()) {

                    // Introduce errors
                    bbic_enc_oligo_vec_gl[x+1].second[2] = !bbic_enc_oligo_vec_gl[x+1].second[2];

                    // Typo Fixed: You were flipping x+2 based on x+1. Usually you flip based on itself.
                    bbic_enc_oligo_vec_gl[x+2].second[5] = !bbic_enc_oligo_vec_gl[x+2].second[5];

                    bbic_enc_oligo_vec_gl[x+7].second[9] = !bbic_enc_oligo_vec_gl[x+7].second[9];
                }
            }
        }

        decode_ldpc();

        // decode_knuth_balance();
        //print encode bbic data
        // cout << endl << endl << endl << endl;
        // cout << "KNUTH BALANCE DATA AFTER DECODE" << endl;
        // cout << endl;
        //
        // for (int i = 0; i < bbic_enc_oligo_vec_gl.size(); i++) {
        //     printNucVectorBool(bbic_enc_oligo_vec_gl[i].second);
        // }

        decode_bbic();

        cout << "----------------" << endl;
        cout << "size_DATA: " << data_vec_gl.size() << endl;
        cout << "size_beta: " << size_beta_gl << endl;
        cout << "K: " << size_k_gl << endl;
        cout << "m: " << size_m_gl << endl;


        //print restored data
        cout << endl << endl << endl << endl;
        cout << "RESTORED DATA" << endl;
        cout << endl;
        printNucVectorBool(restored_vec_gl);
        fflush(stdout); // Force the output immediately

        cout << "----------EQ?----------" << endl;
        if (!areEqual(data_vec_gl, restored_vec_gl, true)) {
            return 1;
        }

        areEqual(data_vec_gl, restored_vec_gl, true);
        parity.clear();
        // exit(0);
    }


}