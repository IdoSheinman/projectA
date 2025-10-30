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

//GLOBAL VAR
vector<bool> data_vec;
vector<pair<int, vector<bool> > > result_vec_bic;
vector<Vec> parity = vector<Vec>();
vector<bool> decoded_vec_bic;
vector<vector<bool> > decoded_vec_bic_ldpc;
int i_s, i_e, size_m, size_k, size_beta;
unsigned long size_n;
bool used_all_data;
int dummy_bit_counter;
int kb_idxs[8];

bool add_dummy_bits();

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


void printNucVectorBool(Vec vec) {
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

int restart_param(int m, int k)
/**
 * Initialize the input parameters.
 * :param m: longest allowed RL value.
 * :param k: As appears in the paper.
 */
{
    size_m = m;
    size_k = k;
    size_beta = (size_k + size_m - 1) / size_m - 1;
    i_s = 0;
    i_e = 0;
    used_all_data = false;
    dummy_bit_counter = 0;
    // BBIC
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

int addEmptyItems(std::vector<bool> &result, int m)
/**
 * Fills result with zeros up to its size.
 */
{
    for (int i = 2 * m + 1 - 4; i <= result.size(); i += 2 * m) {
        result.insert(result.begin() + i, false);
    }
    return 0;
}

int copy_data_without_Q(int ol_i)
/**
 * this writes the raw data into the oligo, while not writing anything
 * for the Q indexes.
 * :param ol_i: the index of the oligo to fill.
 */
{
    auto it = result_vec_bic.begin();
    for (; it != result_vec_bic.end(); it++) {
        if (it->first == ol_i) {
            break;
        }
    }

    if (it == result_vec_bic.end()) {
        return 1;
    }


    for (int l = i_s; l <= i_e; l++) {
        if (l < data_vec.size()) {
            it->second[l - i_s] = data_vec[l];
        } else {
            it->second[l - i_s] = add_dummy_bits();
            dummy_bit_counter++;
        }
    }

    return (addEmptyItems((it->second), size_m));
}


bool is_condition_1_satisfied(vector<bool> oli, int position)
/**
 * Check for condition 1 as appears in the paper in page 768.
 * :param oli: The index of the oligo to check.
 * :param position: The index of the position to check (where the Q value is inserted.
 */
{
    int first_upper_pos = position - (2 * size_m + 1);
    int first_lower_pos = position - (2 * size_m);

    bool curr_upper_bit;
    bool curr_lower_bit;
    int current_streak = 1;
    for (int i = 2; i < 2 * (2 * size_m); i += 2) {
        if (first_lower_pos + i >= oli.size()) {
            return false;
        }
        if (i == 2 * size_m) {
            continue;
        }
        curr_upper_bit = oli[first_upper_pos + i];
        curr_lower_bit = oli[first_lower_pos + i];
        if (i == 2 * size_m + 2) {
            if (curr_upper_bit == oli[first_upper_pos + i - 4] &&
                curr_lower_bit == oli[first_lower_pos + i - 4]) {
                current_streak += 1;
            } else {
                current_streak = 1;
            }
        } else {
            if (curr_upper_bit == oli[first_upper_pos + i - 2] &&
                curr_lower_bit == oli[first_lower_pos + i - 2]) {
                current_streak += 1;
            } else {
                current_streak = 1;
            }
        }
        if (current_streak == size_m) {
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


int add_bbic_bits(int ol_i)
/**
 * This function write into the bbic positions their correct values ( either stops a sequence of repeating ATCG or data )
 * :param ol_i: The  oligo index
 */
{
    int d_i = 0;
    auto it = result_vec_bic.begin();
    for (; it != result_vec_bic.end(); it++) {
        if (it->first == ol_i) {
            break;
        }
    }

    if (it == result_vec_bic.end()) {
        // This shouldnt happen, but for hw purposes not sure return 1 is ideal.
        return 1;
    }

    for (int l = 2 * size_m + 1; l < it->second.size(); l += 2 * size_m) {
        // If need to stop sequence.
        if (is_condition_1_satisfied(it->second, l) && is_condition_2_satisfied(it->second, l)) {
            // TODO: Shouldnt this be l-1
            it->second[l] = !(it->second[l - 2]);
        }
        // Else can fill data
        else {
            d_i++;
            if (i_e + d_i < data_vec.size()) {
                it->second[l] = data_vec[i_e + d_i];
                if (i_e + d_i == data_vec.size() - 1) {
                    used_all_data = true;
                }
            }
            // We want to fill data but there is no more data to fill
            else {
                cout << "ADD DUMMY BITS" << endl;
                it->second[l] = add_dummy_bits(); //need to see what val at the end
                dummy_bit_counter++;
            }
        }
    }
    return d_i;;
}

bool add_dummy_bits()
/**
 * This gives us dummy bits to read if there is no more data, in a way that doesnt create long running A(/C/T/G) only sequences.
 */
{
    used_all_data = true;
    bool dummy_arr[] = {false, false, false, true, true, false, true, true};
    return dummy_arr[dummy_bit_counter % 8];
}

int decode_Q_data(int ol_i)
/**
 * Retrieves the data encoded in a given bbic oligos Q bits, and inserts them into decoded_vec_bic.
 * :param ol_i: The oligo index.
 */
{
    unsigned int first_q = 2 * size_m + 1;
    auto it = result_vec_bic.begin();
    for (; it != result_vec_bic.end(); it++) {
        if (it->first == ol_i) {
            break;
        }
    }

    if (it == result_vec_bic.end()) {
        return 1;
    }
    for (int l = first_q; l <= it->second.size(); l += 2 * size_m) {
        if (!(is_condition_1_satisfied(it->second, l) && is_condition_2_satisfied(it->second, l))) {
            decoded_vec_bic.insert(decoded_vec_bic.end(), 1, it->second[l]);
        }
    }
    return 0;
}


int take_data_without_Q(int ol_i)
/**
 * Retrieves the data encoded in a given bbic oligos non-Q bits, and inserts them into decoded_vec_bic.
 * :param ol_i: The oligo index.
 */
{
    auto it = result_vec_bic.begin();
    for (; it != result_vec_bic.end(); ++it) {
        if (it->first == ol_i) {
            break;
        }
    }

    if (it == result_vec_bic.end()) {
        return 1;
    }
    unsigned int first_q = 2 * size_m + 1;
    for (int l = 0; l < 2 * size_k; l++) {
        if (l < first_q || ((l - first_q) % (2 * size_m) != 0)) {
            decoded_vec_bic.insert(decoded_vec_bic.end(), 1, it->second[l]);
        }
    }

    return 0;
}

void fill_random_data()
/**
 * This function creates random data for us to encode and decode to verify our success.
 */
{
    // Constants for 0 to 1 MB size in bits
    const size_t MAX_BITS = 128 * 8 * 8; // 1 KB in bits

    // Seed random number generator
    random_device rd; // Random device for seed
    mt19937 gen(rd()); // Mersenne Twister RNG
    uniform_int_distribution<size_t> size_dist(0, MAX_BITS); // Size between 0 and 1MB
    uniform_int_distribution<int> bit_dist(0, 1); // Values 0 or 1

    // Generate a random size
    size_t test_size = size_dist(gen) * 8;

    // Resize the vector
    data_vec.resize(test_size);

    // Fill with random values
    for (size_t i = 0; i < test_size; ++i) {
        data_vec[i] = static_cast<bool>(bit_dist(gen));
    }
}



int find_best_idx(int idx)
/**
 * This function finds the best index, for a given idx,
 * to knuth balance it out of the possible values.
 * :param idx: The index of the oligo to find the best index for.
 */
{
    vector<bool> olig = result_vec_bic[idx].second;
    int scores[8] = {0};
    for (int i = 0; i < 8; i++) {
        vector<bool> tmp_olig;
        copy(olig.begin(), olig.end(), back_inserter(tmp_olig));
        for (int j = 0; j <= real_kb_idxs(kb_idxs[i]); j += 2) {
            tmp_olig[j].flip();
        }
        int counter = 0;
        for (int j = 0; j < tmp_olig.size(); j += 2) {
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


void knuth_balance(int idx)
/**
 * Flips bits according to the knuth balance technique and inserts the knuth balance index in the beggining of the
 * given oligo.
 * :param idx: The index of the oligo to balance.
 */
{
    vector<bool> olig = result_vec_bic[idx].second;
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
    for (int i = 0; i <= real_kb_idxs(idx_to_balance); i += 2) {
        result_vec_bic[idx].second[i] = !result_vec_bic[idx].second[i];
    }

    // 4 index bits, depend on oligo size (this should work with 20)
    result_vec_bic[idx].second.emplace(result_vec_bic[idx].second.begin(), idx_to_balance & 1);
    result_vec_bic[idx].second.emplace(result_vec_bic[idx].second.begin(), ((idx_to_balance & 2) >> 1) & 1);
    result_vec_bic[idx].second.emplace(result_vec_bic[idx].second.begin(), ((idx_to_balance & 4) >> 2) & 1);
    result_vec_bic[idx].second.emplace(result_vec_bic[idx].second.begin(), ((idx_to_balance & 8) >> 3) & 1);
}

int encode_bbic()
/**
 * Encodes the data into valid bbic oligos.
 */
{
    int oligo_num = 0;
    int rest;
    int oligo_data_size = 2 * size_k - size_beta;
    while (!used_all_data) {
        i_e = i_s + oligo_data_size - 1;
        vector<bool> oligo(oligo_data_size);
        // Add the pair (i, oligo) to result_vec
        result_vec_bic.emplace_back(oligo_num, oligo);
        if (!copy_data_without_Q(oligo_num)) {
            knuth_balance(oligo_num);
            int added = add_bbic_bits(oligo_num); //add q bits
            i_s = i_e + 1 + added;
        } else {
            return 1;
        }
        oligo_num++;
        rest = signed(data_vec.size() - i_s);
        if (rest <= 0) {
            used_all_data = true;
        }
    }

    return 0;
}

int decode_bbic()
/**
 * Decodes the encoded bbic oligos back into data.
 */
{
    decoded_vec_bic = vector<bool>();
    for (int i = 0; i < decoded_vec_bic_ldpc.size(); i++) {
        int kbt_idx =
                (decoded_vec_bic_ldpc[i][0] << 3) |
                (decoded_vec_bic_ldpc[i][1] << 2) |
                (decoded_vec_bic_ldpc[i][2] << 1) |
                (decoded_vec_bic_ldpc[i][3] << 0);

        int first_q = 2 * size_m + 1;


        // add_Q_data
        vector<bool> Q_data = vector<bool>();
        for (int l = first_q; l <= decoded_vec_bic_ldpc[i].size(); l += 2 * size_m) {
            if (!(is_condition_1_satisfied(decoded_vec_bic_ldpc[i], l) && is_condition_2_satisfied(
                      decoded_vec_bic_ldpc[i], l))) {
                Q_data.insert(Q_data.end(), 1, decoded_vec_bic_ldpc[i][l]);
            }
        }

        for (int j = 4; j <= real_kb_idxs(kbt_idx) + 4; j += 2) {
            decoded_vec_bic_ldpc[i][j] = !decoded_vec_bic_ldpc[i][j];
        }

        // take_data_without_Q
        for (int l = 4; l < 2 * size_k + 4; l++) {
            if (l < first_q || ((l - first_q) % (2 * size_m) != 0)) {
                decoded_vec_bic.insert(decoded_vec_bic.end(), 1, decoded_vec_bic_ldpc[i][l]);
            }
        }


        // add_Q_data
        decoded_vec_bic.insert(decoded_vec_bic.end(), Q_data.begin(), Q_data.end());
    }
    return 0;
}


// Q data
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
    for (int i = 0; i < result_vec_bic.size(); i += 8) {
        int p = 8;
        if (result_vec_bic.size() - i < 8) {
            p = result_vec_bic.size() - i;
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

        for (int j = 0; j < 2 * size_k + 4; j++) {
            bool ignored = false;
            Vec word = {0, 0, 0, 0, 0, 0, 0, 0}; // K=8
            for (int k = 0; k < p; k++) {
                word[k] = result_vec_bic[i + k].second[j];
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
    for (int i = 0; i < result_vec_bic.size(); i += 8) {
        bool tmp[8][24] = {false};

        int p = 8;
        if (result_vec_bic.size() - i < 8) {
            p = result_vec_bic.size() - i;
        }

        for (int j = 0; j < 24; j++) {
            Vec word = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}; // K=8
            for (int k = 0; k < p; k++) {
                word[k] = result_vec_bic[i + k].second[j];
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
            vector<bool> tmp_vec(24, false);
            for (int k = 0; k < 24; k++) {
                tmp_vec[k] = tmp[j][k];
            }
            decoded_vec_bic_ldpc.emplace_back(tmp_vec);
        }
    }
    return 0;
}

int main() {
    for (int i=0; i<50; i++) {
        result_vec_bic.clear();
        data_vec.clear();
        result_vec_bic.clear();
        parity.clear();
        decoded_vec_bic.clear();
        decoded_vec_bic_ldpc.clear();

        fill_random_data(); //create random input data for test (must be even length)
        //restart algo param according to paper
        size_n = data_vec.size();
        restart_param(3, 10);
        decoded_vec_bic = vector<bool>();

        //check param
        // size_m >= because of the knuth balance index
        if (2 * size_m >= size_k || size_m <= 2) {
            cout << "K and M is invalid" << endl;
            cin;
        }


        if (encode_bbic() != 0) {
            return 1;
        }

        cout << endl << endl << endl << endl;
        cout << "ENCODE DATA" << endl;
        cout << endl << endl << endl << endl;

        for (int i = 0; i < result_vec_bic.size(); i++) {
            // printNucVectorBool(result_vec_bic[i].second);
            if (verify_oligo(result_vec_bic[i].second) != 0) {
                cout << "BAD OLIGO" << endl;
                return 1; // BAD OLIGO
            };
        }



        encode_ldpc();

        cout << endl << endl << endl << endl;
        cout << "LDPC PARITY OLIGOS:" << endl;
        cout << endl << endl << endl << endl;

        for (int i = 0; i < size(parity); i++) {
            vector<bool> v = {
                0, 0, 0, 0, 0, 0, 0, 0,
                0, 0, 0, 0, 0, 0, 0, 0,
                0, 0, 0, 0, 0, 0, 0, 0
            }; // K=8

            for (int j = 0; j < 24; j++) {
                v[j] = parity[i][j];
            }
            // printNucVectorBool(v);
            if (verify_oligo(v) != 0) {
                cout << "BAD VERIFY" << endl;
                return 1; // BAD VERIFY
            };

        }

        // Insert errors for testing ldpc error correction
        // Note: There are errors that appear statistically due to the
        //       CG constraint on the parity bits. We insert them
        //       in a position dependent on the data (random data we generate)
        //       This can cause the error to appear in the same column as an
        //       error we insert here, which would result in wrong data
        //       being recovered.
        for (int x=0; x<result_vec_bic.size(); x+=16) {
            result_vec_bic[x+1].second[2] = !result_vec_bic[x+1].second[2];
            result_vec_bic[x+2].second[5] = !result_vec_bic[x+1].second[5];
            result_vec_bic[x+7].second[9] = !result_vec_bic[x+7].second[9];
            result_vec_bic[x+4].second[12] = !result_vec_bic[x+4].second[12];
            result_vec_bic[x+5].second[8] = !result_vec_bic[x+5].second[8];
        }

        decode_ldpc();

        decode_bbic();

        cout << "----------------" << endl;
        cout << "size_DATA: " << size_n << endl;
        cout << "size_beta: " << size_beta << endl;
        cout << "K: " << size_k << endl;
        cout << "m: " << size_m << endl;


        cout << "----------EQ?----------" << endl;
        if (!areEqual(data_vec, decoded_vec_bic, true)) {
            return 1;
        }

        areEqual(data_vec, decoded_vec_bic, true);
    }
    return 0;
}










//garbage
//
// int pad_till_success(vector<bool> data) {
//     int counter = 0;
//     while (!areEqual(data_vec, decoded_vec_bic)) {
//         data_vec = data;
//         size_n = data_vec.size();
//         for (int i = 0; i < counter; ++i) {
//             data_vec.insert(data_vec.end(), 1, 0);
//         }
//
//         decoded_vec_bic = vector<bool>();
//         restart_param(3, 10);
//         if (2 * size_m >= size_k) {
//             cout << "K and M is invalid" << endl;
//             cin;
//         }
//         encode_bic();
//         decode_bic();
//         counter += 1;
//     }
//     return counter;
// }


// int encode_bic()
// /**
//  * Fills result_vec_bic withc valid bbic oligos.
//  */
// {
//     int oligo_num = 1;
//     int rest;
//
//     while (!used_all_data) {
//         i_e = i_s + 2 * size_k - size_beta - 1;
//         vector<bool> oligo(2 * size_k - size_beta);
//         // Add the pair (i, oligo) to result_vec
//         result_vec_bic.emplace_back(oligo_num, oligo);
//
//         if (!copy_data_without_Q(oligo_num)) {
//             int added = add_bbic_bits(oligo_num);
//             i_s = i_e + 1 + added;
//         } else {
//             return 1;
//         }
//         oligo_num++;
//         rest = signed(data_vec.size() - i_s);
//         if (rest <= 0) {
//             used_all_data = true;
//         }
//     }
//
//     return 0;
// }
//
// int decode_bic()
// /**
//  * Retrieves the data encoded results_vec_bicc.
//  * :param ol_i: The oligo index.
//  */
// {
//     for (auto oligo_vec: result_vec_bic) {
//         take_data_without_Q(oligo_vec.first);
//         decode_Q_data(oligo_vec.first);
//     }
//     return 0;
// }
