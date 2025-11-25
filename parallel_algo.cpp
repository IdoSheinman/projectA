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
    //    kb_idxs[0] = 2; // AC
    //    kb_idxs[1] = 3; // AG
    //    kb_idxs[2] = 6; // TC
    //    kb_idxs[3] = 7; // TG

    //   kb_idxs[4] = 8; // CA
    //    kb_idxs[5] = 9; // CT
    //    kb_idxs[6] = 12; // GA
    //    kb_idxs[7] = 13; // GT

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
            it->second[l - idx_s] = data_vec_gl[l];
        } else {
            it->second[l - idx_s] = add_dummy_bits();
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


int  encode_bbic() {
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
            vector<bool> oligo(oligo_data_size);
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
/**
 * Decodes the encoded bbic oligos back into data.
 */
{
    restored_vec_gl = vector<bool>();
    int oligo_num = 0;
    vector<bool> Q_data = vector<bool>();
    while (oligo_num < bbic_enc_oligo_vec_gl.size()) {
        for (int ol_n = oligo_num; ol_n < oligo_num + parallel_number_gl; ol_n++){
            if (ol_n >= bbic_enc_oligo_vec_gl.size()) break;

            // take_data_without_Q
            int first_q = 2 * size_m_gl + 1;
            for (int l = 0; l < 2 * size_k_gl; l++) {
                if (l < first_q || ((l - first_q) % (2 * size_m_gl) != 0)) {
                    restored_vec_gl.insert(restored_vec_gl.end(), 1, bbic_enc_oligo_vec_gl[ol_n].second[l]);//TODO: ERR
                }
            }




            // add_Q_data

            for (int l = first_q; l <= bbic_enc_oligo_vec_gl[ol_n].second.size(); l += 2 * size_m_gl) {
                //TODO: ERR
                if (!(is_condition_1_satisfied(bbic_enc_oligo_vec_gl[ol_n].second, l) && is_condition_2_satisfied(
                bbic_enc_oligo_vec_gl[ol_n].second, l))) {
                    Q_data.insert(Q_data.end(), 1, bbic_enc_oligo_vec_gl[ol_n].second[l]);
                }
            }
        }
        restored_vec_gl.insert(restored_vec_gl.end(), Q_data.begin(), Q_data.end());
        Q_data = vector<bool>();
        oligo_num = oligo_num + parallel_number_gl;
    }
    return 0;
}


int main() {
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
        printNucVectorBool(bbic_enc_oligo_vec_gl[i].second);
    }


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

    exit(0);



}