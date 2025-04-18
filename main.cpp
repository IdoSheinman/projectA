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
vector<vector<bool> > mem;
int current_vec = 0;

vector<bool> data_vec;
vector<pair<int, vector<bool> > > result_vec_bic;
vector<Vec> parity = vector<Vec>();
vector<bool> decoded_vec_bic;
vector<vector<bool>> decoded_vec_bic_ldpc;
int i_s, i_e, size_m, size_k, size_beta;
unsigned long size_n;
bool used_all_data;
int counter;
int kb_idxs[8];

bool add_dummy_bits();

void printVectorBool(const vector<bool> &vec) {
    for (bool val: vec) {
        cout << val << " ";
    }
    cout << endl;
}

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
    /* if (vec1.size() != vec2.size()) {
         cout << "The vectors have different sizes, so differences cannot be computed." << endl;
         return;
     }*/
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

int restart_param(int m, int k) {
    size_m = m;
    size_k = k;
    size_beta = (size_k + size_m - 1) / size_m - 1;
    i_s = 0;
    i_e = 0;
    used_all_data = false;
    counter = 0;
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

int addEmptyItems(std::vector<bool> &result, int m) {
    for (int i = 2 * m + 1-4; i <= result.size(); i += 2 * m) {
        result.insert(result.begin() + i, false);
    }
    //result.resize(2*size_k);
    return 0;
}

int copy_data_without_Q(int ol_i) {
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
            counter++;
        }
    }

    return (addEmptyItems((it->second), size_m));
}


bool is_condition_1_satisfied(vector<bool> oli, int position) {
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

bool is_condition_2_satisfied(vector<bool> oli, int position) {
    return (oli[position - 1] == oli[position - 3]);
}


int add_bbic_bit(int ol_i) {
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
        if (is_condition_1_satisfied(it->second, l) && is_condition_2_satisfied(it->second, l)) {
            // TODO: Shouldnt this be l-1?
            it->second[l] = !(it->second[l - 2]);
        } else {
            d_i++;
            //cout << "ADDING Q BIT IN: " << l << " from: " << i_e+d_i << endl;
            if (i_e + d_i < data_vec.size()) {
                it->second[l] = data_vec[i_e + d_i];
                if (i_e + d_i == data_vec.size() - 1) {
                    used_all_data = true;
                }
            } //if the last???
            else {
                cout << "ERROR: Reached unreachable code" << endl;
                it->second[l] = add_dummy_bits(); //need to see what val at the end
                counter++;
            }
        }
    }
    return d_i;;
}

bool add_dummy_bits() {
    used_all_data = true;
    bool dummy_arr[] = {false, false, false, true, true, false, true, true};
    //data_vec.insert(data_vec.end(),1,dummy_arr[counter % 8]);
    //counter++;
    //size_n = data_vec.size();
    return dummy_arr[counter % 8];
}

int encode_bic() {
    int oligo_num = 1;
    int rest;

    while (!used_all_data) {
        i_e = i_s + 2 * size_k - size_beta - 1;
        vector<bool> oligo(2 * size_k - size_beta);
        // Add the pair (i, oligo) to result_vec
        result_vec_bic.emplace_back(oligo_num, oligo);

        if (!copy_data_without_Q(oligo_num)) {
            int added = add_bbic_bit(oligo_num);
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

int add_Q_data(int ol_i) {
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


int take_data_without_Q(int ol_i) {
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
        //TODO: Skip Q
        if (l < first_q || ((l - first_q) % (2 * size_m) != 0)) {
            decoded_vec_bic.insert(decoded_vec_bic.end(), 1, it->second[l]);
        }
    }

    return 0;
}

int decode_bic() {
    for (auto oligo_vec: result_vec_bic) {
        take_data_without_Q(oligo_vec.first);
        add_Q_data(oligo_vec.first);
    }
    return 0;
}

void fill_random_data() {
    // Constants for 0 to 1 MB size in bits
    const size_t MAX_BITS = 128 * 8 * 8; // 1 KB in bits
    // const size_t MAX_BITS = 8 * 8 * 8; // 1 KB in bits

    // Seed random number generator
    random_device rd; // Random device for seed
    mt19937 gen(rd()); // Mersenne Twister RNG
    uniform_int_distribution<size_t> size_dist(0, MAX_BITS); // Size between 0 and 1MB
    uniform_int_distribution<int> bit_dist(0, 1); // Values 0 or 1

    // Generate a random size
    size_t random_size = 128*8*8*3;//size_dist(gen);

    // Resize the vector
    data_vec.resize(random_size);

    // Fill with random values
    for (size_t i = 0; i < random_size; ++i) {
        data_vec[i] = static_cast<bool>(bit_dist(gen));
    }
}

int pad_till_success(vector<bool> data) {
    int counter = 0;
    while (!areEqual(data_vec, decoded_vec_bic)) {
        data_vec = data;
        size_n = data_vec.size();
        for (int i = 0; i < counter; ++i) {
            data_vec.insert(data_vec.end(), 1, 0);
        }

        decoded_vec_bic = vector<bool>();
        restart_param(3, 10);
        if (2 * size_m >= size_k) {
            cout << "K and M is invalid" << endl;
            cin;
        }
        encode_bic();
        decode_bic();
        counter += 1;
    }
    return counter;
}


int find_best_idx(int idx) {
    vector<bool> olig = result_vec_bic[idx].second;
    int scores[8] = {0};
    for (int i=0; i<8; i++) {
        vector<bool> tmp_olig;
        copy(olig.begin(), olig.end(), back_inserter(tmp_olig));
        for (int j = 0; j < kb_idxs[i]; j+=1) {
            tmp_olig[j].flip();
        }
        int counter = 0;
        for (int i = 0; i < tmp_olig.size(); i+=2) {
            if (tmp_olig[i] == 0) {
                counter += 1;
            } else {
                counter -= 1;
            }
        }
        scores[i] = abs(counter);
    }


    int index = 0;

    for(int i = 0; i < 8; i++)
    {
        if(scores[i] < scores[index])
            index = i;
    }
    return kb_idxs[index];
}


void knuth_balance(int idx) {
    vector<bool> olig = result_vec_bic[idx].second;
    int zeros = 0;
    int ones = 0;
    for (int i = 0; i < olig.size(); i+=2) {
        if (olig[i] == 0) {
            zeros += 1;
        }
        if (olig[i] == 1) {
            ones += 1;
        }
    }

    int too_many_ones = (ones - zeros)/2;
    int tmp = 0;
    int idx_to_balance = 0;
    for (idx_to_balance = 0; idx_to_balance < olig.size(); idx_to_balance+=2) {
        if (olig[idx_to_balance] == 1) {
            tmp += 1;
        }
        else {
            tmp -= 1;
        }

        if (tmp == too_many_ones) {
            break;
        }
    }
    idx_to_balance = find_best_idx(idx);
    for (int i = 0; i <= idx_to_balance; i += 2) {
        result_vec_bic[idx].second[i] = !result_vec_bic[idx].second[i];
    }

    result_vec_bic[idx].second.emplace(result_vec_bic[idx].second.begin(), idx_to_balance & 1);
    result_vec_bic[idx].second.emplace(result_vec_bic[idx].second.begin(), ((idx_to_balance & 2)>>1)&1);
    result_vec_bic[idx].second.emplace(result_vec_bic[idx].second.begin(), ((idx_to_balance & 4)>>2)&1);
    result_vec_bic[idx].second.emplace(result_vec_bic[idx].second.begin(), ((idx_to_balance & 8)>>3)&1);
    // result_vec_bic[idx].second.emplace(result_vec_bic[idx].second.begin(), ((idx_to_balance & 16)>>4)&1);
    // result_vec_bic[idx].second.emplace(result_vec_bic[idx].second.begin(), ((idx_to_balance & 32)>>5)&1);
}

// WIP
int encode_bbic() {
    int oligo_num = 0;
    int rest;
    int oligo_size = 2 * size_k - size_beta;
    while (!used_all_data) {
        i_e = i_s + oligo_size - 1;
        vector<bool> oligo(oligo_size);
        // Add the pair (i, oligo) to result_vec
        result_vec_bic.emplace_back(oligo_num, oligo);
        if (!copy_data_without_Q(oligo_num)) {
            knuth_balance(oligo_num);
            int added = add_bbic_bit(oligo_num);
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

int decode_bbic() {
    decoded_vec_bic = vector<bool>();
    // for (auto oligo_vec: result_vec_bic) {
        for (int i=0; i<result_vec_bic.size(); i++) {
            int kbt_idx =
                  // (result_vec_bic[i].second[0] << 5) |
                  // (result_vec_bic[i].second[1] << 4) |
                  (result_vec_bic[i].second[0] << 3) |
                  (result_vec_bic[i].second[1] << 2) |
                  (result_vec_bic[i].second[2] << 1) |
                  (result_vec_bic[i].second[3] << 0);

            int ol_i = result_vec_bic[i].first;

            // take_data_without_Q(result_vec_bic[i].first);
            auto it = result_vec_bic.begin();
            for (; it != result_vec_bic.end(); ++it) {
                if (it->first == ol_i) {
                    break;
                }
            }

            if (it == result_vec_bic.end()) {
                return 1;
            }
            int first_q = 2 * size_m + 1;


            // add_Q_data(result_vec_bic[i].first);
            vector<bool> Q_data = vector<bool>();
            for (int l = first_q; l <= it->second.size(); l += 2 * size_m) {
                if (!(is_condition_1_satisfied(it->second, l) && is_condition_2_satisfied(it->second, l))) {
                    Q_data.insert(Q_data.end(), 1, it->second[l]);
                }
            }

            //result_vec_bic[i].second.erase(result_vec_bic[i].second.begin(),result_vec_bic[i].second.begin()+6);
            for (int j = 4; j <= kbt_idx+4; j += 2) {
                result_vec_bic[i].second[j] = !result_vec_bic[i].second[j];
            }

            // take_data_without_Q(result_vec_bic[i].first);
            for (int l = 4; l < 2 * size_k + 4; l++) {
                //TODO: Skip Q
                if (l < first_q || ((l - first_q) % (2 * size_m) != 0)) {
                    decoded_vec_bic.insert(decoded_vec_bic.end(), 1, it->second[l]);
                }
            }



            // add_Q_data(result_vec_bic[i].first);
            decoded_vec_bic.insert(decoded_vec_bic.end(), Q_data.begin(), Q_data.end());
        }
    return 0;
}




int IGNORED_POSITIONS[] = {6, 10, 14, 18};
// Algo 5-6
int encode_ldpc() {
    Vec par_1;
    Vec par_2;
    Vec par_3;
    Vec par_4;
    for (int i = 0; i < result_vec_bic.size(); i += 8) {
        int p = 8;
        if (result_vec_bic.size() - i < 8) {
            p = result_vec_bic.size() - i;
        }

        par_1 = {0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0,
        0, 0}; // K=8
        par_2 = {0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0,
        0, 0}; // K=8
        par_3 = {0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0,
        0, 0}; // K=8
        par_4 = {0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0,
        0, 0}; // K=8

        for (int j = 0; j < 24; j++) {
            bool ignored = false;
            Vec word = {0, 0, 0, 0, 0, 0, 0, 0}; // K=8
            for (int k = 0; k < p; k++) {
               word[k] = result_vec_bic[i+k].second[j];
            }
            Vec encoded_data = algorithm5_encode_simple(word, G_example);
            for (int k = 0; k < size(IGNORED_POSITIONS); k++) {
                if (j == IGNORED_POSITIONS[k]) {
                    ignored = true;
                }
            }
            if (ignored) {
                par_1[j] = !par_1[j-2];
                par_2[j] = !par_2[j-2];
                par_3[j] = !par_3[j-2];
                par_4[j] = !par_4[j-2];
            } else {
                par_1[j]= encoded_data[8];
                par_2[j] = encoded_data[9];
                par_3[j] = encoded_data[10];
                par_4[j] = encoded_data[11];
            }
        }

        parity.emplace_back(par_1);
        parity.emplace_back(par_2);
        parity.emplace_back(par_3);
        parity.emplace_back(par_4);
    }

    return 0;
}

int decode_bbic_ldpc() {
    decoded_vec_bic = vector<bool>();
    for (int i=0; i<decoded_vec_bic_ldpc.size(); i++) {
        int kbt_idx =
              (decoded_vec_bic_ldpc[i][0] << 5) +
              (decoded_vec_bic_ldpc[i][1] << 4) +
              (decoded_vec_bic_ldpc[i][2] << 3) +
              (decoded_vec_bic_ldpc[i][3] << 2) +
              (decoded_vec_bic_ldpc[i][4] << 1) +
              (decoded_vec_bic_ldpc[i][5] << 0);
        for (int j = 6; j <= 6 + kbt_idx; j += 2) {
            decoded_vec_bic_ldpc[i][j] = !decoded_vec_bic_ldpc[i][j];
        }
        decoded_vec_bic_ldpc[i].erase(decoded_vec_bic_ldpc[i].begin(),decoded_vec_bic_ldpc[i].begin() + 6);

        for (int k = 0; k < result_vec_bic.size(); k++) {
            if (result_vec_bic[k].first == i) {
                result_vec_bic[k].second = decoded_vec_bic_ldpc[i];
                break;
            }
        }
        take_data_without_Q(result_vec_bic[i].first);
        add_Q_data(result_vec_bic[i].first);
    }
    return 0;
}


int decode_ldpc() {
    for (int i = 0; i < result_vec_bic.size(); i += 8) {
        bool tmp[8][24] = {false} ;

        int p = 8;
        if (result_vec_bic.size() - i < 8) {
            p = result_vec_bic.size() - i;
        }
        for (int j = 0; j < 24; j++) {
            Vec word = {0, 0, 0, 0, 0, 0, 0, 0,    0,0,0,0}; // K=8
            for (int k = 0; k < p; k++) {
                word[k] = result_vec_bic[i+k].second[j];
            }
            for (int k = 0; k < 4; k++) {
                word[8 + k] = parity[k+4*i/8][j];
            }
            Vec t = algorithm6_decode_simple(word, H_example, 5).first;

            if (t.size() == 0) {
                cout << "Failed to decode" << endl;
                return 0;
            }
            for (int k = 0; k < 8; k++) {
                bool is_ignored = false;
                for (int pos : IGNORED_POSITIONS) {
                    if (pos == j) {
                        is_ignored = true;
                        break;
                    }
                }
                if (!is_ignored) {
                    tmp[k][j] = t[k];
                } else {
                    tmp[k][j] = word[k];
                }
            }
            continue2:;
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
    //data_vec = vector<bool>(36, false);
    // data_vec.insert(data_vec.begin(), 1, true); 0,0, 0 ,1, 1 ,0 ,1, 1 ,0,0, 0 ,1, 1 ,0 ,1, 1,0,0, 0 ,1, 1 ,0 ,1, 1,0,0, 0 ,1, 1 ,0 ,1, 1,0,0, 0 ,1, 1 ,0 ,

    //0 0 0 1 0 1 0 1 0 1 0 0 0 1 0 1 0 0 1 1 1 0 0 0 0 0 1 0 0 1 0 0 0 0 0 0 1 1 0
    //data_vec = {0, 0, 0, 1, 0, 1, 0, 1, 0, 1, 0, 0, 0, 1, 0, 1, 0, 0, 1, 1, 1, 0, 0, 0, 0};
    //data_vec = {0, 0, 1, 1, 0, 1, 1, 1, 0, 1, 1, 0, 0, 1, 1, 0, 1, 0, 0, 0, 1, 1, 1, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0, 1, 0, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 1, 1, 1, 1, 0, 0, 1, 0, 0, 1, 1, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0, 1, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 1, 0, 1, 0, 0, 1, 0, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 1, 1, 1, 0, 0, 1, 1, 1, 0, 0, 0, 1, 1, 1, 1, 0, 0, 1, 0, 0, 1, 1, 1, 0, 1, 0, 0, 0, 1, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 1, 0, 1, 0, 0, 1, 1, 1, 0, 0, 0, 1, 1, 1, 1, 1, 1, 0, 1, 1, 1, 0, 1, 1, 1, 1, 1, 1, 1, 1, 0, 1, 1, 0, 1, 0, 0, 0, 1, 1, 0, 0, 0, 0, 1, 0, 1, 0, 1, 1, 1, 0, 0, 1, 0, 1, 1, 1, 0, 1, 0, 1, 1, 0, 0, 1, 0, 1, 1, 1, 1, 1, 1, 0, 1, 1, 0, 1, 1, 0, 0, 0, 0, 1, 1, 0, 1, 0, 0, 1, 0, 1, 0, 1, 0, 0, 0, 1, 0, 0, 1, 0, 1, 1, 0, 1, 1, 1, 0, 1, 0, 0, 1, 0, 0, 0, 1, 1, 0, 1, 0, 1, 1, 0, 1, 1, 0, 0, 1, 1, 1, 1, 1, 0, 1, 0, 0, 0, 0, 1, 0, 1, 1, 1, 1, 1, 0, 1, 0, 0, 1, 1, 0, 1, 0, 0, 1, 1, 0, 0, 1, 0, 1, 1, 0, 1, 1, 1, 1, 0, 0, 1, 1, 1, 0, 1, 1, 1, 1, 0, 1, 1, 1, 1, 1, 0, 0, 1, 0, 0, 1, 1, 1, 1, 0, 1, 0, 0, 0, 0, 1, 0, 1, 1, 1, 1, 0, 0, 0, 1, 0, 1, 1, 1, 0, 1, 0, 0, 0, 1, 1, 1, 0, 1, 0, 1, 1, 1, 0, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 1, 0, 0, 1, 1, 1, 0, 1, 1, 1, 0, 1, 1, 1, 1, 1, 1, 0, 1, 1, 1, 1, 0, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 1, 0, 0, 0, 1, 0, 1, 1, 1, 0, 0, 0, 1, 0, 1, 1, 1, 1, 0, 1, 0, 1, 0, 0, 1, 1, 1, 0, 1, 1, 0, 1, 1, 0, 1, 1, 1, 0, 1, 0, 1, 0, 0, 0, 1, 1, 0, 1, 0, 1, 0, 1, 1, 1, 0, 0, 0, 0, 0, 1, 0, 0, 0, 1, 1, 1, 1, 0, 1, 1, 0, 0, 0, 1, 1, 1, 1, 1, 0, 1, 1, 1, 1, 0, 0, 0, 0, 1, 0, 1, 1, 1, 0, 0, 0, 0, 1, 0, 0, 1, 0, 0, 0, 0, 1, 1, 0, 1, 1, 0, 1, 1, 1, 1, 0, 0, 0, 1, 1, 1, 0, 0, 1, 0, 1, 0, 1, 0, 0, 1, 1, 1, 1, 0, 1, 1, 0, 1, 0, 0, 0, 0, 1, 1, 1, 0, 1, 0, 1, 1, 1, 0, 1, 0, 0, 1, 1, 1, 0, 1, 1, 1, 1, 1, 0, 1, 1, 1, 1, 1, 1, 0, 1, 0, 1, 0, 0, 0, 0, 0, 1, 1, 1, 0, 0, 0, 0, 1, 0, 1, 1, 1, 1, 1, 0, 1, 0, 1, 1, 0, 0, 1, 0, 0, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 1, 1, 0, 1, 1, 1, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 1, 0, 1, 0, 0, 1, 0, 0, 1, 0, 1, 1, 1, 0, 0, 1, 1, 1, 1, 1, 0, 0, 0, 0, 1, 1, 1, 1, 1, 0, 0, 1, 0, 0, 0, 1, 1, 0, 1, 0, 1, 0, 0, 1, 0, 1, 0, 1, 1, 1, 0, 1, 1, 1, 0, 1, 1, 0, 1, 0, 1, 1, 1, 0, 1, 1, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 0, 1, 1, 0, 0, 0, 1, 0, 1, 1, 0, 1, 1, 1, 1, 1, 0, 1, 1, 0, 1, 0, 1, 0, 1, 1, 0, 1, 0, 0, 1, 1, 1, 0, 1, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 1, 1, 1, 0, 1, 0, 0, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 1, 0, 0, 0, 0, 1, 1, 1, 1, 1, 0, 0, 0, 1, 1, 1, 1, 1, 0, 0, 1, 0, 0, 1, 1, 1, 1, 1, 0, 0, 1, 1, 0, 1, 1, 0, 1, 0, 0, 0, 0, 1, 1, 1, 1, 0, 1, 0, 0, 1, 0, 1, 1, 0, 1, 1, 1, 0, 1, 0, 0, 1, 1, 1, 1, 1, 0, 1, 1, 0, 1, 1, 1, 0, 1, 0, 0, 1, 0, 0, 1, 0, 1, 1, 0, 1, 1, 1, 0, 1, 1, 1, 0, 1, 1, 1, 1, 0, 0, 1, 0, 0, 1, 1, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 1, 1, 0, 1, 0, 0, 1, 0, 1, 0, 0, 0, 1, 1, 1, 0, 0, 0, 1, 1, 1, 0, 1, 0, 1, 0, 0, 0, 1, 1, 0, 1, 0, 0, 1, 0, 1, 1, 1, 1, 1, 1, 1, 0, 0, 1, 1, 0, 1, 0, 0, 0, 1, 0, 0, 1, 0, 1, 0, 1, 1, 0, 1, 1, 0, 0, 1, 0, 1, 0, 1, 1, 1, 1, 1, 1, 0};
    for (int i = 0; i < 1; i++) {
        fill_random_data();
        // data_vec = {1, 1, 0, 1, 0, 1, 0, 0, 1, 0, 1, 0, 0, 1, 1, 0, 1, 0, 1, 0, 0, 0, 1, 0, 0, 0, 0, 0, 1, 1, 0, 1, 0, 1, 0, 1, 1, 1, 0, 1, 0, 1, 0, 1, 1, 0, 1, 1, 1, 0, 0, 1, 0, 0, 0, 1, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 1, 1, 1, 0, 1, 0, 0, 0, 1, 1, 0, 0, 0, 1, 0, 0, 0};
        size_n = data_vec.size();

        if (size_n % 2 == 1) {
            data_vec.insert(data_vec.end(), 1, 0);
            size_n = data_vec.size();
        }

        decoded_vec_bic = vector<bool>();
        //if (data_vec.size()%2!=0) {
        //   data_vec.insert(data_vec.end(),1,false);
        //}
        restart_param(3, 10);
        //check param
        if (2 * size_m >= size_k) {
            cout << "K and M is invalid" << endl;
            cin;
        }
        //unsigned long num_i = ((size_n)+(2*size_k-size_beta)-1)/(2*size_k-size_beta);
        encode_bbic();

        cout <<endl<<endl<<endl<<endl;
        cout << "DATA" << endl;
        cout <<endl<<endl<<endl<<endl;

        for (int i = 0; i < result_vec_bic.size(); i++) {
            printNucVectorBool(result_vec_bic[i].second);
        }

        encode_ldpc();
        cout <<endl<<endl<<endl<<endl;
        cout << "PARITY:" << endl;
        cout <<endl<<endl<<endl<<endl;

        for (int i = 0; i < size(parity); i++) {
            vector<bool> v = {0, 0, 0, 0, 0, 0, 0, 0,
                              0, 0, 0, 0, 0, 0, 0, 0,
                              0, 0, 0, 0, 0, 0, 0, 0,
                              0, 0}; // K=8

            for (int j = 0; j < 24; j++) {
                v[j] = parity[i][j];
            }
            printNucVectorBool(v);
        }
        decode_ldpc();

        // decode_bic();
        // decode_bbic_ldpc();
       decode_bbic();
        cout << "----------------" << endl;
        cout << "size_DATA: " << size_n << endl;
        cout << "size_beta: " << size_beta << endl;
        cout << "K: " << size_k << endl;
        cout << "m: " << size_m << endl;
        //cout << "need to be oligo num: " <<num_i<< endl;

        // cout << "----------ORIGINAL----------" << endl;
        // printVectorBool(data_vec);
        // printNucVectorBool(data_vec);
        // cout << "---------------------------" << endl;

        // for (int l=0; l<result_vec_bic.size(); l++) {
        //     cout << "i is: "<<result_vec_bic[l].first << endl;
        //     cout << "vec is: " << endl;
        //     printVectorBool(result_vec_bic[l].second);
        //     printNucVectorBool(result_vec_bic[l].second);
        // }

        // cout << "----------DECODE----------" << endl;
        // printVectorBool(decoded_vec_bic);
        // printNucVectorBool(decoded_vec_bic);

        cout << "----------EQ?----------" << endl;
        if (!areEqual(data_vec, decoded_vec_bic, true)) {
            return 1;
        }

        //cout << "Padded:" << pad_till_success(data_vec) << '\n';
        areEqual(data_vec, decoded_vec_bic, true);
    }

    return 0;
}
