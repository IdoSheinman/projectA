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
using namespace std;

#define NUC_A 0
#define NUC_T 1
#define NUC_C 2
#define NUC_G 3

//GLOBAL VAR
vector<bool> data_vec;
vector<pair<int,vector<bool>>> result_vec_bic;
vector<bool> decoded_vec_bic;
int i_s,i_e, size_m, size_k, size_beta;
unsigned long size_n;
bool used_all_data;
int counter;
bool add_dummy_bits();

void printVectorBool(const vector<bool>& vec) {
    for (bool val : vec) {
        cout << val << " ";
    }
    cout << endl;
}

void printNucVectorBool(vector<bool> vec) {
    auto it = vec.begin();
    for (; (it != vec.end())&&(it+1!= vec.end()); it=it+2) {
        if ((!*it)&&(!*(it+1))) {
            cout << "A" << " ";
        }
        if ((!*it)&&(*(it+1))) {
            cout << "T" << " ";
        }
        if ((*it)&&(!*(it+1))) {
            cout << "C" << " ";
        }
        if (*it&&*(it+1)) {
            cout << "G" << " ";
        }
    }
    cout << endl;
}

bool areEqual(const vector<bool>& vec1, const vector<bool>& vec2, bool verbose=false) {
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

int restart_param (int m, int k) {
    size_m = m;
    size_k = k;
    size_beta = (size_k+size_m-1)/size_m - 1;
    i_s = 0;
    i_e = 0;
    used_all_data = false;
    counter = 0;
    return 0;
}

int addEmptyItems(std::vector<bool>& result, int m) {
    for (int i = 2*m+1; i <= result.size(); i += 2*m) {
        result.insert(result.begin() + i, false);
    }
    //result.resize(2*size_k);
    return 0;
}

int copy_data_without_Q( int ol_i) {
    auto it = result_vec_bic.begin();
    for (; it != result_vec_bic.end(); it++) {
        if (it->first == ol_i) {
            break;
        }
    }

    if (it == result_vec_bic.end()) {
        return 1;
    }


    for (int l=i_s; l<=i_e; l++) {
        if (l < data_vec.size()) {
            it->second[l-i_s] = data_vec[l];
        } else {
            it->second[l-i_s] = add_dummy_bits();
            counter++;
        }
    }

    return (addEmptyItems((it->second), size_m));

}


bool is_condition_1_satisfied(vector<bool> oli, int position) {
    int first_upper_pos = position-(2*size_m+1);
    int first_lower_pos = position-(2*size_m);

    bool curr_upper_bit;
    bool curr_lower_bit;
    int current_streak = 1;
    for (int i = 2; i < 2*(2*size_m); i += 2) {
        if (first_lower_pos+i>=oli.size()) {
            return false;
        }
        if (i == 2*size_m) {
            continue;
        }
        curr_upper_bit = oli[first_upper_pos+i];
        curr_lower_bit = oli[first_lower_pos+i];
        if (i == 2*size_m+2) {
            if (curr_upper_bit == oli[first_upper_pos+i-4] &&
            curr_lower_bit == oli[first_lower_pos+i-4]) {
                current_streak += 1;
            } else {
                current_streak = 1;
            }
        } else {
            if (curr_upper_bit == oli[first_upper_pos+i-2] &&
            curr_lower_bit == oli[first_lower_pos+i-2]) {
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
    return (oli[position-1]==oli[position-3]);
}


int add_bbic_bit (int ol_i) {
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

    for (int l = 2*size_m+1; l < it->second.size(); l += 2*size_m) {
        if (is_condition_1_satisfied(it->second,l)&&is_condition_2_satisfied(it->second,l)) {
           // TODO: Shouldnt this be l-1?
           it->second[l] = !(it->second[l-2]);
       }
       else {
           d_i++;
           //cout << "ADDING Q BIT IN: " << l << " from: " << i_e+d_i << endl;
           if (i_e+d_i<data_vec.size()){
               it->second[l] = data_vec[i_e+d_i];
               if (i_e+d_i==data_vec.size()-1) {
                   used_all_data=true;
               }
           }//if the last???
           else {
               cout << "ERROR: Reached unreachable code"<<endl;
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
int encode_bic () {
    int oligo_num=1;
    int rest;

    while (!used_all_data) {
        i_e= i_s + 2*size_k-size_beta-1;
        vector<bool> oligo(2*size_k-size_beta);
        // Add the pair (i, oligo) to result_vec
        result_vec_bic.emplace_back(oligo_num, oligo);

        if (!copy_data_without_Q(oligo_num)) {
            int added = add_bbic_bit(oligo_num);
           i_s=i_e + 1 + added;
        }
        else {
            return 1;
        }
        oligo_num++;
        rest = signed(data_vec.size()-i_s);
        if (rest<=0) {
            used_all_data = true;
        }
    }

    return 0;
}

int add_Q_data (int ol_i) {
    unsigned int first_q = 2*size_m+1;
    auto it = result_vec_bic.begin();
    for (; it != result_vec_bic.end(); it++) {
        if (it->first == ol_i) {
            break;
        }
    }

    if (it == result_vec_bic.end()) {
        return 1;
    }
    for (int l = first_q; l <= it->second.size(); l += 2*size_m) {
        if (!(is_condition_1_satisfied(it->second,l)&&is_condition_2_satisfied(it->second,l))) {
            decoded_vec_bic.insert(decoded_vec_bic.end(),1, it->second[l]);
        }
    }
    return 0;
}


int take_data_without_Q (int ol_i) {
    auto it = result_vec_bic.begin();
    for (; it != result_vec_bic.end(); ++it) {
        if (it->first == ol_i) {
            break;
        }
    }

    if (it == result_vec_bic.end()) {
        return 1;
    }
    unsigned int first_q = 2*size_m+1;
    for (int l = 0;l<2*size_k;l++) {
        //TODO: Skip Q
        if (l<first_q || ((l-first_q)%(2*size_m) != 0)) {
            decoded_vec_bic.insert(decoded_vec_bic.end(),1, it->second[l]);
        }
    }

    return 0;
}

int decode_bic () {
    for (auto oligo_vec:result_vec_bic) {
        take_data_without_Q(oligo_vec.first) ;
        add_Q_data (oligo_vec.first);
    }
    return 0;
}

void fill_random_data() {
    // Constants for 0 to 1 MB size in bits
    const size_t MAX_BITS = 128*8*8; // 1 KB in bits

    // Seed random number generator
    random_device rd;    // Random device for seed
    mt19937 gen(rd());   // Mersenne Twister RNG
    uniform_int_distribution<size_t> size_dist(0, MAX_BITS); // Size between 0 and 1MB
    uniform_int_distribution<int> bit_dist(0, 1);           // Values 0 or 1

    // Generate a random size
    size_t random_size = size_dist(gen);

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
            data_vec.insert(data_vec.end(),1,0);
        }

        decoded_vec_bic = vector<bool>();
        restart_param(3,10);
        if (2*size_m >= size_k) {
            cout << "K and M is invalid"<< endl;
            cin;
        }
        encode_bic();
        decode_bic();
        counter+=1;
    }
    return counter;
}

int main() {
    //data_vec = vector<bool>(36, false);
   // data_vec.insert(data_vec.begin(), 1, true); 0,0, 0 ,1, 1 ,0 ,1, 1 ,0,0, 0 ,1, 1 ,0 ,1, 1,0,0, 0 ,1, 1 ,0 ,1, 1,0,0, 0 ,1, 1 ,0 ,1, 1,0,0, 0 ,1, 1 ,0 ,

    //0 0 0 1 0 1 0 1 0 1 0 0 0 1 0 1 0 0 1 1 1 0 0 0 0 0 1 0 0 1 0 0 0 0 0 0 1 1 0
    //data_vec = {0, 0, 0, 1, 0, 1, 0, 1, 0, 1, 0, 0, 0, 1, 0, 1, 0, 0, 1, 1, 1, 0, 0, 0, 0};
    //data_vec = {0, 0, 1, 1, 0, 1, 1, 1, 0, 1, 1, 0, 0, 1, 1, 0, 1, 0, 0, 0, 1, 1, 1, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0, 1, 0, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 1, 1, 1, 1, 0, 0, 1, 0, 0, 1, 1, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0, 1, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 1, 0, 1, 0, 0, 1, 0, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 1, 1, 1, 0, 0, 1, 1, 1, 0, 0, 0, 1, 1, 1, 1, 0, 0, 1, 0, 0, 1, 1, 1, 0, 1, 0, 0, 0, 1, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 1, 0, 1, 0, 0, 1, 1, 1, 0, 0, 0, 1, 1, 1, 1, 1, 1, 0, 1, 1, 1, 0, 1, 1, 1, 1, 1, 1, 1, 1, 0, 1, 1, 0, 1, 0, 0, 0, 1, 1, 0, 0, 0, 0, 1, 0, 1, 0, 1, 1, 1, 0, 0, 1, 0, 1, 1, 1, 0, 1, 0, 1, 1, 0, 0, 1, 0, 1, 1, 1, 1, 1, 1, 0, 1, 1, 0, 1, 1, 0, 0, 0, 0, 1, 1, 0, 1, 0, 0, 1, 0, 1, 0, 1, 0, 0, 0, 1, 0, 0, 1, 0, 1, 1, 0, 1, 1, 1, 0, 1, 0, 0, 1, 0, 0, 0, 1, 1, 0, 1, 0, 1, 1, 0, 1, 1, 0, 0, 1, 1, 1, 1, 1, 0, 1, 0, 0, 0, 0, 1, 0, 1, 1, 1, 1, 1, 0, 1, 0, 0, 1, 1, 0, 1, 0, 0, 1, 1, 0, 0, 1, 0, 1, 1, 0, 1, 1, 1, 1, 0, 0, 1, 1, 1, 0, 1, 1, 1, 1, 0, 1, 1, 1, 1, 1, 0, 0, 1, 0, 0, 1, 1, 1, 1, 0, 1, 0, 0, 0, 0, 1, 0, 1, 1, 1, 1, 0, 0, 0, 1, 0, 1, 1, 1, 0, 1, 0, 0, 0, 1, 1, 1, 0, 1, 0, 1, 1, 1, 0, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 1, 0, 0, 1, 1, 1, 0, 1, 1, 1, 0, 1, 1, 1, 1, 1, 1, 0, 1, 1, 1, 1, 0, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 1, 0, 0, 0, 1, 0, 1, 1, 1, 0, 0, 0, 1, 0, 1, 1, 1, 1, 0, 1, 0, 1, 0, 0, 1, 1, 1, 0, 1, 1, 0, 1, 1, 0, 1, 1, 1, 0, 1, 0, 1, 0, 0, 0, 1, 1, 0, 1, 0, 1, 0, 1, 1, 1, 0, 0, 0, 0, 0, 1, 0, 0, 0, 1, 1, 1, 1, 0, 1, 1, 0, 0, 0, 1, 1, 1, 1, 1, 0, 1, 1, 1, 1, 0, 0, 0, 0, 1, 0, 1, 1, 1, 0, 0, 0, 0, 1, 0, 0, 1, 0, 0, 0, 0, 1, 1, 0, 1, 1, 0, 1, 1, 1, 1, 0, 0, 0, 1, 1, 1, 0, 0, 1, 0, 1, 0, 1, 0, 0, 1, 1, 1, 1, 0, 1, 1, 0, 1, 0, 0, 0, 0, 1, 1, 1, 0, 1, 0, 1, 1, 1, 0, 1, 0, 0, 1, 1, 1, 0, 1, 1, 1, 1, 1, 0, 1, 1, 1, 1, 1, 1, 0, 1, 0, 1, 0, 0, 0, 0, 0, 1, 1, 1, 0, 0, 0, 0, 1, 0, 1, 1, 1, 1, 1, 0, 1, 0, 1, 1, 0, 0, 1, 0, 0, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 1, 1, 0, 1, 1, 1, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 1, 0, 1, 0, 0, 1, 0, 0, 1, 0, 1, 1, 1, 0, 0, 1, 1, 1, 1, 1, 0, 0, 0, 0, 1, 1, 1, 1, 1, 0, 0, 1, 0, 0, 0, 1, 1, 0, 1, 0, 1, 0, 0, 1, 0, 1, 0, 1, 1, 1, 0, 1, 1, 1, 0, 1, 1, 0, 1, 0, 1, 1, 1, 0, 1, 1, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 0, 1, 1, 0, 0, 0, 1, 0, 1, 1, 0, 1, 1, 1, 1, 1, 0, 1, 1, 0, 1, 0, 1, 0, 1, 1, 0, 1, 0, 0, 1, 1, 1, 0, 1, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 1, 1, 1, 0, 1, 0, 0, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 1, 0, 0, 0, 0, 1, 1, 1, 1, 1, 0, 0, 0, 1, 1, 1, 1, 1, 0, 0, 1, 0, 0, 1, 1, 1, 1, 1, 0, 0, 1, 1, 0, 1, 1, 0, 1, 0, 0, 0, 0, 1, 1, 1, 1, 0, 1, 0, 0, 1, 0, 1, 1, 0, 1, 1, 1, 0, 1, 0, 0, 1, 1, 1, 1, 1, 0, 1, 1, 0, 1, 1, 1, 0, 1, 0, 0, 1, 0, 0, 1, 0, 1, 1, 0, 1, 1, 1, 0, 1, 1, 1, 0, 1, 1, 1, 1, 0, 0, 1, 0, 0, 1, 1, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 1, 1, 0, 1, 0, 0, 1, 0, 1, 0, 0, 0, 1, 1, 1, 0, 0, 0, 1, 1, 1, 0, 1, 0, 1, 0, 0, 0, 1, 1, 0, 1, 0, 0, 1, 0, 1, 1, 1, 1, 1, 1, 1, 0, 0, 1, 1, 0, 1, 0, 0, 0, 1, 0, 0, 1, 0, 1, 0, 1, 1, 0, 1, 1, 0, 0, 1, 0, 1, 0, 1, 1, 1, 1, 1, 1, 0};
    //for (int i=0; i<1; i++) {
        fill_random_data();
        //data_vec = {1, 1, 0, 1, 0, 1, 0, 0, 1, 0, 1, 0, 0, 1, 1, 0, 1, 0, 1, 0, 0, 0, 1, 0, 0, 0, 0, 0, 1, 1, 0, 1, 0, 1, 0, 1, 1, 1, 0, 1, 0, 1, 0, 1, 1, 0, 1, 1, 1, 0, 0, 1, 0, 0, 0, 1, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 1, 1, 1, 0, 1, 0, 0, 0, 1, 1, 0, 0, 0, 1, 0, 0, 0};
        size_n = data_vec.size();

        if (size_n%2 == 1) {
            data_vec.insert(data_vec.end(), 1, 0);
            size_n = data_vec.size();
        }

        decoded_vec_bic = vector<bool>();
        //if (data_vec.size()%2!=0) {
        //   data_vec.insert(data_vec.end(),1,false);
        //}
        restart_param(3,10);
        //check param
        if (2*size_m >= size_k) {
            cout << "K and M is invalid"<< endl;
            cin;
        }
        //unsigned long num_i = ((size_n)+(2*size_k-size_beta)-1)/(2*size_k-size_beta);
        encode_bic();
        decode_bic();

        cout << "----------------"<< endl;
        cout << "size_DATA: " <<size_n<< endl;
        cout << "size_beta: " <<size_beta<< endl;
        cout << "K: " <<size_k<< endl;
        cout << "m: " <<size_m<< endl;
        //cout << "need to be oligo num: " <<num_i<< endl;
        cout << "----------ORIGINAL----------"<< endl;
        printVectorBool(data_vec);
        printNucVectorBool(data_vec);
        cout << "---------------------------"<< endl;
        for (int l=0; l<result_vec_bic.size(); l++) {
            cout << "i is: "<<result_vec_bic[l].first << endl;
            cout << "vec is: " << endl;
            printVectorBool(result_vec_bic[l].second);
            printNucVectorBool(result_vec_bic[l].second);
        }

        cout << "----------DECODE----------"<< endl;
        printVectorBool(decoded_vec_bic);
        printNucVectorBool(decoded_vec_bic);

        cout << "----------EQ?----------"<< endl;
        if (!areEqual(data_vec, decoded_vec_bic, true)) {
         return 1;
        }

        //cout << "Padded:" << pad_till_success(data_vec) << '\n';
        areEqual(data_vec, decoded_vec_bic, true);
    //}
    return 0;
}
