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

void printVectorBool(const vector<bool>& vec) {
    for (bool val : vec) {
        cout << val << " ";
    }
    cout << endl;
}

void printNucVectorBool(vector<bool> vec) {
    auto it = vec.begin();
    for (; (it != vec.end())&&(it-1!= vec.end()); it=it+2) {
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

int restart_param (int m, int k) {
    size_m = m;
    size_k = k;
    size_beta = (size_k+size_m-1)/size_m - 1;
    i_s = 0;
    i_e = 0;
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
        it->second[l-i_s] = data_vec[l];
    }

    return (addEmptyItems((it->second), size_m));

}


bool is_condition_1_satisfied(vector<bool> oli, int position) {
    bool curr_upper_bit = oli[position-(2*size_m+1)];
    bool curr_lower_bit = oli[position-2*size_m];
    int current_streak = 1;
    for (int pos = position-(2*size_m-1); (pos < position+2*size_m && (pos+1) < size_n); pos+=2) {
        if (pos+1 == position) continue;
        if ((curr_upper_bit != oli[pos])||(curr_lower_bit != oli[pos+1])) {
            current_streak = 1;
            curr_upper_bit = oli[pos];
            curr_lower_bit = oli[pos+1];
        }
        else {
            current_streak++;
            if (current_streak == size_m) return true;
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
        return 1;
    }

    for (int l = 2*size_m+1; l <= it->second.size(); l += 2*size_m) {
       if (is_condition_1_satisfied(it->second,l)&&is_condition_2_satisfied(it->second,l)) {
           it->second[l] = !(it->second[l-2]);
           }
       else {
           d_i++;
           if (i_e+d_i<data_vec.size()){
               it->second[l] = data_vec[i_e+d_i];
           }//if the last???
           else {
               cout << "ERROR: Reached unreachable code"<<endl;
               it->second[l] = false; //need to see what val at the end
               d_i--;
           }
       }
    }
    return d_i;;
}
int add_dummy_bits(int add_num) {
    bool dummy_arr[] = {false, false, false, true, true, false, true, true};
    int counter = 0;
    while (add_num > 0) {
        data_vec.insert(data_vec.end(),1,dummy_arr[counter % 8]);
        counter++;
        add_num--;
    }
    size_n = data_vec.size();
    return 0;
}
int encode_bic () {
    int oligo_num=1, temp_d;
    bool do_loop = true;
    int rest = signed(data_vec.size()-i_e);

    while (rest>=2*size_k) {
        i_e= i_s + 2*size_k-size_beta-1;
        vector<bool> oligo(2*size_k-size_beta);
        // Add the pair (i, oligo) to result_vec
        result_vec_bic.emplace_back(oligo_num, oligo);

        if (!copy_data_without_Q(oligo_num)) {
           i_s=i_e + 1 + add_bbic_bit(oligo_num);
        }
        else {
            return 1;
        }
        oligo_num++;
        rest = signed(data_vec.size()-i_s);
    }

    if (rest>0) {
        i_e = static_cast<int>(size_n) - 1; //???
        add_dummy_bits(2*size_k-rest);
        vector<bool> oligo(2*size_k-size_beta);//check???
        // Add the pair (i, oligo) to result_vec
        result_vec_bic.emplace_back(oligo_num, oligo);

        if (!copy_data_without_Q(oligo_num)) {
            i_s=i_e + 1 + add_bbic_bit(oligo_num);
        }
        else {
            return 1;
        }

    }

    return 0;
}

int take_data_without_Q (int ol_i) {
    auto it = result_vec_bic.begin();
    for (; it != result_vec_bic.end(); it++) {
        if (it->first == ol_i) {
            break;
        }
    }

    if (it == result_vec_bic.end()) {
        return 1;
    }
    for (int l = i_s;l<2*size_k;l++) {
        //TODO: Skip Q
        decoded_vec_bic [l] = it->second[l-i_s];
    }
}

int decode_bic () {
    for (auto oligo_vec:result_vec_bic) {
        take_data_without_Q(oligo_vec.first) ;


    }
}


int main() {
    //data_vec = vector<bool>(36, false);
    //data_vec.insert(data_vec.end(), 1, true);

    data_vec = {0,0, 0 ,1, 1 ,0 ,1, 1 ,0,0, 0 ,1, 1 ,0 ,1, 1,0,0, 0 ,1, 1 ,0 ,1, 1,0,0, 0 ,1, 1 ,0 ,1, 1,0,0, 0 ,1, 1 ,0 ,1, 1,0,0, 0 ,1, 1 ,0 ,1, 1,0,0, 0 ,1, 1 ,0 ,1, 1};
    size_n = data_vec.size();

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


    cout << "----------------"<< endl;
    cout << "size_DATA: " <<size_n<< endl;
    cout << "size_beta: " <<size_beta<< endl;
    cout << "K: " <<size_k<< endl;
    cout << "m: " <<size_m<< endl;
    //cout << "need to be oligo num: " <<num_i<< endl;
    cout << "----------------"<< endl;
    printVectorBool(data_vec);
    printNucVectorBool(data_vec);
    for (int l=0; l<result_vec_bic.size(); l++) {
        cout << "i is: "<<result_vec_bic[l].first << endl;
        cout << "vec is: " << endl;
        printVectorBool(result_vec_bic[l].second);
        printNucVectorBool(result_vec_bic[l].second);
    }


    return 0;
}
