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
vector<pair<int,vector<bool>>> result_vec;
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
    result.resize(2*size_k);
    return 0;
}

int copy_data_without_Q( int ol_i) {
    auto it = result_vec.begin();
    for (; it != result_vec.end(); it++) {
        if (it->first == ol_i) {
            break;
        }
    }

    if (it == result_vec.end()) {
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
    for (int pos = position-(2*size_m-1); (pos <= position+2*size_m && pos < size_n); pos+=2) {
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
    auto it = result_vec.begin();
    for (; it != result_vec.end(); it++) {
        if (it->first == ol_i) {
            break;
        }
    }

    if (it == result_vec.end()) {
        return 1;
    }

    for (int l = 2*size_m+1; l <= it->second.size(); l += 2*size_m) {
       if (is_condition_1_satisfied(it->second,l)&&is_condition_2_satisfied(it->second,l)) {
           it->second[l] = !(it->second[l-2]);
           }
       else {
           d_i++;
           it->second[l] = data_vec[i_e+d_i]; //if the last???
       }
    }
    return d_i;;
}
int add_dummy_bits(int add_num) {
    int counter = 0;
    while (add_num > 0) {
        switch (counter%8) {
            case 0:
            case 1:
            case 2:
            case 5:
                data_vec.insert(data_vec.end(),1,0);
            break;
            case 3:
            case 4:
            case 6:
            case 7:
                data_vec.insert(data_vec.end(),1,1);
            break;
            default:
                break;
        }
        counter++;
        add_num--;
    }
    return 0;
}
int encode_bic () {
    int i=1, temp_d;
    bool do_loop = true;


    while (signed(data_vec.size()-i_e)>=2*size_k-size_beta) {
        i_e= i_s + 2*size_k-size_beta-1;
        vector<bool> oligo(2*size_k);
        // Add the pair (i, oligo) to result_vec
        result_vec.emplace_back(i, oligo);

        if (!copy_data_without_Q(i)) {
           i_s=i_e + 1 + add_bbic_bit(i);
        }
        else {
            return 1;
        }
        i++;
    }

    if (signed(data_vec.size()-i_e)>0) {
        add_dummy_bits(2*size_k-size_beta-(data_vec.size()-i_e));
    }

    while (signed(data_vec.size()-i_e)>=2*size_k-size_beta) {
        i_e= i_s + 2*size_k-size_beta-1;
        vector<bool> oligo(2*size_k);
        // Add the pair (i, oligo) to result_vec
        result_vec.emplace_back(i, oligo);

        if (!copy_data_without_Q(i)) {
            i_s=i_e + 1 + add_bbic_bit(i);
        }
        else {
            return 1;
        }
        i++;
    }



    return 0;
}




int main() {
    data_vec = vector<bool>(36, false);
    data_vec.insert(data_vec.end(), 1, true);

    //data_vec = {0,1, 0 ,1, 1 ,1 ,0, 1 ,0 ,0, 1,0 ,1 ,0 ,0 ,0, 0, 1, 0 ,1 ,0};
    size_n = data_vec.size();

    //if (data_vec.size()%2!=0) {
     //   data_vec.insert(data_vec.end(),1,false);
    //}
    restart_param(4,10);
    //check param
    if (2*size_m >= size_k) {
        cout << "K and M is invalid"<< endl;
    }
    unsigned long num_i = ((size_n)+(2*size_k-size_beta)-1)/(2*size_k-size_beta);
    encode_bic();


    cout << "----------------"<< endl;
    cout << "size_DATA: " <<size_n<< endl;
    cout << "size_beta: " <<size_beta<< endl;
    cout << "K: " <<size_k<< endl;
    cout << "m: " <<size_m<< endl;
    cout << "need to be oligo num: " <<num_i<< endl;
    cout << "----------------"<< endl;
    printVectorBool(data_vec);
    printNucVectorBool(data_vec);
    for (int l=0; l<result_vec.size(); l++) {
        cout << "i is: "<<result_vec[l].first << endl;
        cout << "vec is: " << endl;
        printVectorBool(result_vec[l].second);
        printNucVectorBool(result_vec[l].second);
    }


    return 0;
}
