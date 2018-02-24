#ifndef STRING_TOOLS_HH
#define STRING_TOOLS_HH

#include <vector>
#include <string>
#include <algorithm>
#include <utility>
#include <set>
#include "Interval.hh"
#include "sdsl/rank_support_v.hpp"
#include "BPR_tools.hh"

using namespace std;

// Computes the right-side empirical entropy of W in T
double entropy_brute(string W, string T){
    vector<int64_t> right_counts(256);
    int64_t right_extensions = 0;
    for(int64_t i = 0; i <= (int64_t)T.size() - (int64_t)W.size() - 1; i++){
        // Loop upper bound has -1, because we are not interested in an occcurrence at the
        // last position of T, because there is no information about what follows on the right
        // Casts to int64_t needed to prevent underflow if the upper bound is negative
        if(T.substr(i,W.size()) == W){
            right_counts[T[i+W.size()]]++;
            right_extensions++;
        }
    }
    double entropy = 0;
    for(int64_t i = 0; i < 256; i++){
        if(right_counts[i] != 0){
            double prob = ((double)right_counts[i] / right_extensions);
            entropy -= prob * log2(prob);
        }
    }
    return entropy;
}

// Counts the number of times W occurs in T
// If W is the empty string, returns T.size()
int64_t count_brute(string W, string T){
    int64_t ans = 0;
    for(int64_t i = 0; i < T.size(); i++){
        if(T.substr(i,W.size()) == W){
            ans++;
        }
    }
    return ans;
}

set<string> get_all_substrings(const string& T){
    set<string> all;
    for(int64_t start = 0; start < T.size(); start++){
        for(int64_t length = 0; length < T.size(); length++){
            all.insert(T.substr(start,length));
        }
    }
    return all;
}

set<char> left_extensions(const string& W, const string& T){
    set<char> ans;
    for(int64_t i = 0; i < T.size(); i++){
        if(T.substr(i,W.size()) == W){
            if(i == 0) ans.insert('$');
            else ans.insert(T[i-1]);
        }
    }
    return ans;
}

string left_saturate(string W, const string& T){
    while(true){
        set<char> left = left_extensions(W,T);
        if(left.size() == 1) W = *left.begin() + W;
        else return W;
    }
}

// W is a context iff f(W)H(W) - \sum_{a \in \Sigma} f(aW) H(aW) >= threshold
// where f is the number of occurrences of W and H(W) is the empirical entropy of W
// to the right, i.e. if P(b) = f(Wb) / f(W), then
// H(W) = - \sum_{b \in \Sigma} P(b) log P(b)
vector<string> get_contexts_entropy_brute(const string& T, double threshold, int64_t length_bound){
    set<char> alphabet;
    for(char c : T) alphabet.insert(c);
    assert(alphabet.find('$') == alphabet.end());
    set<string> contexts = {""};
    for(string W : get_all_substrings(T)){
        int64_t f_W = count_brute(W,T);
        if(T.substr(T.size()-W.size()) == W) 
            f_W--; // Don't count the last occurrence since it has no right context
        if(f_W <= 1) continue;
        double lefthand = f_W * entropy_brute(W,T);
        for(char a : alphabet){
            int64_t f_aW = count_brute(a + W, T);
            if(T.substr(T.size()-(a+W).size()) == a + W) 
                f_aW--; // Don't count the last occurrence since it has no right context 
            lefthand -= f_aW * entropy_brute(a + W, T); 
        }
        if(W.size() <= length_bound && lefthand >= threshold)
            contexts.insert(W);
    }
    
    vector<string> ans(contexts.begin(), contexts.end());
    return ans;
}

vector<string> get_contexts_entropy_brute(const string& T, double threshold){
    return get_contexts_entropy_brute(T,threshold, (int64_t)1 << 18);
}


// A context is a substring aW which *does not contain any dollars*,
// such that there exists a character b such that a, W and b satisfy equations 2,3 and 4.
// The end sentinel '$' can not be used as the character b i.e. the dollar symbol is *not* 
// counted as a right extension.
// In equation 2, the length of T is taken to *not* include the dollar at the end of T.
// The number of occurrences of the empty string is |T|.
// Root is always marked, so that the neastest marked ancestor query in scoring always succeeds.

// Based on equations 2,3 and 4 of the paper A Framework for Space-efficient String Kernels:
// 2) f(aW) / (|T| - |aW| + 1) > \tau_1
// 3) f(aWb) / f(aW) > \tau_2
// 4) (f(aWb) / f(aW)) / (f(Wb) / f(W)) \in (0..\tau_3] \cup [\tau_4..\infty]

vector<string> get_contexts_formulas234_brute(string T, double tau1, double tau2, double tau3, double tau4, int64_t length_bound){
    set<char> alphabet;
    for(char c : T) alphabet.insert(c);
    assert(alphabet.find('$') == alphabet.end());
    set<string> contexts;
    contexts.insert(""); // Root is always marked
    for(int64_t aW_start = 0; aW_start < T.size(); aW_start++){
        for(int64_t aW_end = aW_start; aW_end < T.size(); aW_end++){
            int64_t aW_length = aW_end - aW_start + 1;
            char a = T[aW_start];
            string W = T.substr(aW_start+1,aW_length-1);
            double f_aW = count_brute(a+W,T);
            double eq2 = f_aW / ((T.size()) - (1 + W.size()) + 1);
            if(eq2 < tau1) continue;
            for(char b : alphabet){
                double f_aWb = count_brute(a+W+b,T);
                double f_Wb = count_brute(W+b,T);
                double f_W = count_brute(W,T);
                double eq3 = f_aWb / f_aW ;
                double eq4_numerator = f_aWb / f_aW;
                double eq4_denominator = f_Wb / f_W;
                if(aW_length <= length_bound && eq3 >= tau2 && eq4_denominator != 0 && 
                  (eq4_numerator / eq4_denominator <= tau3 || eq4_numerator / eq4_denominator >= tau4)){
                    contexts.insert(a+W);
                    break;
                }
            }
        }
    }
    vector<string> ans(contexts.begin(), contexts.end());
    return ans;
}

vector<string> get_contexts_formulas234_brute(string T, double tau1, double tau2, double tau3, double tau4){
    return get_contexts_formulas234_brute(T,tau1,tau2,tau3,tau4,(int64_t)1 << 18);
}

vector<string> get_contexts_KL_brute(string T, double threshold, int64_t length_bound){
    set<char> alphabet;
    for(char c : T) alphabet.insert(c);
    assert(alphabet.find('$') == alphabet.end());
    set<string> contexts;
    contexts.insert(""); // Root is always marked
    for(int64_t aW_start = 0; aW_start < T.size(); aW_start++){
        for(int64_t aW_end = aW_start; aW_end < T.size(); aW_end++){
            int64_t aW_length = aW_end - aW_start + 1;
            char a = T[aW_start];
            string W = T.substr(aW_start+1,aW_length-1);
            double f_aW = count_brute(a+W,T);
            double f_W = count_brute(W,T);
            double KL_divergence = 0;
            for(char b : alphabet){
                double f_aWb = count_brute(a+W+b,T);
                double f_Wb = count_brute(W+b,T);
                if(f_aW != 0 && f_W != 0 && f_aWb != 0 && f_Wb != 0)
                    KL_divergence += f_aWb * log((f_aWb / f_aW) / (f_Wb / f_W));
            }
            if(aW_length <= length_bound && KL_divergence >= threshold)
                contexts.insert(a + W);
        }
    }
    vector<string> ans(contexts.begin(), contexts.end());
    return ans;
}

vector<string> get_contexts_KL_brute(string T, double threshold){
    return get_contexts_KL_brute(T, threshold, (int64_t) 1 << 18);
}

vector<string> get_contexts_p_norm_brute(string T, double p, double threshold, int64_t length_bound){
    set<char> alphabet;
    for(char c : T) alphabet.insert(c);
    assert(alphabet.find('$') == alphabet.end());
    set<string> contexts;
    contexts.insert(""); // Root is always marked
    for(int64_t aW_start = 0; aW_start < T.size(); aW_start++){
        for(int64_t aW_end = aW_start; aW_end < T.size(); aW_end++){
            int64_t aW_length = aW_end - aW_start + 1;
            char a = T[aW_start];
            string W = T.substr(aW_start+1,aW_length-1);
            double f_aW = count_brute(a+W,T);
            if(f_aW == 0) continue;
            double f_W = count_brute(W,T);
            double p_norm = 0;
            for(char b : alphabet){
                double f_aWb = count_brute(a+W+b,T);
                double f_Wb = count_brute(W+b,T);
                p_norm += pow(fabs(f_aWb / f_aW - f_Wb / f_W),p);
            }
            p_norm = f_aW * pow(p_norm, 1.0/p);
            if(aW_length <= length_bound && p_norm >= threshold)
                contexts.insert(a + W);
        }
    }
    vector<string> ans(contexts.begin(), contexts.end());
    return ans;
}

vector<string> get_contexts_p_norm_brute(string T, double p, double threshold){
    return get_contexts_p_norm_brute(T,p,threshold,(int64_t) 1 << 18);
}


vector<int64_t> get_suffix_array(string& s){
    int64_t n = s.size();
    vector<int64_t> v(n);
    for(int64_t i = 0; i < n; i++) v[i] = (int64_t) s[i];
    for(int64_t L = 1; L <= n; L *= 2){
        vector<pair<pair<int64_t,int64_t>, int64_t> > pairs(n);
        for(int64_t i = 0; i < n; i++){
            pairs[i] = {{v[i], (i + L < n) ? v[i+L] : -1}, i};
        }
        sort(pairs.begin(), pairs.end());
        int64_t rank = 0;
        for(int64_t i = 0; i < n; i++){
            if(i > 0 && pairs[i].first != pairs[i-1].first) rank++;
            v[pairs[i].second] = rank;
        }
   }
   vector<int64_t> SA(n);
   for(int64_t i = 0; i < n; i++) SA[v[i]] = i;
   return SA;
}

// slow, for debug purposes
string colex_range_to_string(string S, int64_t left, int64_t right){
    //assert(left <= right);
    string S_rev(S.rbegin(), S.rend());
    S_rev += "$";
    
    vector<int64_t> SA_rev = get_suffix_array(S_rev);
    
    if(left == right){
        string label = S_rev.substr(SA_rev[left]);
        reverse(label.begin(), label.end());
        return label;
    }
    
    int64_t i1 = SA_rev[left];
    int64_t i2 = SA_rev[right];
    int64_t lcp = 0;
    while(S_rev[i1] == S_rev[i2]){
        lcp++;
        i1++; i2++;
    }
    string label = S_rev.substr(SA_rev[left],lcp);  // SA_rev[lex.right] would give the same
    reverse(label.begin(), label.end());
    return label;

}

//slow, for debug purposes
string bpr_node_to_string(string S, int64_t open, int64_t close, sdsl::rank_support_v<10,2>& rev_st_rs_10){
    Interval colex = bpr_interval_to_leaf_interval(Interval(open,close), rev_st_rs_10);
    return colex_range_to_string(S,colex.left,colex.right);

}

vector<string> all_binary_strings_up_to(int64_t k){
    vector<string> ans;
    for(int64_t length = 1; length <= k; length++){
        for(int64_t mask = 0; mask < (1 << length); mask++){
            string s = "";
            for(int64_t i = 0; i < length; i++){
                if(mask & (1 << i)) s += 'a';
                else s += 'b';
            }
            ans.push_back(s);
        }
    }
    return ans;
}

string get_random_string(int64_t length, int64_t alphabet_size){
    string s;
    for(int64_t i = 0; i < length; i++){
        s.push_back('a' + rand() % alphabet_size);
    }
    return s;
}

#endif
