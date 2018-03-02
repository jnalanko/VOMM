#ifndef MS_ENUMERATOR_TESTS_HH
#define MS_ENUMERATOR_TESTS_HH

#include <iostream>
#include <vector>
#include <string>
#include <utility>
#include "enumerate_MS_intervals.hh"
#include "brute_tools.hh"
#include <algorithm>

using namespace std;

vector<Interval> enumerate_MS_intervals_naive(string S, string T, vector<string>& matches){
    string T_rev(T.rbegin(), T.rend());
    T_rev += "$";
    vector<int64_t> T_rev_SA = get_suffix_array(T_rev);
    
    vector<Interval> answers;
    
    for(int64_t i = 0; i < (int64_t)S.size(); i++){
        int64_t length = 0;
        for(int64_t j = 0; j < (int64_t)T.size(); j++){
            for(int64_t k = 0; k <= min(i,j); k++){
                if(S[i-k] == T[j-k]) length = max(length, k+1);
                else break;
            }
        }
        string match = S.substr(i-length+1,length);
        
        // Find the interval of the match in the reverse SA of T
        string match_rev(match.rbegin(), match.rend());
        matches.push_back(match_rev);

        Interval I(1e9,-1e9);
        for(int64_t lex = 0; lex < (int64_t)T.size(); lex++){
            string candidate = T_rev.substr(T_rev_SA[lex],length);
            if(candidate == match_rev){
                I.left = min(I.left, (int64_t)lex);
                I.right = max(I.right, (int64_t)lex);
            }
        }
        
        answers.push_back(I);
    }
    
    return answers;
}

void MS_Enumerator_tests(){
    string T = "ababbbabaaabababababababbaabababababababababababaaabababaa";
    string S = "abbababaaaababbababababababaabaababaa";
    vector<string> matches;
    vector<Interval> intervals = enumerate_MS_intervals(S,T);
    vector<Interval> naive_intervals = enumerate_MS_intervals_naive(S,T,matches);
    assert(intervals.size() == naive_intervals.size());
    for(int64_t i = 0; i < intervals.size(); i++){
        //cout << i << " " << matches[i] << " " << intervals[i].toString() << " " << naive_intervals[i].toString() << endl;
        assert(intervals[i] == naive_intervals[i]);
    }
    cout << "MS enumerator test OK" << endl;
}

#endif
