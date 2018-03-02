#ifndef ENUMERATE_MS_INTERVALS
#define ENUMERATE_MS_INTERVALS

#include <iostream>
#include <string>
#include "String_Depth_Support.hh"
#include "Parent_Support.hh"
#include "Interfaces.hh"
#include "Basic_bitvector.hh"
#include <vector>

// Matching statistics from S to T
vector<Interval> enumerate_MS_intervals(std::string S, std::string T){
    
    vector<Interval> result;
    
    BD_BWT_index<> index((uint8_t*)T.c_str());
    std::shared_ptr<Basic_bitvector> rev_st_bpr = make_shared<Basic_bitvector>(get_rev_st_topology(index));
    rev_st_bpr->init_bps_support();
    rev_st_bpr->init_select_10_support();
    rev_st_bpr->init_rank_10_support();
    Parent_Support PS(rev_st_bpr);
    
    std::vector<bool> T_alphabet(256);
    for(char c : T){
        T_alphabet[c] = true;
    }
    
    Interval_pair I(0,index.size()-1,0,index.size()-1); // Interval pair of the empty string.
    
    for(int64_t i = 0; i < S.size(); i++){
        
        // Loop invariant: I is the Interval pair of the longest match
        // that ends at S[i-1]

        if(T_alphabet[S[i]] == false){
            // S[i] is not in the alphabet of T
            I = Interval_pair(0,index.size()-1,0,index.size()-1); // Reset to the empty string
        } else{
            Interval_pair I_right = index.right_extend(I, S[i]);
            while(I_right.reverse.size() == 0){
                I = Interval_pair(Interval(0,0), PS.lex_parent(I.reverse)); // Don't care about the forward interval
                I_right = index.right_extend(I, S[i]);
                // This loop terminates, because at least the empty string
                // is right-extensible by S[i]
            }
            I = I_right;
        }
        
        result.push_back(I.reverse);
    }
    
    return result;
}


#endif
