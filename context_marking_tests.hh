#ifndef CONTEXT_MARKING_TESTS
#define CONTEXT_MARKING_TESTS

#include <iostream>
#include <string>
#include <vector>
#include <set>
#include <cmath>
#include <algorithm>
#include "BPR_tools.hh"
#include "brute_tools.hh"
#include "context_marking.hh"
#include "globals.hh"

using namespace std;

void test_mark_contexts_entropy(string text, double threshold){
    //cout << "Running entropy marking with threshold " << threshold << ", text " << text << endl;
    BD_BWT_index<> index((uint8_t*)text.c_str());
    sdsl::bit_vector rev_st_bpr = get_rev_st_topology(index);
    sdsl::select_support_mcl<10,2> rev_st_ss_10(&rev_st_bpr);
    sdsl::rank_support_v<10,2> rev_st_rs_10(&rev_st_bpr);
    sdsl::bp_support_g<> rev_st_bps(&rev_st_bpr);
    
    Full_Topology_Mapper mapper(&rev_st_bps, &rev_st_ss_10, &rev_st_rs_10);
    SLT_Iterator iterator(&index);
    Entropy_Formula F(threshold);
    sdsl::bit_vector contexts = F.get_rev_st_context_marks(&index, rev_st_bpr.size(), iterator, mapper);

    // Find substrings corresponding to contexts: the substring is the LCP of the first and last leaf of the node interval
    string text_rev(text.rbegin(), text.rend());
    text_rev += '$';
    vector<int64_t> SA_rev = get_suffix_array(text_rev);
    vector<string> context_strings;    

    // Find labels of all marked strings
    for(int64_t i = 0; i < contexts.size(); i++){
        if(contexts[i] == 1 && rev_st_bpr[i] == 1){
            int64_t open = i;
            int64_t close = rev_st_bps.find_close(open);
            Interval colex = bpr_interval_to_leaf_interval(Interval(open,close), rev_st_rs_10);
            context_strings.push_back(colex_range_to_string(text,colex.left,colex.right));
        }
    }
    
    sort(context_strings.begin(), context_strings.end());

    assert(context_strings == get_contexts_entropy_brute(text, threshold));
    //cerr << "Entropy marking test OK (threshold " << threshold << ")" << endl;
    
}


void test_mark_contexts_formulas_234(string text, double tau1, double tau2, double tau3, double tau4){
    //cerr << "Formulas 2,3,4 test: thresholds " << tau1 << " " << tau2 << " " << tau3 <<
    //" " << tau4 << ", text: " << text + ")" << endl;
    BD_BWT_index<> index((uint8_t*)text.c_str());
    sdsl::bit_vector rev_st_bpr = get_rev_st_topology(index);
    sdsl::select_support_mcl<10,2> rev_st_ss_10(&rev_st_bpr);
    sdsl::rank_support_v<10,2> rev_st_rs_10(&rev_st_bpr);
    sdsl::bp_support_g<> rev_st_bps(&rev_st_bpr);
    
    SLT_Iterator iterator(&index);
    Full_Topology_Mapper mapper(&rev_st_bps, &rev_st_ss_10, &rev_st_rs_10);
    EQ234_Formula F(tau1,tau2,tau3,tau4);
    sdsl::bit_vector contexts = F.get_rev_st_context_marks(&index, rev_st_bpr.size(), iterator, mapper);
    
    vector<string> context_strings;
    for(int64_t i = 0; i < contexts.size(); i++){
        if(contexts[i] == 1 && rev_st_bpr[i] == 1){
            int64_t open = i;
            int64_t close = rev_st_bps.find_close(open);
            string label = bpr_node_to_string(text,open,close,rev_st_rs_10);
            context_strings.push_back(label);
        }
    }
    
    sort(context_strings.begin(), context_strings.end());
    vector<string> brute = get_contexts_formulas234_brute(text,tau1,tau2,tau3,tau4);
    for(int64_t i = 0; i < brute.size(); i++) brute[i] = left_saturate(brute[i],text);
    sort(brute.begin(), brute.end());
    assert(context_strings == brute);
}

void test_mark_contexts_KL(string text, double threshold){
    //cerr << "KL marking test, threshold: " << threshold << ", text: " << text << endl;
    BD_BWT_index<> index((uint8_t*)text.c_str());
    sdsl::bit_vector rev_st_bpr = get_rev_st_topology(index);
    sdsl::select_support_mcl<10,2> rev_st_ss_10(&rev_st_bpr);
    sdsl::rank_support_v<10,2> rev_st_rs_10(&rev_st_bpr);
    sdsl::bp_support_g<> rev_st_bps(&rev_st_bpr);
    
    SLT_Iterator iterator(&index);
    Full_Topology_Mapper mapper(&rev_st_bps, &rev_st_ss_10, &rev_st_rs_10);
    KL_Formula F(threshold);
    sdsl::bit_vector contexts = F.get_rev_st_context_marks(&index, rev_st_bpr.size(), iterator, mapper);
    
    vector<string> context_strings;
    for(int64_t i = 0; i < contexts.size(); i++){
        if(contexts[i] == 1 && rev_st_bpr[i] == 1){
            int64_t open = i;
            int64_t close = rev_st_bps.find_close(open);
            string label = bpr_node_to_string(text,open,close,rev_st_rs_10);
            context_strings.push_back(label);
        }
    }
    
    sort(context_strings.begin(), context_strings.end());
    vector<string> brute = get_contexts_KL_brute(text,threshold);
    for(int64_t i = 0; i < brute.size(); i++) brute[i] = left_saturate(brute[i],text);
    sort(brute.begin(), brute.end());
    assert(context_strings == brute);
}

void test_mark_contexts_p_norm(string text, double p, double threshold){
    //cerr << "p_norm marking test, p: " << p << ", threshold: " << threshold << ", text: " << text << endl;
    BD_BWT_index<> index((uint8_t*)text.c_str());
    sdsl::bit_vector rev_st_bpr = get_rev_st_topology(index);
    sdsl::select_support_mcl<10,2> rev_st_ss_10(&rev_st_bpr);
    sdsl::rank_support_v<10,2> rev_st_rs_10(&rev_st_bpr);
    sdsl::bp_support_g<> rev_st_bps(&rev_st_bpr);
    
    SLT_Iterator iterator(&index);
    Full_Topology_Mapper mapper(&rev_st_bps, &rev_st_ss_10, &rev_st_rs_10);
    pnorm_Formula F(p,threshold);
    sdsl::bit_vector contexts = F.get_rev_st_context_marks(&index, rev_st_bpr.size(), iterator, mapper);
    
    vector<string> context_strings;
    for(int64_t i = 0; i < contexts.size(); i++){
        if(contexts[i] == 1 && rev_st_bpr[i] == 1){
            int64_t open = i;
            int64_t close = rev_st_bps.find_close(open);
            string label = bpr_node_to_string(text,open,close,rev_st_rs_10);
            context_strings.push_back(label);
        }
    }
    
    sort(context_strings.begin(), context_strings.end());
    vector<string> brute = get_contexts_p_norm_brute(text,p,threshold);
    for(int64_t i = 0; i < brute.size(); i++) brute[i] = left_saturate(brute[i],text);
    sort(brute.begin(), brute.end());
    assert(context_strings == brute);
}


void test_mark_contexts_entropy_all(){
    cerr << "Running entropy marking tests" << endl;
	
    string text = "abababbababbbabababababbabbbabaabbabababbabababaaabaa";
    for(double threshold = 0.1; threshold <= 3; threshold += 0.1){
        test_mark_contexts_entropy(text,threshold);
    }
    
    text = "abcabcbabcabcbabcacbbbcbcbbbabcbabcbabcbbbabbacbbcbaaababcbbbaacbcbaabcbabbabca";
    for(double threshold = 0.1; threshold <= 3; threshold += 0.1){
        test_mark_contexts_entropy(text,threshold);
    }
}

void test_mark_contexts_formulas_234_all(){
    cerr << "Running formula 2,3,4 marking tests" << endl;
    test_mark_contexts_formulas_234("aaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa",
                                    0.136848, 0.152268, 0.834861, 1.05035);
    srand(134121901);
    //string text = "abababbababbbabababababbabbbabaabbabababbabababaaababab";
    for(int64_t repeat = 0; repeat < 100; repeat++){
        int64_t alphabet_size = 2 + (rand()%2); // 2 or 3
        string text = get_random_string(100,alphabet_size);
        double tau1 = ((double)rand() / RAND_MAX) / 5;
        double tau2 = ((double)rand() / RAND_MAX);
        double tau3 = 1 - ((double)rand() / RAND_MAX) / 5;
        double tau4 = 1 + ((double)rand() / RAND_MAX) / 5;
        test_mark_contexts_formulas_234(text,tau1,tau2,tau3,tau4);
    }
}

void test_mark_contexts_KL_all(){
    cerr << "Running KL marking tests" << endl;
    test_mark_contexts_KL("aaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa",0.2);
    srand(124121901);
    //string text = "abababbababbbabababababbabbbabaabbabababbabababaaababab";
    for(int64_t repeat = 0; repeat < 100; repeat++){
        int64_t alphabet_size = 2 + (rand()%2); // 2 or 3
        string text = get_random_string(100,alphabet_size);
        double threshold = ((double)rand() / RAND_MAX);
        test_mark_contexts_KL(text,threshold);
    }
}

void test_mark_contexts_p_norm_all(){
    cerr << "Running p norm marking tests" << endl;
    test_mark_contexts_p_norm("aaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa",
                                    2,0.5);
    srand(1341211);
    //string text = "abababbababbbabababababbabbbabaabbabababbabababaaababab";
    for(int64_t repeat = 0; repeat < 100; repeat++){
        int64_t alphabet_size = 2 + (rand()%2); // 2 or 3
        string text = get_random_string(100,alphabet_size);
        int64_t p = 1 + (rand() % 5);
        double threshold = ((double)rand() / RAND_MAX);
        test_mark_contexts_p_norm(text,p,threshold);
    }
}



#endif
