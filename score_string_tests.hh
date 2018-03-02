#ifndef SCORE_STRING_TESTS
#define SCORE_STRING_TESTS

#include "score_string.hh"
#include "BWT_iteration.hh"
#include <vector>
#include <string>
#include <set>

using namespace std;

// Empirical probability that W is followed by c in S
double get_emission_prob(string T, string W, char c){
    double total = 0;
    double matching = 0;
    for(int64_t i = 0; i < T.size(); i++){
        if(T.substr(i,W.size()) == W){
            if(i + W.size() < T.size() && T[i+W.size()] == c) 
                matching++;
            total++;
        }
    }
    //cout << "get " << W << " " << c << " " << matching << " " << total << endl;
    return matching / total;
}

double score_string_brute_given_contexts(string S, string T, set<string> contexts, double escape_prob){
    double score = 0;
    //cout << S << endl;
    //cout << T << endl;
    for(int64_t i = 0; i < S.size(); i++){
        // Try candidate contexts in decreasing order of length
        for(int64_t context_start = 0; context_start <= i; context_start++){
            string W = S.substr(context_start, i-context_start);
            if(contexts.find(W) != contexts.end()){
                double prob = get_emission_prob(T,W,S[i]);
                if(prob != 0) {
                    score += log2(get_emission_prob(T,W,S[i]));
                }
                else{
                    score += log2(escape_prob);
                }
                break;
            }
            assert(i != context_start); // Empty string must always be a context
        }
    }
    return score;   
}

double score_string_brute_depth_bounded_given_contexts(string S, string T, set<string> contexts, double escape_prob, int64_t depth_bound){
    double score = 0;
    for(int64_t i = 0; i < S.size(); i++){
        // Try candidate contexts in decreasing order of length
        for(int64_t context_start = 0; context_start <= i; context_start++){
            string W = S.substr(context_start, i-context_start);
            if(W.size() > depth_bound) continue;
            if(contexts.find(W) != contexts.end()){
                double prob = get_emission_prob(T,W,S[i]);
                if(prob != 0) {
                    score += log2(get_emission_prob(T,W,S[i]));
                }
                else{
                    score += log2(escape_prob);
                }
                break;
            }
            assert(i != context_start); // Empty string must always be a context
        }
    }
    return score;   
}

double score_string_entropy_brute(string S, string T, double threshold, double escape_prob){
    vector<string> contexts_vec = get_contexts_entropy_brute(T, threshold);
    set<string> contexts(contexts_vec.begin(), contexts_vec.end());
    return score_string_brute_given_contexts(S,T,contexts,escape_prob);
}

double score_string_eq_234_brute(string S, string T, double t1, double t2, double t3, double t4, double escape_prob){
    vector<string> context_vec = get_contexts_formulas234_brute(T,t1,t2,t3,t4);
    set<string> contexts(context_vec.begin(), context_vec.end());
    return score_string_brute_given_contexts(S,T,contexts,escape_prob);
}

double score_string_KL_brute(string S, string T, double threshold, double escape_prob){
    vector<string> context_vec = get_contexts_KL_brute(T,threshold);
    set<string> contexts(context_vec.begin(), context_vec.end());
    return score_string_brute_given_contexts(S,T,contexts,escape_prob);
}

double score_string_p_norm_brute(string S, string T,  double p, double threshold, double escape_prob){
    vector<string> context_vec = get_contexts_p_norm_brute(T,p,threshold);
    set<string> contexts(context_vec.begin(), context_vec.end());
    return score_string_brute_given_contexts(S,T,contexts,escape_prob);
}

double score_string_entropy_non_brute(string S, string T, bool maxreps_only, double threshold, double escape_prob){
    if(maxreps_only){
        typedef SLT_Iterator slt_iterator_t;
        typedef Rev_ST_Maxrep_Iterator rev_slt_iterator_t;
        typedef Basic_Scorer scorer_t;
        typedef Maxrep_Pruned_Updater updater_t ;
        
        slt_iterator_t slt_it;
        rev_slt_iterator_t rev_slt_it;
        
        Entropy_Formula formula(threshold);
        Global_Data G; 
        build_model(G, T, formula, slt_it, rev_slt_it, slt_it, true, false);
        scorer_t scorer(escape_prob, true);
        updater_t updater;
        return score_string(S, G, scorer, updater);
    }
    
    else{
        typedef SLT_Iterator slt_iterator_t;
        typedef Rev_ST_Iterator rev_slt_iterator_t;
        typedef Basic_Scorer scorer_t;
        typedef Basic_Updater updater_t;
        
        slt_iterator_t slt_it;
        rev_slt_iterator_t rev_slt_it;
        Entropy_Formula formula(threshold);
        Global_Data G;
        build_model(G, T, formula, slt_it, rev_slt_it, slt_it, true, false);
        scorer_t scorer(escape_prob, true);
        updater_t updater;
        return score_string(S, G, scorer, updater);    
    }
    
}

void test_entropy(string S, string T, double threshold, double escape){
    //typedef Rev_ST_Iterator<BD_BWT_index<> > FULL;
    //typedef Rev_ST_Maxrep_Iterator<BD_BWT_index<> > MR;
    
    //Rev_ST_Iterator<BD_BWT_index<> > rev_slt_it_full;
    //Rev_ST_Maxrep_Iterator<BD_BWT_index<> > rev_slt_it_mr;

    double brute = score_string_entropy_brute(S,T,threshold,escape);
    double non_brute_1 = score_string_entropy_non_brute(S,T,false,threshold,escape);
    double non_brute_2 = score_string_entropy_non_brute(S,T,true,threshold,escape);
    assert(fabs(brute-non_brute_1) < 1e-6);
    assert(fabs(brute-non_brute_2) < 1e-6);
}

void test_serialization(){
    cerr << "Running serialization tests" << endl;
    
    srand(11112223);
    string T = get_random_string(300,3);
    string S = get_random_string(300,3);
    double threshold = 0.2;
    double escape_prob = 0.05;
         
    SLT_Iterator slt_it;
    Rev_ST_Maxrep_Iterator rev_slt_it;
    Entropy_Formula formula(threshold);
    Global_Data G1;
    build_model(G1, T, formula, slt_it, rev_slt_it, slt_it, true, false);
    G1.store_all_to_disk("models","test");
    Global_Data G2;
    G2.load_all_from_disk("models","test",true);
    
    Basic_Scorer scorer(escape_prob, true);
    Maxrep_Pruned_Updater updater;
    assert(score_string(S, G2, scorer, updater) == score_string_entropy_brute(S,T,threshold,escape_prob));
}

void test_RLE(){
    cerr << "Running RLE tests" << endl;
    
    srand(1412341);
    string T = get_random_string(300,3);
    string S = get_random_string(300,3);
    double threshold = 0.2;
    double escape_prob = 0.05;
         
    SLT_Iterator slt_it;
    Rev_ST_Maxrep_Iterator rev_slt_it;
    Entropy_Formula formula(threshold);
    Global_Data G_RLE;
    build_model(G_RLE, T, formula, slt_it, rev_slt_it, slt_it, true, false);
    Global_Data G_non_RLE;
    build_model(G_non_RLE, T, formula, slt_it, rev_slt_it, slt_it, false, false);
    
    Basic_Scorer scorer(escape_prob, true);
    Maxrep_Pruned_Updater updater;
    assert(score_string(S, G_RLE, scorer, updater) == score_string(S, G_non_RLE, scorer, updater));
}

void score_string_random_tests(int64_t number){
    cerr << "Running random score string tests for all context types" << endl;
    srand(1231231290);
    for(int64_t i = 0; i < number; i++){
        string S = get_random_string(100,3);
        string T = get_random_string(100,3);
        //string T = "cacabcbaacbbcaaabacbaabbaaacabccbaaacbaabccbcbabbccaacbbabbaababcbccbcacbbbbcbcbbabcabaaacbbbccabaaa"; // DEBUG
        //cout << "Testing " << T << endl;
        double threshold = rand() / (double)RAND_MAX;
        double escape = rand() / (double)RAND_MAX;
        //cout << T << endl;
                
        test_entropy(S, T, threshold, escape);
        
        // Try depth-bounded entropy scoring
        {
            int64_t depth_bound = 3;
            vector<string> contexts_vec = get_contexts_entropy_brute(T, threshold);
            set<string> contexts(contexts_vec.begin(), contexts_vec.end());
            double brute = score_string_brute_depth_bounded_given_contexts(S,T,contexts,escape,depth_bound);

            
            typedef Depth_Bounded_SLT_Iterator slt_it_t;
            typedef Rev_ST_Depth_Bounded_Maxrep_Iterator rev_slt_it_t;
            slt_it_t slt_it(depth_bound);
            rev_slt_it_t rev_slt_it(depth_bound);
            
            slt_it_t context_it(depth_bound);
            Entropy_Formula formula(threshold);
            
            Global_Data G; 
            build_model(G, T, formula, slt_it, rev_slt_it, context_it, true, false);
                                
            Basic_Scorer scorer(escape, true);
            Maxrep_Depth_Bounded_Updater updater;

            double non_brute = score_string(S, G, scorer, updater);
            assert(fabs(brute - non_brute) < 1e-6);
        }
        
        
        // Try depth-bounded KL scoring
        {
            int64_t depth_bound = 3;
            double threshold = rand() / (double)RAND_MAX;
            double escape = rand() / (double)RAND_MAX;
            vector<string> contexts_vec = get_contexts_KL_brute(T, threshold);
            set<string> contexts(contexts_vec.begin(), contexts_vec.end());
            double brute = score_string_brute_depth_bounded_given_contexts(S,T,contexts,escape,depth_bound);
            
            typedef Depth_Bounded_SLT_Iterator slt_it_t;
            typedef Rev_ST_Depth_Bounded_Maxrep_Iterator rev_slt_it_t;
            slt_it_t slt_it(depth_bound);
            rev_slt_it_t rev_slt_it(depth_bound);
            
            slt_it_t context_it(depth_bound-1); // -1 because contexts are of form aW
            KL_Formula formula(threshold);
            
            Global_Data G;
            build_model(G, T, formula, slt_it, rev_slt_it, context_it, true, false);
                
            int64_t asd = 0; // debug start
            for(string C : contexts) if(C.size() <= depth_bound) asd++;
                
            Basic_Scorer scorer(escape, false);
            Maxrep_Depth_Bounded_Updater updater;

            double non_brute = score_string(S, G, scorer, updater);
            assert(fabs(brute - non_brute) < 1e-6);
        }
                
        // Try formula 2,3,4 scoring
        {
            double t1 = ((double)rand() / RAND_MAX) / 5;
            double t2 = ((double)rand() / RAND_MAX);
            double t3 = 1 - ((double)rand() / RAND_MAX) / 5;
            double t4 = 1 + ((double)rand() / RAND_MAX) / 5;
            double escape = rand() / (double)RAND_MAX;
            double brute = score_string_eq_234_brute(S,T,t1,t2,t3,t4,escape);

            typedef SLT_Iterator slt_it_t;
            typedef Rev_ST_Maxrep_Iterator rev_slt_it_t;
            
            slt_it_t slt_it;
            rev_slt_it_t rev_slt_it;
            EQ234_Formula formula(t1,t2,t3,t4);

            
            Global_Data G; 
            build_model(G, T, formula, slt_it, rev_slt_it, slt_it, true, false);

            Basic_Scorer scorer(escape, false);
            Maxrep_Pruned_Updater updater;
            double non_brute = score_string(S, G, scorer, updater);
            //double non_brute = score_string_new2(S,T,NONE,formula,true,escape);
            
            assert(fabs(brute - non_brute) < 1e-6);
            //cerr << "Random string formula234 scoring test OK (" << t1 << " " << t2 << " " << t3 << " " << t4 <<
            //" " << escape << ")" << endl;
        }
        
        // Try KL scoring
        {
            double threshold = rand() / (double)RAND_MAX;
            double escape = rand() / (double)RAND_MAX;
            double brute = score_string_KL_brute(S,T,threshold,escape);
            
            typedef SLT_Iterator slt_it_t;
            typedef Rev_ST_Maxrep_Iterator rev_slt_it_t;
            
            slt_it_t slt_it;
            rev_slt_it_t rev_slt_it;
            
            KL_Formula formula(threshold);
            
            Global_Data G; 
            build_model(G, T, formula, slt_it, rev_slt_it, slt_it, true, false);

            Basic_Scorer scorer(escape, false);
            Maxrep_Pruned_Updater updater;
            double non_brute = score_string(S, G, scorer, updater);
            
            assert(fabs(brute - non_brute) < 1e-6);
            //cerr << "Random string KL scoring test OK. Threshold: " << threshold << ", escape: " << escape << endl;
        }
        
        // Try p-norm scoring
        {
            double threshold = rand() / (double)RAND_MAX;
            double escape = rand() / (double)RAND_MAX;
            int64_t p = 1 + (rand() % 5);
            double brute = score_string_p_norm_brute(S,T,p,threshold,escape);
            
            typedef SLT_Iterator slt_it_t;
            typedef Rev_ST_Maxrep_Iterator rev_slt_it_t;

            slt_it_t slt_it;
            rev_slt_it_t rev_slt_it;
            
            pnorm_Formula formula(p,threshold);
            
            Global_Data G;
            build_model(G, T, formula, slt_it, rev_slt_it, slt_it, true, false);

            Basic_Scorer scorer(escape, false);
            Maxrep_Pruned_Updater updater;
            double non_brute = score_string(S, G, scorer, updater);
            
            assert(fabs(brute - non_brute) < 1e-6);
            //cerr << "Random string p-norm scoring test OK. Threshold: " << threshold << ", escape: " << escape << endl;
        }
        

    }
}

void score_string_tests(){
    //score_string_testcase1();
    score_string_random_tests(100);
    
}

#endif
