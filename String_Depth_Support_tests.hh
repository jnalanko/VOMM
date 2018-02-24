#ifndef String_Depth_Support_tests_HH
#define String_Depth_Support_tests_HH


#include "brute_tools.hh"
#include "String_Depth_Support.hh"
#include "Maxreps.hh"
#include "Precalc.hh"
#include "context_marking.hh"
#include "globals.hh"
#include "score_string.hh"
#include <utility>

class String_Depth_Support_Tester{
    
public:
    
    void test_rev_st_bpr(){
        {
            string text = "mississippi"; // Note that this is a bad testcase, because its topology is for the reverse ST and the SDL
            BD_BWT_index<> index((uint8_t*)text.c_str());
            sdsl::bit_vector rev_st_bpr = get_rev_st_topology(index);
            
            sdsl::bit_vector rev_st_topology_correct = {1,1,0,1,1,0,1,0,1,1,0,1,0,0,0,1,0,1,1,0,1,0,0,1,1,1,0,1,0,0,1,1,0,1,0,0,0,0};
            assert(rev_st_bpr == rev_st_topology_correct);
        }
        
        {
            string text = "abracabra";
            BD_BWT_index<> index((uint8_t*)text.c_str());
            sdsl::bit_vector rev_st_bpr = get_rev_st_topology(index);
            
            sdsl::bit_vector rev_st_topology_correct = {1,1,0,1,1,0,1,0,1,1,0,1,0,0,0,1,1,0,1,0,0,1,0,1,1,0,1,0,0,0}; // st of reverse
            assert(rev_st_bpr == rev_st_topology_correct);
        }
        
        cerr << "ST Topology building test OK" << endl;
    }
    
    
    
    void test_slt_bpr(){
        string text = "mississippi"; // Note that this is a bad testcase, because its topology is same for the reverse ST and the SDL
        BD_BWT_index<> index((uint8_t*)text.c_str());
        sdsl::bit_vector slt_bpr = get_slt_topology(index);
        
        sdsl::bit_vector slt_topology_correct = {1,1,1,1,1,0,0,0,0,1,0,1,0,0};
        assert(slt_bpr == slt_topology_correct);
        cerr << "SLT Topology building test OK" << endl;
    }
    
    /*
     
     void test_lca(string text){
     BD_BWT_index<> index((uint8_t*)text.c_str());
     String_Depth_Support SDS(index);
     
     // Test LCA for all pairs of leaves
     for(int64_t leaf_A = 0; leaf_A < index.size(); leaf_A++){
     for(int64_t leaf_B = leaf_A+1; leaf_B < index.size(); leaf_B++){
     Interval A = SDS.get_st_leaf_bpr(leaf_A);
     Interval B = SDS.get_st_leaf_bpr(leaf_B);
     Interval LCA_bps = SDS.LCA(A,B);
     Interval LCA_correct = LCA_naive(SDS.st_bpr,A,B);
     assert(LCA_bps == LCA_correct);
     }
     }
     cerr << "LCA test OK" << endl;
     }
     
     */
    
    
    void test_mark_maximal_both(){
        string input = "mississippi";
        BD_BWT_index<> index((uint8_t*)input.c_str());
        
        sdsl::bit_vector slt_bpr = get_slt_topology(index);
        sdsl::bit_vector rev_st_bpr = get_rev_st_topology(index);
        
        sdsl::select_support_mcl<10,2> rev_st_ss_10; // find the i-th leaf in the bpr
        sdsl::rank_support_v<10,2> rev_st_rs_10;
        sdsl::bp_support_g<> rev_st_bps; // enclose
        
        sdsl::util::init_support(rev_st_ss_10, &rev_st_bpr);
        sdsl::util::init_support(rev_st_rs_10, &rev_st_bpr);
        sdsl::util::init_support(rev_st_bps, &rev_st_bpr);
        
        Full_Topology_Mapper mapper(&rev_st_bps, &rev_st_ss_10, &rev_st_rs_10);
        
        pair<sdsl::bit_vector,sdsl::bit_vector> marks
        = get_rev_st_and_slt_maximal_marks(index,
                                           rev_st_bpr.size(),
                                           mapper,
                                           slt_bpr);
        
        sdsl::bit_vector rev_st_maximal_marks = marks.first;
        sdsl::bit_vector slt_maximal_marks = marks.second;
        
        sdsl::bit_vector rev_st_correct = {1,0,0,1,0,0,0,0,1,0,0,0,0,0,0,0,0,1,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
        sdsl::bit_vector slt_correct = {1,1,0,0,1,0,0,0,0,1,0,1,0,0};

        assert(rev_st_maximal_marks == rev_st_correct);
        assert(slt_maximal_marks == slt_correct);
    }
    
    /*sdsl::bit_vector mark_contexts_entropy(BD_BWT_index<>& index, int64_t rev_st_bpr_length, double threshold,
     sdsl::select_support_mcl<10,2>& rev_st_ss_10,
     sdsl::bp_support_g<>& rev_st_bps){*/
    void test_mark_contexts_entropy(string text){
        BD_BWT_index<> index((uint8_t*)text.c_str());
        sdsl::bit_vector rev_st_bpr = get_rev_st_topology(index);
        sdsl::select_support_mcl<10,2> rev_st_ss_10(&rev_st_bpr);
        sdsl::rank_support_v<10,2> rev_st_rs_10(&rev_st_bpr);
        sdsl::bp_support_g<> rev_st_bps(&rev_st_bpr);
        Full_Topology_Mapper mapper(&rev_st_bps, &rev_st_ss_10, &rev_st_rs_10);
        SLT_Iterator iterator(&index);
        
        Entropy_Formula F(0.5);
        sdsl::bit_vector contexts = F.get_rev_st_context_marks(&index, rev_st_bpr.size(), iterator, mapper);
        //sdsl::bit_vector contexts = mark_contexts_entropy(index, rev_st_bpr.size(), 0.5, rev_st_ss_10, rev_st_bps);
        
        // Check that the contexts are a subset of the maximal repeats
        sdsl::bit_vector maximals = get_rev_st_maximal_marks(index, rev_st_bpr.size(),mapper);
        
        for(int64_t i = 0; i < rev_st_bpr.size(); i++){
            // Contexts must be a subset of maximal repeats
            // In the contexts vector, both open and close parenthesis are marked
            // In the maximal marks vector, only open parentheses are marked
            if(contexts[i] && rev_st_bpr[i] == 1) assert(maximals[i]);
        }
        
        cout << contexts << endl;
    }
    
    void test_string_depth(string input){
        
        // Init
        
        BD_BWT_index<> index((uint8_t*)input.c_str());
        
        sdsl::bit_vector slt_bpr = get_slt_topology(index);
        sdsl::bit_vector rev_st_bpr = get_rev_st_topology(index);
        
        sdsl::select_support_mcl<10,2> rev_st_ss_10; // find the i-th leaf in the bpr
        sdsl::rank_support_v<10,2> rev_st_rs_10;
        sdsl::bp_support_g<> rev_st_bps; // enclose
        
        sdsl::util::init_support(rev_st_ss_10, &rev_st_bpr);
        sdsl::util::init_support(rev_st_rs_10, &rev_st_bpr);
        sdsl::util::init_support(rev_st_bps, &rev_st_bpr);
        Full_Topology_Mapper mapper(&rev_st_bps, &rev_st_ss_10, &rev_st_rs_10);
        
        pair<sdsl::bit_vector,sdsl::bit_vector> marks
        = get_rev_st_and_slt_maximal_marks(index,
                                           rev_st_bpr.size(),
                                           mapper,
                                           slt_bpr);
        
        sdsl::bit_vector rev_st_maximal_marks = marks.first;
        sdsl::bit_vector slt_maximal_marks = marks.second;
        
        String_Depth_Support SDS(rev_st_bpr,
                                 slt_bpr,
                                 rev_st_maximal_marks,
                                 slt_maximal_marks);
        
        // Do the testing
        
        std::vector<std::pair<std::string, Interval_pair> > maxreps = find_maxreps(index); // For reference
        
        int64_t nmarked = SDS.rev_st_maximal_marks_rs->rank(rev_st_bpr.size());
        
        assert(nmarked == maxreps.size());
        
        for(auto rep : maxreps){
            string label = rep.first;
            Interval_pair IP = rep.second;
            assert(SDS.string_depth(IP.reverse) == label.size());
        }
        
        //cout << "String depth test OK" << endl;
    }
    
public:
    
    
    Interval parent_naive(sdsl::bit_vector& v, Interval A){
        if(A.left == 0) return A; // Root
        Interval P(A.left, A.right);
        
        int64_t left_excess = 0;
        while(left_excess >= 0){
            P.left--;
            left_excess += (v[P.left] == 0) ? 1 : -1;
        }
        
        int64_t right_excess = 0;
        while(right_excess >= 0){
            P.right++;
            right_excess += (v[P.right] == 1) ? 1 : -1;
        }
        
        return P;
    }
    
    Interval LCA_naive(sdsl::bit_vector& v, Interval A, Interval B){
        assert(A.right < B.left);
        vector<Interval> A_parents;
        do{
            A_parents.push_back(A);
            A = parent_naive(v,A);
        } while(parent_naive(v,A) != A_parents.back());
        
        vector<Interval> B_parents;
        do{
            B_parents.push_back(B);
            B = parent_naive(v,B);
        } while(parent_naive(v,B) != B_parents.back());
        
        Interval ans;
        while(A_parents.back() == B_parents.back()){
            ans = A_parents.back();
            A_parents.pop_back();
            B_parents.pop_back();
        }
        
        return ans;
    }
    
    int64_t tree_depth_naive(sdsl::bit_vector v, Interval A){
        int64_t depth = 0;
        while(parent_naive(v,A) != A){ // While A is not the root
            depth++;
            A = parent_naive(v,A);
        }
        return depth;
    }
    
    int64_t rank_naive(sdsl::bit_vector v, int64_t k){
        int64_t ans = 0;
        for(int64_t i = 0; i < k; i++){
            ans += v[i];
        }
        return ans;
    }
    
    /*
     void wtf(){
     BD_BWT_index<> index((uint8_t*)"abababbabbabbabaabbaabbabaaabababbbababbbabababababab");
     String_Depth_Support SDS(index);
     std::vector<std::pair<std::string, Interval_pair> > maxreps = find_maxreps(index);
     cout << "st marks " << SDS.st_maximal_marks << endl;
     cout << "slt marks " << SDS.slt_maximal_marks << endl;
     cout << "st bpr " << SDS.st_bpr << endl;
     cout << "slt bpr " << SDS.slt_bpr << endl;
     
     int64_t st_nmarked = 0;
     int64_t slt_nmarked = 0;
     for(int64_t i = 0; i < SDS.st_maximal_marks.size(); i++){
     st_nmarked += SDS.st_maximal_marks[i];
     }
     for(int64_t i = 0; i < SDS.slt_maximal_marks.size(); i++){
     slt_nmarked += SDS.slt_maximal_marks[i];
     }
     
     cout << "st nmarked " << st_nmarked << endl;
     cout << "slt nmarked " << slt_nmarked << endl;
     cout << "Actual number of maximal: " << maxreps.size() << endl;
     
     
     //assert(nmarked == maxreps.size());
     }
     
     */
    
};

void String_Depth_Support_tests(){
    
    srand(3248923);
    
    String_Depth_Support_Tester SDST;
    
    cerr << "Running rev st bpr building test" << endl;
    SDST.test_rev_st_bpr();
    cerr << "Running slt bpr building test" << endl;
    SDST.test_slt_bpr();
    cerr << "Running mark maximal both test" << endl;
    SDST.test_mark_maximal_both();
    SDST.test_mark_contexts_entropy("aaaaaaabaaaaabaaaaaaaabaaaaabaaaaabaaaabbaaaaaabaabaaaaaabaaabaaaabaaaabaaaabaaaaaabababaaaaaaa");
    
    cerr << "Running string depth support tests" << endl;
    SDST.test_string_depth("abababbabbabbabaabbaabbabaaabababbbababbbabababababab");
    SDST.test_string_depth("abbabbabaabbaabbabaaabababbbababbbababa");
    SDST.test_string_depth("aaaaaaabaaaaabaaaaaaaabaaaaabaaaaabaaaabbaaaaaabaabaaaaaabaaabaaaabaaaabaaaabaaaaaabababaaaaaaa");
    SDST.test_string_depth("abracabra");
    for(int64_t i = 0; i < 100; i++){
        SDST.test_string_depth(get_random_string(1000,3));
    }
    //cout << "All tests OK" << endl;
}

#endif

