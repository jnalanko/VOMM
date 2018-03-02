#ifndef String_Depth_Support_tests_HH
#define String_Depth_Support_tests_HH


#include "brute_tools.hh"
#include "BPR_tools.hh"
#include "String_Depth_Support.hh"
#include "Maxreps.hh"
#include "Precalc.hh"
#include "context_marking.hh"
#include "globals.hh"
#include "score_string.hh"
#include <utility>

class String_Depth_Support_Tester{
    
public:
    
    void test_mark_maximal_both(){
        string input = "mississippi";
        BD_BWT_index<> index((uint8_t*)input.c_str());
        
        sdsl::bit_vector rev_st_bpr_sdsl = get_rev_st_topology(index);
        sdsl::bit_vector slt_bpr_sdsl = get_slt_topology(index);
        
        std::shared_ptr<Basic_bitvector> rev_st_bpr = make_shared<Basic_bitvector>(rev_st_bpr_sdsl);
        rev_st_bpr->init_rank_10_support();
        rev_st_bpr->init_select_10_support();
        rev_st_bpr->init_bps_support();
        
        Full_Topology_Mapper mapper(rev_st_bpr);
        
        pair<sdsl::bit_vector,sdsl::bit_vector> marks
        = get_rev_st_and_slt_maximal_marks(index,
                                           rev_st_bpr->size(),
                                           mapper,
                                           slt_bpr_sdsl);
        
        sdsl::bit_vector rev_st_maximal_marks = marks.first;
        sdsl::bit_vector slt_maximal_marks = marks.second;
        
        sdsl::bit_vector rev_st_correct = {1,0,0,1,0,0,0,0,1,0,0,0,0,0,0,0,0,1,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
        sdsl::bit_vector slt_correct = {1,1,0,0,1,0,0,0,0,1,0,1,0,0};

        assert(rev_st_maximal_marks == rev_st_correct);
        assert(slt_maximal_marks == slt_correct);
    }
        
    void test_string_depth(string input){
        
        // Init
        
        BD_BWT_index<> index((uint8_t*)input.c_str());
        
        sdsl::bit_vector slt_bpr_sdsl = get_slt_topology(index);
        sdsl::bit_vector rev_st_bpr_sdsl = get_rev_st_topology(index);
        
        std::shared_ptr<Basic_bitvector> rev_st_bpr = make_shared<Basic_bitvector>(rev_st_bpr_sdsl);
        rev_st_bpr->init_rank_10_support();
        rev_st_bpr->init_select_10_support();
        rev_st_bpr->init_bps_support();
        
        std::shared_ptr<Basic_bitvector> slt_bpr = make_shared<Basic_bitvector>(slt_bpr_sdsl);
        slt_bpr->init_bps_support();
        slt_bpr->init_rank_support();
        
        Full_Topology_Mapper mapper(rev_st_bpr);
        
        pair<sdsl::bit_vector,sdsl::bit_vector> marks
        = get_rev_st_and_slt_maximal_marks(index,
                                           rev_st_bpr->size(),
                                           mapper,
                                           slt_bpr_sdsl);
        
        sdsl::bit_vector rev_st_maximal_marks_sdsl = marks.first;
        std::shared_ptr<Basic_bitvector> rev_st_maximal_marks = make_shared<Basic_bitvector>(rev_st_maximal_marks_sdsl);
        rev_st_maximal_marks->init_rank_support();
        
        sdsl::bit_vector slt_maximal_marks_sdsl = marks.second;
        std::shared_ptr<Basic_bitvector> slt_maximal_marks = make_shared<Basic_bitvector>(slt_maximal_marks_sdsl);
        slt_maximal_marks->init_select_support();
        
        String_Depth_Support SDS(rev_st_bpr,
                                 slt_bpr,
                                 rev_st_maximal_marks,
                                 slt_maximal_marks);
        
        // Do the testing
        
        std::vector<std::pair<std::string, Interval_pair> > maxreps = find_maxreps(index); // For reference
        
        int64_t nmarked = SDS.rev_st_maximal_marks->rank(rev_st_bpr->size());
        
        assert(nmarked == maxreps.size());
        
        for(auto rep : maxreps){
            string label = rep.first;
            Interval_pair IP = rep.second;
            int64_t open = enclose_leaves(IP.reverse.left, IP.reverse.right, rev_st_bpr->ss_10, rev_st_bpr->bps).left;
            assert(SDS.string_depth(open) == label.size());
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
    
};

void String_Depth_Support_tests(){
    
    srand(3248923);
    
    String_Depth_Support_Tester SDST;
    
    cerr << "Running mark maximal both test" << endl;
    SDST.test_mark_maximal_both();
    
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

