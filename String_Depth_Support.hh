#ifndef STRING_DEPTH_SUPPORT_HH
#define STRING_DEPTH_SUPPORT_HH

#include "BPR_tools.hh"
#include <iostream>
#include "BD_BWT_index/include/BD_BWT_index.hh"
#include "Interfaces.hh"
#include "sdsl/rank_support_v.hpp"
#include "sdsl/select_support_mcl.hpp"
#include "sdsl/bp_support_g.hpp"
#include "sdsl/util.hpp"
#include "globals.hh"
#include <string>
#include <set>
#include <algorithm>
#include <numeric>
#include <vector>


using namespace std;

// 1 represents an open parenthesis, 0 a closed parenthesis
// This class should contain only pointer data members, so it can be copied easily
class String_Depth_Support{
    public:
    String_Depth_Support() {};
    String_Depth_Support(std::shared_ptr<Bitvector> rev_st_bpr, 
                         std::shared_ptr<Bitvector> slt_bpr, 
                         std::shared_ptr<Bitvector> rev_st_maximal_marks,
                         std::shared_ptr<Bitvector> slt_maximal_marks)
      : rev_st_bpr(rev_st_bpr), slt_bpr(slt_bpr), rev_st_maximal_marks(rev_st_maximal_marks), slt_maximal_marks(slt_maximal_marks) {}; // Assuming these have all the required support structures

    public:

    // Utilities
   
    std::shared_ptr<Bitvector> rev_st_bpr;
    std::shared_ptr<Bitvector> slt_bpr;
    std::shared_ptr<Bitvector> rev_st_maximal_marks;
    std::shared_ptr<Bitvector> slt_maximal_marks;
    
    int64_t string_depth(int64_t open){
        return slt_tree_depth(preorder_rank_to_slt_bpr(rev_st_bpr_to_preorder_rank(open)));
    }

    int64_t rev_st_bpr_to_preorder_rank(int64_t open){
        // Rank opening parenthesis withing marked nodes
        return rev_st_maximal_marks->rank(open+1) - 1; // -1 to make indexing start from zero
    }

    int64_t preorder_rank_to_slt_bpr(int64_t rank){
        // Select opening parenthesis within marked nodes
        int64_t open = slt_maximal_marks->select(rank+1); // select indexing is from 1
        return open;
    }

    int64_t slt_tree_depth(int64_t open){
        int64_t ones = slt_bpr->rank(open+1);
        int64_t zeros = open + 1 - ones;
        return ones - zeros - 1; // -1: root has depth zero
    }
        
};

class String_Depth_Support_Store_All{
    
    public:
    
    std::vector<int64_t>* depths;
    
    String_Depth_Support_Store_All() {}
    String_Depth_Support_Store_All(vector<int64_t>* depths) : depths(depths) {}
    
    int64_t string_depth(int64_t open_paren){
        return (*depths)[open_paren];
    }
};



#endif
