#ifndef STRING_DEPTH_SUPPORT_HH
#define STRING_DEPTH_SUPPORT_HH

#include "BPR_tools.hh"
#include <iostream>
#include "BD_BWT_index/include/BD_BWT_index.hh"
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

// todo: int64_t -> size_t

// 1 represents an open parenthesis, 0 a closed parenthesis

using namespace std;

// Does this:
// Reverse ST interval of a maximal node -> string depth
//
// In more detail
// ST interval -> ST BPR interval -> ST preorder rank -> SLT preorder rank -> SLT BPR interval -> SLT node depth

//  ST interval             
//       |                         
//       v                         
//  ST BPR interval         SLT BPR interval  ---> SLT node depth --> ST string depth
//       |                         ^
//       v                         |
//  ST preorder rank  --->  SLT preorder rank
//
// This class should contain only pointer data members, so it can be copied easily
class String_Depth_Support{
    public:
    String_Depth_Support() {};
    String_Depth_Support(
        sdsl::select_support_mcl<10,2>* rev_st_ss_10,
        sdsl::bp_support_g<>* rev_st_bps,
        sdsl::rank_support_v<1>* rev_st_maximal_marks_rs,
        sdsl::select_support_mcl<1>* slt_maximal_marks_ss,
        sdsl::bp_support_g<>* slt_bps);
    String_Depth_Support(sdsl::bit_vector& rev_st_bpr,
                         sdsl::bit_vector& slt_bpr,
                         sdsl::bit_vector& rev_st_maximal_marks,
                         sdsl::bit_vector& slt_maximal_marks);
    
    // The most important function: string_depth
    int64_t string_depth(Interval I); // Lex interval of a maximal interval the reverse st -> string depth
    int64_t string_depth(int64_t open_paren);

    Interval rev_st_interval_to_bpr(Interval I); // Map left and right into BPR and take LCA = select '()' in st bpr and enclose
    int64_t rev_st_bpr_to_preorder_rank(Interval I); // Rank opening parenthesis wihtin marked nodes
    Interval preorder_rank_to_slt_bpr(int64_t rank); // Select opening parenthesis within marked nodes
    int64_t slt_tree_depth(Interval bpr_node); // Excess in bpr

    public:

    // Utilities
    Interval get_rev_st_leaf_bpr(int64_t leaf_rank);
    Interval LCA(Interval A, Interval B); // In ST
   
    // Support structures for the suffix tree
    sdsl::select_support_mcl<10,2>* rev_st_ss_10; // find the i-th leaf in the bpr
    sdsl::bp_support_g<>* rev_st_bps; // enclose
    sdsl::rank_support_v<1>* rev_st_maximal_marks_rs; // for preorder ranks

    // Support structures for the suffix link tree
    sdsl::select_support_mcl<1>* slt_maximal_marks_ss; // for preorder selects
    sdsl::bp_support_g<>* slt_bps; // for node depth in slt

};

int64_t String_Depth_Support::string_depth(Interval I){
    return slt_tree_depth(preorder_rank_to_slt_bpr(rev_st_bpr_to_preorder_rank(rev_st_interval_to_bpr(I))));
}

int64_t String_Depth_Support::string_depth(int64_t open_paren){
    int64_t close = rev_st_bps->find_close(open_paren);
    return slt_tree_depth(
               preorder_rank_to_slt_bpr(
                   rev_st_bpr_to_preorder_rank(
                       Interval(open_paren,close)
                   )
               )
           );
}

int64_t String_Depth_Support::rev_st_bpr_to_preorder_rank(Interval I){
    // Rank opening parenthesis withing marked nodes
    return rev_st_maximal_marks_rs->rank(I.left+1) - 1; // -1 to make indexing start from zero
}

Interval String_Depth_Support::preorder_rank_to_slt_bpr(int64_t rank){
    // Select opening parenthesis within marked nodes
    int64_t open = slt_maximal_marks_ss->select(rank+1); // select indexing is from 1
    int64_t close = slt_bps->find_close(open);
    return Interval(open, close);
}

int64_t String_Depth_Support::slt_tree_depth(Interval bpr_node){
    return slt_bps->excess(bpr_node.left) - 1; // minus one to set the root at depth zero
}


String_Depth_Support::String_Depth_Support(
    sdsl::select_support_mcl<10,2>* rev_st_ss_10,
    sdsl::bp_support_g<>* rev_st_bps,
    sdsl::rank_support_v<1>* rev_st_maximal_marks_rs,
    sdsl::select_support_mcl<1>* slt_maximal_marks_ss,
    sdsl::bp_support_g<>* slt_bps) : 
        rev_st_ss_10(rev_st_ss_10),
        rev_st_bps(rev_st_bps),
        rev_st_maximal_marks_rs(rev_st_maximal_marks_rs),
        slt_maximal_marks_ss(slt_maximal_marks_ss),
        slt_bps(slt_bps)
    {}

// Constructor for the lazy (and for tests)
String_Depth_Support::String_Depth_Support(sdsl::bit_vector& rev_st_bpr,
                                           sdsl::bit_vector& slt_bpr,
                                           sdsl::bit_vector& rev_st_maximal_marks,
                                           sdsl::bit_vector& slt_maximal_marks){
                                               
    // TODO: These leak memory. Can't just delete in destructor because of the other constructor.
    this->rev_st_ss_10 = new sdsl::select_support_mcl<10,2>();
    this->rev_st_bps = new sdsl::bp_support_g<>();
    this->rev_st_maximal_marks_rs = new sdsl::rank_support_v<1>();
    this->slt_maximal_marks_ss = new sdsl::select_support_mcl<1>();
    this->slt_bps = new sdsl::bp_support_g<>();
    
    sdsl::util::init_support(*this->rev_st_maximal_marks_rs, &rev_st_maximal_marks);
    sdsl::util::init_support(*this->slt_maximal_marks_ss, &slt_maximal_marks);
    sdsl::util::init_support(*this->rev_st_ss_10, &rev_st_bpr);
    sdsl::util::init_support(*this->rev_st_bps, &rev_st_bpr);
    sdsl::util::init_support(*this->slt_bps, &slt_bpr);

}


Interval String_Depth_Support::get_rev_st_leaf_bpr(int64_t leaf_rank){
    int64_t close = rev_st_ss_10->select(leaf_rank+1); // Indexing starts from 1, hence the +1
    return Interval(close-1,close);
}


Interval String_Depth_Support::rev_st_interval_to_bpr(Interval I){
    if(I.left == I.right) return get_rev_st_leaf_bpr(I.left);
    Interval L = get_rev_st_leaf_bpr(I.left);
    Interval R = get_rev_st_leaf_bpr(I.right);
    return LCA(L,R);
}

// A must be to the left of B and B can not be nested inside A
Interval String_Depth_Support::LCA(Interval A, Interval B){
    assert(A.right < B.left);
    int64_t open = rev_st_bps->double_enclose(A.left, B.left);
    int64_t close = rev_st_bps->find_close(open);
    return Interval(open,close);
}

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
