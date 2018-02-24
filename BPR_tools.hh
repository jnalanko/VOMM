#ifndef BPR_TOOLS_HH
#define BPR_TOOLS_HH

#include <cassert>
#include "BD_BWT_index/include/BD_BWT_index.hh"
#include "sdsl/rank_support_v.hpp"
#include "sdsl/select_support_mcl.hpp"
#include "sdsl/bp_support_g.hpp"
#include "sdsl/util.hpp"
#include <iostream>
 
Interval LCA(Interval A, Interval B, sdsl::bp_support_g<>& bps){
    assert(A.right < B.left);
    int64_t open = bps.double_enclose(A.left, B.left);
    int64_t close = bps.find_close(open);
    return Interval(open,close);
}

Interval find_leaf_in_bpr(int64_t leaf_rank, sdsl::select_support_mcl<10,2>& ss_10){
    int64_t close = ss_10.select(leaf_rank+1); // Indexing starts from 1, hence the +1
    return Interval(close-1,close);
}

Interval enclose_leaves(int64_t leaf1, int64_t leaf2, sdsl::select_support_mcl<10,2>& ss_10, sdsl::bp_support_g<>& bps){
    if(leaf1 == leaf2) return find_leaf_in_bpr(leaf1,ss_10);
    Interval L = find_leaf_in_bpr(leaf1,ss_10);
    Interval R = find_leaf_in_bpr(leaf2,ss_10);
    return LCA(L,R,bps);
}

Interval bpr_interval_to_leaf_interval(Interval I, sdsl::rank_support_v<10,2>& bpr_rs_10){
    int64_t leaves_before = bpr_rs_10.rank(I.left);
    int64_t leaves_inside = bpr_rs_10.rank(I.right+1);
    return Interval(leaves_before, leaves_inside-1); //-1: make end of range inclusive
}

#endif
