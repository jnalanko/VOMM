#ifndef PRECALC_HH
#define PRECALC_HH

#include <iostream>
#include "BD_BWT_index/include/BD_BWT_index.hh"
#include "BD_BWT_index/include/Iterators.hh"
#include "sdsl/rank_support_v.hpp"
#include "sdsl/select_support_mcl.hpp"
#include "sdsl/bp_support_g.hpp"
#include "sdsl/util.hpp"
#include "Maxreps.hh"
#include "BWT_iteration.hh"
#include "Interfaces.hh"
#include <string>
#include <set>
#include <algorithm>
#include <numeric>
#include <vector>

/*
 * The BWT iterators in these functions are supposed to iterate all nodes that will
 * be included in the topology
 */

vector<bool> counters_to_bpr(vector<int64_t>& counters_open, vector<int64_t>& counters_close, bool fill_in_runs){
    
    int64_t last_pos = 0;
    int64_t last_type = -1;
    vector<bool> bpr;
    
    // Fill in artificial nodes so that every compact range of pruned nodes is
    // represented by a node
    
    if(fill_in_runs){ // Todo: delete
        for(int64_t i = 0; i < counters_close.size(); i++){
            for(int64_t j = 0; j < counters_open[i]; j++){
                if(last_type == 1 && i - last_pos >= 1){
                    counters_open[last_pos]++;
                    counters_close[i-1]++;
                }
                last_type = 1;
                last_pos = i;
            }
            for(int64_t j = 0; j < counters_close[i]; j++){
                if(last_type == 0 && i - last_pos >= 1){
                    counters_open[last_pos]++;
                    counters_close[i-1]++;
                }
                last_type = 0;
                last_pos = i;
            }
        }
    }

    // Build the BPR
    for(int64_t i = 0; i < counters_close.size(); i++){
        for(int64_t j = 0; j < counters_open[i]; j++){
            bpr.push_back(1);
        }
        for(int64_t j = 0; j < counters_close[i]; j++){
            bpr.push_back(0);
        }
    }
    
    return bpr;
}

sdsl::bit_vector get_slt_topology(BIBWT& index, Iterator& iterator){
    
    vector<int64_t> counters_open(index.size());
    vector<int64_t> counters_close(index.size());
    
    iterator.init();
    while(iterator.next()){
        counters_open[iterator.get_top().intervals.reverse.left]++;
        counters_close[iterator.get_top().intervals.reverse.right]++;
    }
    
    vector<bool> topology = counters_to_bpr(counters_open, counters_close, false);
    sdsl::bit_vector topology_sdsl(topology.size());
    for(int64_t i = 0; i < topology.size(); i++)
        topology_sdsl[i] = topology[i];
    return topology_sdsl;
}

sdsl::bit_vector get_slt_topology(BIBWT& index){
    
    SLT_Iterator iterator(&index);
    return get_slt_topology(index, iterator);
}


struct Rev_st_topology{
    sdsl::bit_vector bpr;
    sdsl::bit_vector pruning_marks;
};

void fill_in_counters(vector<int64_t>& counters_open, vector<int64_t>& counters_close){
    // Fill in artificial nodes so that every compact range of pruned nodes is
    // represented by a node. i.e. in between every ( (, every ) ) and every ) (, add a node
    
    int64_t n = counters_open.size();
    int64_t last_pos = 0;
    int64_t last_type = -1;
    int64_t debug = 0;
    
    for(int64_t i = 0; i < n; i++){
        if(counters_open[i] > 0){ // (
            if(last_type == 1 && i - last_pos >= 1){ // case ( (, gap size >= 1
                counters_open[last_pos]++;
                counters_close[last_pos]++;
                debug++;
            }
            if(last_type == 0 && i - last_pos >= 2){ // case ) (, gap size >= 2
                counters_open[last_pos+1]++;
                counters_close[last_pos+1]++;
                debug++;
            }
            last_type = 1;
            last_pos = i;
        }
        if(counters_close[i] > 0){ // )
            if(last_type == 0 && i - last_pos >= 1){ // case ) ), gap size >= 1
                counters_open[last_pos+1]++;
                counters_close[last_pos+1]++;
                debug++;
            }
            last_type = 0;
            last_pos = i;
        }
    }
    
}

Rev_st_topology get_rev_st_bpr_and_pruning(BIBWT& index, Iterator& iterator){
    int64_t n = index.size(); // index.size()?
    vector<int64_t> counters_open(n,0);
    vector<int64_t> counters_close(n,0);
    iterator.init();
    while(iterator.next()){
        Interval_pair I = iterator.get_top().intervals;
        counters_open[I.reverse.left]++;
        counters_close[I.reverse.right]++;
    }
    

    vector<bool> bpr; // don't know the length yet
    sdsl::bit_vector pruning_marks(n,0);

    //fill_in_counters(counters_open, counters_close); // Don't need because we always have maxrep left extensions in
    
    // Build the BPR and pruning marks
    for(int64_t i = 0; i < n; i++){
        for(int64_t j = 0; j < counters_open[i]; j++){
            bpr.push_back(1);
        }
        for(int64_t j = 0; j < counters_close[i]; j++){
            bpr.push_back(0);
        }
        
        if(counters_open[i] > 0){
            pruning_marks[i] = 1;
        }
    }
    
    // Into a sdsl bit vector
    sdsl::bit_vector bpr_sdsl(bpr.size());
    for(int64_t i = 0; i < bpr.size(); i++)
        bpr_sdsl[i] = bpr[i];
    
    return {bpr_sdsl, pruning_marks};
}

sdsl::bit_vector get_rev_st_topology(BIBWT& index){
    Rev_ST_Iterator iterator(&index);
    return get_rev_st_bpr_and_pruning(index, iterator).bpr;
}

template <typename topology_mapper_t>
sdsl::bit_vector get_rev_st_maximal_marks(BIBWT& index, int64_t rev_st_bpr_length, Iterator& it, topology_mapper_t& mapper){
    sdsl::bit_vector marks(rev_st_bpr_length,0);
    it.init();
    while(it.next()){
        Interval_pair I = it.get_top().intervals;
        if(index.is_left_maximal(I) && index.is_right_maximal(I)){
            marks[mapper.leaves_to_node(I.reverse)] = 1;
        }
    }
    return marks;
}

template <typename topology_mapper_t>
sdsl::bit_vector get_rev_st_maximal_marks(BIBWT& index, int64_t rev_st_bpr_length, topology_mapper_t& mapper){
    Rev_ST_Iterator iterator(&index);
    return get_rev_st_maximal_marks(index, rev_st_bpr_length, iterator, mapper);
}

sdsl::bit_vector get_slt_maximal_marks(BIBWT& index, sdsl::bit_vector& slt_bpr, Iterator& it){
    sdsl::bit_vector marks(slt_bpr.size(),0);
    
    sdsl::select_support_mcl<1> slt_bpr_ss;
    sdsl::util::init_support(slt_bpr_ss, &slt_bpr);
    
    int64_t preorder_rank = 0;
    
    it.init();
    while(it.next()){
        preorder_rank++;
        Interval_pair I = it.get_top().intervals;
        if(index.is_left_maximal(I) && index.is_right_maximal(I)){
            marks[slt_bpr_ss.select(preorder_rank)] = 1;
        }
    }
    
    return marks;
}

sdsl::bit_vector get_slt_maximal_marks(BIBWT& index, sdsl::bit_vector& slt_bpr){
    SLT_Iterator iterator(&index);
    return get_slt_maximal_marks(index,slt_bpr,iterator);
}




// Returns pair rev_st_marks, slt_marks
// SLT iterator need to be left-to-right depth first order
template <typename topology_mapper_t>
std::pair<sdsl::bit_vector,sdsl::bit_vector> get_rev_st_and_slt_maximal_marks(
                                                    BIBWT& index,
                                                    int64_t rev_st_bpr_length,
                                                    Iterator& it,
                                                    topology_mapper_t& mapper,
                                                    sdsl::bit_vector& slt_bpr){
    
    sdsl::bit_vector marks_rev_st(rev_st_bpr_length,0);
    sdsl::bit_vector marks_slt(slt_bpr.size(),0);

    sdsl::select_support_mcl<1> slt_bpr_ss;
    sdsl::util::init_support(slt_bpr_ss, &slt_bpr);
    
    int64_t preorder_rank = 0;

    it.init();
    while(it.next()){
        preorder_rank++;
        Interval_pair I = it.get_top().intervals;
        if(index.is_left_maximal(I) && index.is_right_maximal(I)){
            marks_rev_st[mapper.leaves_to_node(I.reverse)] = 1;
            marks_slt[slt_bpr_ss.select(preorder_rank)] = 1;
        }
    }
    
    return {marks_rev_st, marks_slt};
}

// Returns pair rev_st_marks, slt_marks
template <typename topology_mapper_t>
std::pair<sdsl::bit_vector,sdsl::bit_vector> get_rev_st_and_slt_maximal_marks(
                                                    BIBWT& index,
                                                    int64_t rev_st_bpr_length,
                                                    topology_mapper_t& mapper,
                                                    sdsl::bit_vector& slt_bpr){
    
    SLT_Iterator iterator(&index);
    return get_rev_st_and_slt_maximal_marks(index,rev_st_bpr_length,iterator,mapper,slt_bpr);
}

/*
// Marks the leftmost leaf every iterated node
template <typename bwt_t, typename rev_slt_iterator_t>
sdsl::bit_vector mark_rev_st_pruning(bwt_t& index, rev_slt_iterator_t it){
    sdsl::bit_vector marks(index.size(), 0);
    assert(marks.size() > 0);
    it.init();
    while(it.next()){
        marks[it.top.intervals.reverse.left] = 1;
    }
    return marks;
} */

template <typename bwt_t, typename mapper_t>
std::vector<int64_t> get_all_rev_st_string_depths(bwt_t& index, int64_t rev_st_bpr_length, Iterator& it, mapper_t mapper){
    (void)index;
    vector<int64_t> v(rev_st_bpr_length);
    it.init();
    while(it.next()){
        int64_t depth = it.get_top().depth;
        assert(depth > 0); // Full iterator does not know depths and so returns -1 for all.
        v[mapper.leaves_to_node(it.get_top().intervals.reverse)] = depth;
    }
    return v;
}

#endif
