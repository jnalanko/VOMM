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
#include "Counters.hh"

/*
 * The BWT iterators in these functions are supposed to iterate all nodes that will
 * be included in the topology
 */


vector<bool> counters_to_bpr(Counters& counters){
    
    vector<bool> bpr;
    
    // Build the BPR
    for(int64_t i = 0; i < counters.size(); i++){
        for(int64_t j = 0; j < counters.get_open(i); j++){
            bpr.push_back(1);
        }
        for(int64_t j = 0; j < counters.get_close(i); j++){
            bpr.push_back(0);
        }
    }
    
    return bpr;
}

class Build_SLT_BPR_Callback : public Iterator_Callback{
  
public:
    
    Succinct_Counters counters;
    sdsl::bit_vector bpr_sdsl;
    bool enabled;
    
    Build_SLT_BPR_Callback() : enabled(true) {}
    
    void init(BIBWT& index){
        if(!enabled) return;
        counters.init(index.size());
    }
    
    void enable(){ enabled = true;}
    void disable(){ enabled = false;}
    
    void callback(const Iterator::Stack_frame& top){
        if(!enabled) return;
        counters.increment_open(top.intervals.reverse.left);
        counters.increment_close(top.intervals.reverse.right);
    }
    
    void finish(){
        if(!enabled) return;
        vector<bool> bpr = counters_to_bpr(counters);
        bpr_sdsl.resize(bpr.size());
        for(int64_t i = 0; i < bpr.size(); i++)
            bpr_sdsl[i] = bpr[i];
        counters.free_memory();
    }
    
    sdsl::bit_vector get_result(){
        return bpr_sdsl;
    }
};


struct Rev_st_topology{
    sdsl::bit_vector bpr;
    sdsl::bit_vector pruning_marks;
};

class Build_REV_ST_BPR_And_Pruning_Callback : public Iterator_Callback{
  
public:
    
    Succinct_Counters counters;
    sdsl::bit_vector bpr_sdsl;
    sdsl::bit_vector pruning;
    
    void init(BIBWT& index){
        counters.init(index.size());
        pruning = sdsl::bit_vector(index.size(),0);
        
        // Don't know size of bpr_sdsl yet
    }
    
    virtual void callback(const Iterator::Stack_frame& top){
        counters.increment_open(top.intervals.reverse.left);
        counters.increment_close(top.intervals.reverse.right);
    }
    
    virtual void finish(){
        vector<bool> bpr = counters_to_bpr(counters);
        
        // Copy to bpr_sdsl
        bpr_sdsl.resize(bpr.size());
        for(int64_t i = 0; i < bpr.size(); i++)
            bpr_sdsl[i] = bpr[i];
        
        // Compute pruning marks
        for(int64_t i = 0; i < counters.size(); i++){            
            if(counters.get_open(i) > 0){
                pruning[i] = 1;
            }
        }
        
        counters.free_memory();
    }
    
    Rev_st_topology get_result(){
        return {bpr_sdsl, pruning};
    }
};

class Rev_ST_Maximal_Marks_Callback : public Iterator_Callback{
  
public:
    
    sdsl::bit_vector marks;
    Topology_Mapper* mapper;
    
    void init(BIBWT& index, int64_t rev_st_bpr_length, Topology_Mapper& mapper){
        (void) index;
        marks = sdsl::bit_vector(rev_st_bpr_length,0);
        this->mapper = &mapper;
    }
    
    virtual void callback(const Iterator::Stack_frame& top){
        if(top.is_maxrep){
            marks[mapper->leaves_to_node(top.intervals.reverse)] = 1;
        }
    }
    
    virtual void finish(){}
    
    sdsl::bit_vector get_result(){
        return marks;
    }
};

class SLT_Maximal_Marks_Callback : public Iterator_Callback{
    
public:
    
    sdsl::bit_vector marks;
    sdsl::select_support_mcl<1> slt_bpr_ss;
    int64_t preorder_rank;
    bool enabled;
    
    SLT_Maximal_Marks_Callback() : enabled(true) {}
    
    void enable() {enabled = true;}
    void disable() {enabled = false;}
    
    void init(BIBWT& index, sdsl::bit_vector& slt_bpr){
        if(!enabled) return;
        (void) index;
        marks = sdsl::bit_vector(slt_bpr.size(),0);
        sdsl::util::init_support(slt_bpr_ss, &slt_bpr);
        preorder_rank = 0;
    }
    
    virtual void callback(const Iterator::Stack_frame& top){
        if(!enabled) return;
        preorder_rank++;
        if(top.is_maxrep){ // todo: is_maxrep from the stack frame
            marks[slt_bpr_ss.select(preorder_rank)] = 1;
        }
    }
    
    virtual void finish(){}
    
    sdsl::bit_vector get_result(){
        return marks;
    }
};

class Store_Depths_Callback : public Iterator_Callback{
// Needs a left-to-right DFS traversal
public:
    
    vector<int64_t> depths;
    int64_t maxreps_seen;
    int64_t maxdepth;
    bool enabled;
    
    void init(){
        maxreps_seen = 0;
        maxdepth = 0;
    }
    
    void enable(){
        enabled = true;
    }
    void disable(){ 
        enabled = false;
    }
    
    virtual void callback(const Iterator::Stack_frame& top){
        if(!enabled) return;
        if(top.is_maxrep){
            depths.push_back(top.depth);
            maxdepth = max(maxdepth, top.depth);
            maxreps_seen++;
        }
    }
    
    virtual void finish(){
        if(!enabled) return;
    }
    
    sdsl::int_vector<0> get_result(){
        if(!enabled){
            // Return empty vector
            // Set width to 1 in case width 0 is forbidden
            sdsl::int_vector<0> empty(0,0,1);
            return empty;
        } else{
            sdsl::int_vector<0> v(depths.size(), 0, sdsl::bits::hi(maxdepth)+1);
            for(int64_t i = 0; i < depths.size(); i++){
                v[i] = depths[i];
            }
            return v;
        }
    }
};

void iterate_with_callbacks(Iterator& iterator, Iterator_Callback* cb){
    
    iterator.init();
    while(iterator.next()){
        cb->callback(iterator.get_top());
    }
    
    cb->finish();
}

void iterate_with_callbacks(Iterator& iterator, vector<Iterator_Callback*>& callbacks){
    
    iterator.init();
    while(iterator.next()){
        for(Iterator_Callback* cb : callbacks){
            cb->callback(iterator.get_top());
        }
    }
    
    for(Iterator_Callback* cb : callbacks){
        cb->finish();
    }
}

sdsl::bit_vector get_slt_topology(BIBWT& index, Iterator& iterator){
    
    Build_SLT_BPR_Callback callback;
    callback.init(index);
    iterate_with_callbacks(iterator, &callback);
    return callback.get_result();
}

sdsl::bit_vector get_slt_topology(BIBWT& index){
    
    SLT_Iterator iterator(&index);
    return get_slt_topology(index, iterator);
}

Rev_st_topology get_rev_st_bpr_and_pruning(BIBWT& index, Iterator& iterator){
    Build_REV_ST_BPR_And_Pruning_Callback callback;
    callback.init(index);
    iterate_with_callbacks(iterator, &callback);
    return callback.get_result();
    
    //return {bpr_sdsl, pruning_marks};
}

sdsl::bit_vector get_rev_st_topology(BIBWT& index){
    Rev_ST_Iterator iterator(&index);
    return get_rev_st_bpr_and_pruning(index, iterator).bpr;
}

// Takes in an iterator that visits a set of nodes V such that
// (left-extensions of maxreps) \subseteq V \subseteq all nodes
template <typename topology_mapper_t>
sdsl::bit_vector get_rev_st_maximal_marks(BIBWT& index, int64_t rev_st_bpr_length, Iterator& it, topology_mapper_t& mapper){
    Rev_ST_Maximal_Marks_Callback callback;
    callback.init(index,rev_st_bpr_length,mapper);
    iterate_with_callbacks(it, &callback);
    return callback.get_result();
}

template <typename topology_mapper_t>
sdsl::bit_vector get_rev_st_maximal_marks(BIBWT& index, int64_t rev_st_bpr_length, topology_mapper_t& mapper){
    Rev_ST_Iterator iterator(&index);
    return get_rev_st_maximal_marks(index, rev_st_bpr_length, iterator, mapper);
}

sdsl::bit_vector get_slt_maximal_marks(BIBWT& index, sdsl::bit_vector& slt_bpr, Iterator& it){
    SLT_Maximal_Marks_Callback callback;
    callback.init(index, slt_bpr);
    iterate_with_callbacks(it, &callback);
    return callback.get_result();
}

sdsl::bit_vector get_slt_maximal_marks(BIBWT& index, sdsl::bit_vector& slt_bpr){
    SLT_Iterator iterator(&index);
    return get_slt_maximal_marks(index,slt_bpr,iterator);
}


template <typename bwt_t>
sdsl::int_vector<0> get_all_rev_st_maxrep_string_depths(bwt_t& index){
    SLT_Iterator it(&index);
    Store_Depths_Callback cb;
    cb.init();
    iterate_with_callbacks(it, &cb);
    return cb.get_result();
}

#endif
