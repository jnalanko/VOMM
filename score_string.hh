#ifndef SCORE_STRING_HH
#define SCORE_STRING_HH

#include "BPR_tools.hh"
#include "Parent_Support.hh"
#include "LMA_Support.hh"
#include "String_Depth_Support.hh"
#include "globals.hh"
#include "BWT_iteration.hh"
#include "Precalc.hh"
#include "context_marking.hh"
#include "Interfaces.hh"
#include "BPR_Colex_mapping.hh"
#include "UniBWT.h"
#include "RLEBWT.hh"
#include "Basic_bitvector.hh"
#include "RLE_bitvector.hh"
#include <stack>
#include <vector>
#include <memory>

class Full_Topology_Mapper;

//enum Pruning_Type {NONE,MAXREPS,LENGTH};
//enum Context_Type {ENTROPY,EQ234,PNORM,KL}; // Should be made into a class so it can take arbitrary parameters for thresholds

sdsl::bit_vector get_rev_st_bpr_context_only(Global_Data* G){ // Todo: move to precalc.hh?
    int64_t nMarked = G->rev_st_context_marks->rank(G->rev_st_context_marks->size());
    
    // Build bpr for marked only

    sdsl::bit_vector rev_st_bpr_context_only(nMarked);
    int64_t k = 0;
    for(int64_t i = 0; i < G->rev_st_bpr->size(); i++){
        if(G->rev_st_context_marks->at(i)){
            rev_st_bpr_context_only[k] = G->rev_st_bpr->at(i);
            k++;
        }
    }
    return rev_st_bpr_context_only;
}

template<typename String_Depth_Support_t, typename Parent_Support_t, typename LMA_Support_t>
class Topology_Algorithms : public Topology{

public:

    Global_Data* data; // Todo: maybe don't need to store this pointer
    Topology_Mapper* mapper;
    Parent_Support_t PS;
    String_Depth_Support_t SDS;
    LMA_Support_t LMAS;

    typedef int64_t node_t;
    
    Topology_Algorithms(Global_Data* data, Topology_Mapper* mapper, String_Depth_Support_t SDS, Parent_Support_t PS, LMA_Support_t LMAS) :
        data(data),
        mapper(mapper),
        PS(PS),
        SDS(SDS),
        LMAS(LMAS)
    {}

    // Operations in topology:
    int64_t rev_st_string_depth(node_t node){
        return SDS.string_depth(node);
    }
    
    node_t rev_st_parent(node_t node){
        int64_t close = data->rev_st_bpr->find_close(node);
        return PS.parent(Interval(node,close)).left;
    }
    
    node_t rev_st_lma(node_t node){
        return LMAS.LMA(node);
    }
    
    // Mapping between colex intervals and topology nodes
    node_t leaves_to_node(Interval I){
        return mapper->leaves_to_node(I);
    };
    
    Interval node_to_leaves(node_t node){
        return mapper->node_to_leaves(node);
    }
};

double main_loop(string& S, Global_Data& data, Topology& topo_alg, Scoring_Function& scorer, Loop_Invariant_Updater& updater){
    double logprob = 0;
    Interval I(0,data.revbwt->size()-1); // todo: or: index.empty_string()
    int64_t string_depth = 0;
    for(int64_t i = 0; i < S.size(); i++){
        // Compute probability of S[i]
        logprob += log2(scorer.score(I, string_depth, S[i], topo_alg, *data.revbwt));
        
        // Update I and string_depth
        pair<Interval, int64_t> new_values = updater.update(I, string_depth, S[i], data, topo_alg, *data.revbwt);
        I = new_values.first;
        string_depth = new_values.second;
    }
    return logprob;
}

template<typename T> void init_support(T&, Global_Data*);

template<> void init_support<LMA_Support>(LMA_Support& LMAS, Global_Data* G){
    LMAS = LMA_Support(G->rev_st_context_marks,
             G->rev_st_bpr_context_only);

}

template<> void init_support<Parent_Support>(Parent_Support& PS, Global_Data* G){
    PS = Parent_Support(G->rev_st_bpr);
}

template<> void init_support<String_Depth_Support>(String_Depth_Support& SDS, Global_Data* G){
    SDS = String_Depth_Support(G->rev_st_bpr,G->slt_bpr,G->rev_st_maximal_marks,G->slt_maximal_marks);
}

template<> void init_support<String_Depth_Support_Store_All>(String_Depth_Support_Store_All& SDS, Global_Data* G){
    SDS = String_Depth_Support_Store_All(&G->string_depths);
}

template<> void init_support<Full_Topology_Mapper>(Full_Topology_Mapper& mapper, Global_Data* G){
    mapper = Full_Topology_Mapper(G->rev_st_bpr);
}

template<> void init_support<Pruned_Topology_Mapper>(Pruned_Topology_Mapper& mapper, Global_Data* G){
    mapper = Pruned_Topology_Mapper(G->rev_st_bpr, G->pruning_marks);
}

class Maxrep_Pruned_Updater : public Loop_Invariant_Updater {

public:

    virtual pair<Interval, int64_t> update(Interval I, int64_t d, char c, Global_Data& data, Topology& topology, BWT& index){
        bool parent_taken = false;
        bool first_iteration = true;
        
        // Take parents until right-extension succeeds
        while(index.search(I,c).size() == 0){
            if(first_iteration){
                // Go up the nearest non-pruned node (maxrep)
                I = topology.node_to_leaves(topology.leaves_to_node(I));
                first_iteration = false;
                parent_taken = true;
                continue;
            }
            
            int64_t node = topology.leaves_to_node(I); // Map to topology
            node = topology.rev_st_parent(node); // Take parent
            I = topology.node_to_leaves(node); // Map back to colex interval
            parent_taken = true;
            if(I.size() == 0) break;
        }
        
        if(parent_taken){
            // Need to recalculate depth
            int64_t bpr_pos = topology.leaves_to_node(I);
            if(data.rev_st_maximal_marks->at(bpr_pos)){
                d = topology.rev_st_string_depth(bpr_pos);
            } else{
                int64_t parent = topology.rev_st_parent(bpr_pos);
                assert(data.rev_st_maximal_marks->at(parent)); // Should be at a maxrep
                d = topology.rev_st_string_depth(parent) + 1; // One character left-extension of a maxrep
            }
            
        }
                
        return {index.search(I,c), d+1};
    }
};

class Basic_Updater : public Loop_Invariant_Updater {

public:

    virtual pair<Interval, int64_t> update(Interval I, int64_t d, char c, Global_Data& data, Topology& topology, BWT& index){
        (void) data; // Not needed. Make the compiler happy.
        bool parent_taken = false;
        while(index.search(I,c).size() == 0){
            int64_t node = topology.leaves_to_node(I); // Map to topology
            node = topology.rev_st_parent(node); // Take parent
            I = topology.node_to_leaves(node); // Map back to colex interval
            parent_taken = true;
            if(I.size() == 0) break;
        }
        
        if(parent_taken){
            // Need to recalculate depth
            assert(data.rev_st_maximal_marks->at(topology.leaves_to_node(I))); // Should be at a maxrep
            d = topology.rev_st_string_depth(topology.leaves_to_node(I)); // Guaranteed to be at a maxrep
        }
        
        return {index.search(I,c), d+1};
    }
};



class Maxrep_Depth_Bounded_Updater : public Loop_Invariant_Updater { // The same as non-depth bounded?

public:
    
    virtual pair<Interval, int64_t> update(Interval I, int64_t d, char c, Global_Data& data, Topology& topology, BWT& index){
        bool parent_taken = false;
        bool first_iteration = true;
        while(index.search(I,c).size() == 0){
            if(first_iteration){
                // Go up the the nearest non-pruned node
                I = topology.node_to_leaves(topology.leaves_to_node(I));
                first_iteration = false;
                parent_taken = true; // Not necessarily true
                continue;
            }
            
            int64_t node = topology.leaves_to_node(I); // Map to topology
            node = topology.rev_st_parent(node); // Take parent
            I = topology.node_to_leaves(node); // Map back to colex interval
            parent_taken = true;
            if(I.size() == 0) break;
        }
                
        if(parent_taken){
            // Need to recalculate depth
            int64_t bpr_pos = topology.leaves_to_node(I);
            if(data.rev_st_maximal_marks->at(bpr_pos)){
                d = topology.rev_st_string_depth(bpr_pos);
            } else{
                int64_t parent = topology.rev_st_parent(bpr_pos);
                assert(data.rev_st_maximal_marks->at(parent)); // Should be at a maxrep
                d = topology.rev_st_string_depth(parent) + 1; // One character left-extension of a maxrep
            }
            
        }
        
        return {index.search(I,c), d+1};
    }
};


class Basic_Scorer : public Scoring_Function {

public:

    double escape_prob;
    bool check_depth;

    Basic_Scorer(double escape_prob, bool check_depth) : escape_prob(escape_prob), check_depth(check_depth) {}

    virtual double score(Interval I, int64_t d, char c, Topology& topology, BWT& index){
        int64_t node = topology.leaves_to_node(I);
        if(check_depth && topology.rev_st_string_depth(node) > d){
            // Inside an edge -> lex interval I represents the node at the end that is further
            // away from the root -> need to go to the edge closest to the root first
            // before taking the lowest marked ancestor, because otherwise the lowest common
            // ancestor might give us the node that is further away from the root in case it is marked.
            // This thing is not necessary to do if contexts are left extensions of maxreps
            node = topology.rev_st_parent(node);
        }
        node = topology.rev_st_lma(node);
        I = topology.node_to_leaves(node);

        // Compute the probability of S[i]
        Interval R = index.search(I, c);
        
        if(R.size() == 0){
           return escape_prob;
        } else{
            // don't count in the dollar the context in the interval of the empty string, hence -1
            return (double)R.size() / min(I.size(), index.size()-1);
        }
    }

};

// All components of the model will be stored into G
void build_model(Global_Data& G, string& T, Context_Formula& context_formula,
                        Iterator& slt_it, Iterator& rev_slt_it, Iterator& context_candidate_iterator,
                        bool run_length_coding, bool compute_string_depths){
        
    
    BD_BWT_index<> bibwt((uint8_t*)T.c_str());
    
    slt_it.set_index(&bibwt);
    rev_slt_it.set_index(&bibwt);
    
    Rev_st_topology RSTT = get_rev_st_bpr_and_pruning(bibwt, rev_slt_it);
    G.rev_st_bpr = std::shared_ptr<Bitvector>(new Basic_bitvector(RSTT.bpr));
    
    G.rev_st_bpr->init_rank_10_support();
    G.rev_st_bpr->init_select_10_support();
    G.rev_st_bpr->init_bps_support();
    
    if(run_length_coding){
        G.pruning_marks = std::shared_ptr<Bitvector>(new RLE_bitvector(RSTT.pruning_marks));
    } else{
        G.pruning_marks = std::shared_ptr<Bitvector>(new Basic_bitvector(RSTT.pruning_marks));
    }
    
    G.pruning_marks->init_rank_support();
    G.pruning_marks->init_select_support();
    
    Pruned_Topology_Mapper mapper(G.rev_st_bpr, G.pruning_marks);
            
    G.rev_st_maximal_marks = std::shared_ptr<Bitvector>(new Basic_bitvector(get_rev_st_maximal_marks(bibwt, G.rev_st_bpr->size(), rev_slt_it, mapper)));
    G.rev_st_maximal_marks->init_rank_support();
    
    G.rev_st_context_marks = std::shared_ptr<Bitvector>(new Basic_bitvector(context_formula.get_rev_st_context_marks(&bibwt, G.rev_st_bpr->size(), context_candidate_iterator, mapper)));
    G.rev_st_context_marks->init_rank_support();
    G.rev_st_context_marks->init_select_support();
    
    G.rev_st_bpr_context_only = std::shared_ptr<Bitvector>(new Basic_bitvector(get_rev_st_bpr_context_only(&G)));
    G.rev_st_bpr_context_only->init_bps_support();
    
    sdsl::bit_vector sdsl_slt_bpr = get_slt_topology(bibwt, slt_it);
    if(run_length_coding){
        G.slt_bpr = std::shared_ptr<Bitvector>(new RLE_bitvector(sdsl_slt_bpr));
    } else{
        G.slt_bpr = std::shared_ptr<Bitvector>(new Basic_bitvector(sdsl_slt_bpr));
    }
    G.slt_bpr->init_rank_support();
    
    if(run_length_coding){
        G.slt_maximal_marks = std::shared_ptr<Bitvector>(new RLE_bitvector(get_slt_maximal_marks(bibwt, sdsl_slt_bpr, slt_it)));
    } else{
        G.slt_maximal_marks = std::shared_ptr<Bitvector>(new Basic_bitvector(get_slt_maximal_marks(bibwt, sdsl_slt_bpr, slt_it)));
    }
    G.slt_maximal_marks->init_rank_support();
    G.slt_maximal_marks->init_select_support();
    
    if(compute_string_depths){
        cerr << "Not implemented error: compute string depths" << endl;
        throw(std::runtime_error("Not implemented error: compute string depths"));
        //G.string_depths = CA.get_rev_st_string_depths(G.bibwt, G.rev_st_bpr.size(), mapper);
    }
    
    // Store reverse BWT for scoring. Todo: reuse already computed bibwt
    string T_reverse(T.rbegin(), T.rend());
    if(run_length_coding){
        BWT* bwt = new RLEBWT<>((uint8_t*)T_reverse.c_str());
        G.revbwt = std::move(std::unique_ptr<BWT>(bwt)); // transfer ownership to G
    } else{
        BWT* bwt = new Basic_BWT<>((uint8_t*)T_reverse.c_str());
        G.revbwt = std::move(std::unique_ptr<BWT>(bwt)); // transfer ownership to G
    }
    
    
}

template <typename index_t = BD_BWT_index<>,
          typename String_Depth_Support_t = String_Depth_Support,
          typename Parent_Support_t = Parent_Support,
          typename LMA_Support_t = LMA_Support>
double score_string(string& S, Global_Data& G, Scoring_Function& scorer, Loop_Invariant_Updater& updater){

    assert(G.revbwt != nullptr);
    Pruned_Topology_Mapper mapper; // Also works for non-pruned topology
    init_support(mapper, &G);
    
    String_Depth_Support_t SDS;
    Parent_Support_t PS;
    LMA_Support_t LMAS;
    
    init_support<String_Depth_Support_t>(SDS, &G);
    init_support<Parent_Support_t>(PS, &G);
    init_support<LMA_Support_t>(LMAS, &G);
    
    //typedef Topology_Algorithms<String_Depth_Support,Parent_Support,LMA_Support> topo_alg_t; 
    Topology_Algorithms<String_Depth_Support,Parent_Support,LMA_Support> topology(&G, &mapper, SDS, PS, LMAS);
    
    //Basic_Scorer<topo_alg_t, BD_BWT_index<> > scorer(escape, !aW_contexts);
    //Maxrep_Pruned_Updater<topo_alg_t, BD_BWT_index<> > updater;
    
    return main_loop(S,G,topology,scorer,updater);
}

// double main_loop(string& S, bwt_t& index, topo_alg_t& topo_alg, scorer_t scorer, updater_t updater){
    
    
#endif
