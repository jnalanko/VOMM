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
#include "logging.hh"
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

class Topology_Algorithms : public Topology{

public:

    Global_Data* data; // Todo: maybe don't need to store this pointer
    Topology_Mapper* mapper;
    Parent_Support PS;
    String_Depth_Support* SDS;
    LMA_Support LMAS;

    typedef int64_t node_t;
    
    Topology_Algorithms(Global_Data* data, Topology_Mapper* mapper, String_Depth_Support* SDS, Parent_Support PS, LMA_Support LMAS) :
        data(data),
        mapper(mapper),
        PS(PS),
        SDS(SDS),
        LMAS(LMAS)
    {}

    // Operations in topology:
    int64_t rev_st_string_depth(node_t node){
        return SDS->string_depth(node);
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

template<typename inputstream_t>
double main_loop(inputstream_t& S, Global_Data& data, Topology& topo_alg, Scoring_Function& scorer, Loop_Invariant_Updater& updater){
    double logprob = 0;
    Interval I(0,data.revbwt->size()-1); // todo: or: index.empty_string()
    int64_t string_depth = 0;
    char c;
    
    while(S.getchar(c)){
        
        // Compute log-probability of c
        logprob += scorer.score(I, string_depth, c, topo_alg, *data.revbwt, data);
        
        // Update I and string_depth
        pair<Interval, int64_t> new_values = updater.update(I, string_depth, c, data, topo_alg, *data.revbwt);
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
            if(I.size() == index.size()) return {I,0};
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
            if(I.size() == index.size()) return {I,0};
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
            if(I.size() == index.size()) return {I,0};
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
    bool maxrep_contexts;

    Basic_Scorer(double escape_prob, bool maxrep_contexts) : escape_prob(escape_prob), maxrep_contexts(maxrep_contexts) {}

    virtual double score(Interval I, int64_t d, char c, Topology& topology, BWT& index, Global_Data& G){
        int64_t node = topology.leaves_to_node(I);
        if(maxrep_contexts && G.rev_st_maximal_marks->at(node) == 1 && topology.rev_st_string_depth(node) > d){
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
            return log2(escape_prob);
        } else{
            // don't count in the dollar the context in the interval of the empty string, hence -1
            return log2((double)R.size()) - log2(min(I.size(), index.size()-1));
        }
    }

};

class Recursive_Scorer : public Scoring_Function {

public:

    double escape_prob;
    bool maxrep_contexts;
    
    int64_t get_context_depth(int64_t open, Topology& topology){
        if(maxrep_contexts) return topology.rev_st_string_depth(open);
        else{
            if(open == 0) return 0; // Root is a special case: It's not a left-extension of a maxrep
            return topology.rev_st_string_depth(topology.rev_st_parent(open)) + 1; // One-char left-extension of a maxrep
        }        
    }

    Recursive_Scorer(double escape_prob, bool maxrep_contexts) 
    : escape_prob(escape_prob), maxrep_contexts(maxrep_contexts) {}

    virtual double score(Interval I, int64_t d, char c, Topology& topology, BWT& index, Global_Data& G){
        int64_t node = topology.leaves_to_node(I);
        if(maxrep_contexts && G.rev_st_maximal_marks->at(node) == 1 && topology.rev_st_string_depth(node) > d){
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
            // Did not find c
            if(d == 0) return log2(escape_prob); // Root
            
            int64_t depth1 = get_context_depth(node, topology);
            node = topology.rev_st_lma(topology.rev_st_parent(node));
            I = topology.node_to_leaves(node);
            int64_t depth2 = get_context_depth(node, topology);
            assert(depth1 != depth2);
            int64_t distance_travelled = depth1 - depth2;
            double ancestor_score = score(I, depth2, c, topology, index, G);
            for(int64_t i = 0; i < distance_travelled; i++) ancestor_score += log2(escape_prob);
            return ancestor_score;
        } else{
            // Found c
            return log2(1 - escape_prob) + log2(R.size()) - log2(min(I.size(), index.size()-1));
            // Don't count in the dollar the context in the interval of the empty string, hence -1
        }
    }

};

class Input_Stream{
// getchar function for a string
public:

    std::string S;
    int64_t pos;
    
    Input_Stream(string& S) : S(S), pos(0) {}
    
    bool getchar(char& c){
        if(pos == S.size()){
            return false;
        } else{
            c = S[pos];
            pos++;
            return true;
        }
    }
};

/**  
 * Scores a string S using the simple method described in the paper:
 *
 * "Probabilistic suffix array: efficient modeling and prediction of protein families"
 * by Jie Lin, Donald Adjeroh and Bing-Hua Jiang.
 *
 * @return the base-2 logarithm of the total probability of S.
 */
template <typename input_stream_t> double score_string_lin(input_stream_t& S, Global_Data& G) {
    const int64_t BWT_SIZE = G.revbwt->size();
    const Interval LARGEST_INTERVAL(0,G.revbwt->size()-1);
    const Interval DUMMY_INTERVAL(-1,-1);
    
    char c;
    int64_t sizeFrom, sizeTo, node;
    double logSizeFrom, logSizeTo, out;
    Pruned_Topology_Mapper mapper(G.rev_st_bpr,G.pruning_marks); // Also works for non-pruned topology
    Parent_Support PS(G.rev_st_bpr);
    Interval I_W, I_Wc;
    
    out=0.0;
    I_W=LARGEST_INTERVAL;
    sizeFrom=BWT_SIZE;
    logSizeFrom=log2(sizeFrom-1);  // We don't want to count in the final dollar
    while (S.getchar(c)) {
        // Finding the BWT interval of the longest suffix of W that is followed by c
        node=-1;
        while (true) {
            I_Wc=G.revbwt->search(I_W,c);
            sizeTo=I_Wc.size();
            if (sizeTo==0) {
                if (sizeFrom==BWT_SIZE) {  // c does not occur in the text
                    I_Wc=DUMMY_INTERVAL;
                    break;
                }
                if (node==-1) node=mapper.leaves_to_node(I_W);  // Map to topology
                node=PS.parent(node);  // Take parent
                I_W=mapper.node_to_leaves(node);  // Map back to revbwt
                sizeFrom = I_W.size();
            } else break;
        }
        
        if (I_Wc==DUMMY_INTERVAL) continue;
        if(node!=1) logSizeFrom = log2(min(BWT_SIZE-1,sizeFrom));
        
        // Cumulating the probability
        logSizeTo=log2(sizeTo);
        out+=logSizeTo-logSizeFrom;
        
        // Next iteration
        logSizeFrom=logSizeTo;
        I_W=I_Wc;
        sizeFrom=sizeTo;
    }
    return out;
}

// Input stream must have a function getchar(char& c), which returns
// false it the end of the stream was reached
template <typename input_stream_t>
double score_string(input_stream_t& S, Global_Data& G, Scoring_Function& scorer, Loop_Invariant_Updater& updater){
    
    assert(G.revbwt != nullptr);
    Pruned_Topology_Mapper mapper; // Also works for non-pruned topology
    init_support(mapper, &G);
    
    std::shared_ptr<String_Depth_Support> SDS;
    if(G.have_slt()){
        SDS = make_shared<String_Depth_Support_SLT>(G.rev_st_bpr,G.slt_bpr,G.rev_st_maximal_marks,G.slt_maximal_marks);
    } else{
        SDS = make_shared<String_Depth_Support_Store_All>(G.string_depths, G.rev_st_maximal_marks);
    }
    
    Parent_Support PS;
    LMA_Support LMAS;
    
    init_support<Parent_Support>(PS, &G);
    init_support<LMA_Support>(LMAS, &G);
    
    Topology_Algorithms topology(&G, &mapper, SDS.get(), PS, LMAS);
    
    return main_loop(S,G,topology,scorer,updater);
}

template <typename index_t = BD_BWT_index<>,
          typename String_Depth_Support_t = String_Depth_Support,
          typename Parent_Support_t = Parent_Support,
          typename LMA_Support_t = LMA_Support>
double score_string(string& S, Global_Data& G, Scoring_Function& scorer, Loop_Invariant_Updater& updater){

    Input_Stream is(S);
    return score_string<Input_Stream>(is,G,scorer,updater);
}

int score_string_main(int argc, char** argv);
    
#endif
