#ifndef BPR_COLEX_MAPPING
#define BPR_COLEX_MAPPING

#include "Interfaces.hh"


class Full_Topology_Mapper : public Topology_Mapper{
    
public:
    
    typedef int64_t node_t;
    
    sdsl::bp_support_g<>* rev_st_bps;
    sdsl::select_support_mcl<10,2>* rev_st_ss_10;
    sdsl::rank_support_v<10,2>* rev_st_rs_10;
    
    Full_Topology_Mapper() {}
    Full_Topology_Mapper(sdsl::bp_support_g<>* rev_st_bps,
                         sdsl::select_support_mcl<10,2>* rev_st_ss_10,
                         sdsl::rank_support_v<10,2>* rev_st_rs_10):
        rev_st_bps(rev_st_bps),
        rev_st_ss_10(rev_st_ss_10),
        rev_st_rs_10(rev_st_rs_10) {}
    
    node_t leaves_to_node(Interval leaves){
        return enclose_leaves(leaves.left,leaves.right,*rev_st_ss_10,*rev_st_bps).left; // TODO: make self-contained
    }
    
    Interval node_to_leaves(node_t node){ // TODO: make self-contained
        int64_t close = rev_st_bps->find_close(node);
        return bpr_interval_to_leaf_interval(Interval(node, close), *rev_st_rs_10);
    }
    
    node_t find_close(node_t open){
        return rev_st_bps->find_close(open);
    }
    
};

class Pruned_Topology_Mapper : public Topology_Mapper{
// Also works if the topology is not pruned
    
public:
    
    typedef int64_t node_t;
    
    sdsl::bp_support_g<>* rev_st_bps;
    sdsl::select_support_mcl<10,2>* rev_st_ss_10;
    sdsl::rank_support_v<10,2>* rev_st_rs_10;
    sdsl::rank_support_v<1>* pruning_marks_rs;
    sdsl::select_support_mcl<1>* pruning_marks_ss;
    
    Pruned_Topology_Mapper() {}
    Pruned_Topology_Mapper(sdsl::bp_support_g<>* rev_st_bps,
                         sdsl::select_support_mcl<10,2>* rev_st_ss_10,
                         sdsl::rank_support_v<10,2>* rev_st_rs_10,
                         sdsl::rank_support_v<1>* pruning_marks_rs,
                         sdsl::select_support_mcl<1>* pruning_marks_ss):
        rev_st_bps(rev_st_bps),
        rev_st_ss_10(rev_st_ss_10),
        rev_st_rs_10(rev_st_rs_10),
        pruning_marks_rs(pruning_marks_rs),
        pruning_marks_ss(pruning_marks_ss){}
    
    node_t leaves_to_node(Interval leaves){
            
        // Note: if pruning_marks is all ones, rank is the identity function
        int64_t leftmost_1 = pruning_marks_rs->rank(leaves.left+1); 
        int64_t leftmost_2 = pruning_marks_rs->rank(leaves.right+1);
        
        node_t leaf1 = rev_st_ss_10->select(leftmost_1) - 1; // Select gives the closing paren, so -1 to get the open
        node_t leaf2 = rev_st_ss_10->select(leftmost_2) - 1; // Select gives the closing paren, so -1 to get the open
        if(leaf1 == leaf2) return leaf1;
        
        return rev_st_bps->double_enclose(leaf1, leaf2);
    }
    
    Interval node_to_leaves(node_t node){
        // Number of terminal nodes srictly before the node
        int64_t x1 = rev_st_rs_10->rank(node);
        
        // Number of terminal nodes before and inside the node
        int64_t x2 = rev_st_rs_10->rank(find_close(node)+1); // +1: make the rank inclusive
        
        // Note: if pruning marks is all ones, then the select function is f(x) = x-1
        int64_t left = pruning_marks_ss->select(x1+1);
        
        // Note: if pruning marks is all ones, then this is the length of the bwt
        int64_t total_leaves = pruning_marks_rs->rank(pruning_marks_rs->size());
        int64_t right;
        if(x2 == total_leaves){
            // Rightmost leaf is the last leaf
            right = pruning_marks_rs->size()-1; // Last colex position
        } else{
            // Rightmost leaf is not the last leaf. The colex-interval should extend
            // up to but not including the next pruning-marked node
            right = pruning_marks_ss->select(x2+1)-1;
        }
        
        return Interval(left, right);
    }
    
    node_t find_close(node_t open){
        return rev_st_bps->find_close(open);
    }
    
};

#endif