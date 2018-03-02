#ifndef BPR_COLEX_MAPPING
#define BPR_COLEX_MAPPING

#include "Interfaces.hh"


class Full_Topology_Mapper : public Topology_Mapper{
    
public:
    
    typedef int64_t node_t;
    
    std::shared_ptr<Bitvector> rev_st_bpr;
    
    /*std::shared_ptr<sdsl::bp_support_g<>> rev_st_bps;
    std::shared_ptr<Select_10_support> rev_st_ss_10;
    std::shared_ptr<Rank_10_support> rev_st_rs_10;*/
    
    Full_Topology_Mapper() {}
    Full_Topology_Mapper(std::shared_ptr<Bitvector> rev_st_bpr) : rev_st_bpr(rev_st_bpr) {}
    /*Full_Topology_Mapper(std::shared_ptr<sdsl::bp_support_g<>> rev_st_bps,
                         std::shared_ptr<Select_10_support> rev_st_ss_10,
                         std::shared_ptr<Rank_10_support> rev_st_rs_10):
        rev_st_bps(rev_st_bps),
        rev_st_ss_10(rev_st_ss_10),
        rev_st_rs_10(rev_st_rs_10) {}*/
        
    Interval LCA(Interval A, Interval B){
        assert(A.right < B.left);
        int64_t open = rev_st_bpr->double_enclose(A.left, B.left);
        int64_t close = rev_st_bpr->find_close(open);
        return Interval(open,close);
    }

    Interval find_leaf_in_bpr(int64_t leaf_rank){
        int64_t close = rev_st_bpr->select_10(leaf_rank+1); // Indexing starts from 1, hence the +1
        return Interval(close-1,close);
    }

    Interval enclose_leaves(int64_t leaf1, int64_t leaf2){
        if(leaf1 == leaf2) return find_leaf_in_bpr(leaf1);
        Interval L = find_leaf_in_bpr(leaf1);
        Interval R = find_leaf_in_bpr(leaf2);
        return LCA(L,R);
    }
    
    Interval bpr_interval_to_leaf_interval(Interval I){
        int64_t leaves_before = rev_st_bpr->rank_10(I.left);
        int64_t leaves_inside = rev_st_bpr->rank_10(I.right+1);
        return Interval(leaves_before, leaves_inside-1); //-1: make end of range inclusive
    }
    
    node_t leaves_to_node(Interval leaves){
        return enclose_leaves(leaves.left,leaves.right).left;
    }
    
    Interval node_to_leaves(node_t node){
        int64_t close = rev_st_bpr->find_close(node);
        return bpr_interval_to_leaf_interval(Interval(node, close));
    }
    
    node_t find_close(node_t open){
        return rev_st_bpr->find_close(open);
    }
    
};

class Pruned_Topology_Mapper : public Topology_Mapper{
// Also works if the topology is not pruned
    
public:
    
    typedef int64_t node_t;
    
    std::shared_ptr<Bitvector> rev_st_bpr;
    std::shared_ptr<Bitvector> pruning_marks;
    
    Pruned_Topology_Mapper() {}
    Pruned_Topology_Mapper(std::shared_ptr<Bitvector> rev_st_bpr, std::shared_ptr<Bitvector> pruning_marks) :
        rev_st_bpr(rev_st_bpr), pruning_marks(pruning_marks) {}

    
    node_t leaves_to_node(Interval leaves){
            
        // Note: if pruning_marks is all ones, rank is the identity function
        int64_t leftmost_1 = pruning_marks->rank(leaves.left+1); 
        int64_t leftmost_2 = pruning_marks->rank(leaves.right+1);
        
        node_t leaf1 = rev_st_bpr->select_10(leftmost_1) - 1; // Select gives the closing paren, so -1 to get the open
        node_t leaf2 = rev_st_bpr->select_10(leftmost_2) - 1; // Select gives the closing paren, so -1 to get the open
        if(leaf1 == leaf2) return leaf1;
        
        return rev_st_bpr->double_enclose(leaf1, leaf2);
    }
    
    Interval node_to_leaves(node_t node){
        // Number of terminal nodes srictly before the node
        int64_t x1 = rev_st_bpr->rank_10(node);
        
        // Number of terminal nodes before and inside the node
        int64_t x2 = rev_st_bpr->rank_10(find_close(node)+1); // +1: make the rank inclusive
        
        // Note: if pruning marks is all ones, then the select function is f(x) = x-1
        int64_t left = pruning_marks->select(x1+1);
        
        // Note: if pruning marks is all ones, then this is the length of the bwt
        int64_t total_leaves = pruning_marks->rank(pruning_marks->size());
        int64_t right;
        if(x2 == total_leaves){
            // Rightmost leaf is the last leaf
            right = pruning_marks->size()-1; // Last colex position
        } else{
            // Rightmost leaf is not the last leaf. The colex-interval should extend
            // up to but not including the next pruning-marked node
            right = pruning_marks->select(x2+1)-1;
        }
        
        return Interval(left, right);
    }
    
    node_t find_close(node_t open){
        return rev_st_bpr->find_close(open);
    }
    
};

#endif