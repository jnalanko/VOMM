    #ifndef TOPOLOGY_HH
#define TOPOLOGY_HH

class Standard_Topology{
public:
    
    typedef int64_t node_t; // Open parenthesis
    
    Standard_Topology(sdsl::bit_vector tree_bpr) : bpr(tree_bpr){
        sdsl::util::init_support(bpr_rs_10, &bpr);
        sdsl::util::init_support(bpr_ss_10, &bpr);
        sdsl::util::init_support(bpr_bps, &bpr);
    }
    
    sdsl::bit_vector bpr;
    sdsl::rank_support_v<10,2> bpr_rs_10;
    sdsl::select_support_mcl<10,2> bpr_ss_10;
    sdsl::bp_support_g<> bpr_bps;
    
    Interval LCA(Interval A, Interval B){
        assert(A.right < B.left);
        int64_t open = bpr_bps.double_enclose(A.left, B.left);
        int64_t close = bpr_bps.find_close(open);
        return Interval(open,close);
    }
    
    Interval find_leaf_in_bpr(int64_t leaf_rank){
        int64_t close = bpr_ss_10.select(leaf_rank+1); // Indexing starts from 1, hence the +1
        return Interval(close-1,close);
    }
    
    node_t leaves_to_node(Interval I){
        int64_t leaf1 = I.left;
        int64_t leaf2 = I.right;
        if(leaf1 == leaf2) return find_leaf_in_bpr(leaf1);
        Interval L = find_leaf_in_bpr(leaf1);
        Interval R = find_leaf_in_bpr(leaf2);
        return LCA(L,R,bps);
    }
    
    Interval node_to_leaves(node_t node){
        int64_t open = node;
        int64_t close = bpr_bps.find_close(open);
        int64_t leaves_before = bpr_rs_10.rank(open);
        int64_t leaves_inside = bpr_rs_10.rank(close+1);
        return Interval(leaves_before, leaves_inside-1); //-1: make end of range inclusive
    }
};

class Pruned_Topology{
public:
    
    typedef int64_t node_t; // Open parenthesis
    
    // BPR representation of the part of the suffix tree that
    // contains only maximal nodes
    sdsl::bit_vector bpr;
    
    // lex_marks[i] == 1 iff leaf i from left to right is the
    // leftmost child of a node after pruning (explain better)
    sdsl::bit_vector lex_marks;
    
    Interval leaves_to_bpr(Interval I){
        // 
    }
    
    Interval bpr_to_leaves(Interval I){
        //
    }
};

#endif
