#ifndef PARENT_SUPPORT_HH
#define PARENT_SUPPORT_HH

#include <iostream>
#include "BD_BWT_index/include/BD_BWT_index.hh"
#include "sdsl/rank_support_v.hpp"
#include "sdsl/select_support_mcl.hpp"
#include "sdsl/bp_support_g.hpp"
#include "sdsl/util.hpp"
#include "globals.hh"
#include "Interfaces.hh"

// The purpose of this class is to provide functions to map a lexicographic
// range into the lexicographic range of the parent
class Parent_Support{
    
public:
    
    std::shared_ptr<Bitvector> bpr;
    
    Parent_Support() {};
    Parent_Support(std::shared_ptr<Bitvector> bpr) :
        bpr(bpr) {}
        
    Interval lex_parent(Interval I){ // TODO: delete
        return topology_to_lex(parent(lex_to_topology(I)));
    }

    Interval lex_to_topology(Interval I){ // Takes a lex interval. TODO: DELETE (not the responsibility of this class)
        if(I.left == I.right) return get_leaf_bpr(I.left);

        Interval L = get_leaf_bpr(I.left);
        Interval R = get_leaf_bpr(I.right);
        return LCA(L,R);
    }
    
    Interval parent(Interval I){ // Takes a topology interval
        if(I.left == 0) return I; // Root
        else{
            int64_t open = bpr->enclose(I.left);
            int64_t close = bpr->find_close(open);
            return Interval(open,close);
        }
    } 
    
    Interval topology_to_lex(Interval I){ // Takes a topology interval. TODO: DELETE (not the responsibility of this class)
        int64_t left = bpr->rank_10(I.left);
        int64_t right = bpr->rank_10(I.right + 1) - 1;
        return Interval(left,right);
    }
        
private:
    
    Interval get_leaf_bpr(int64_t leaf_rank){
        int64_t close = bpr->select_10(leaf_rank+1); // Indexing starts from 1, hence the +1
        return Interval(close-1,close);
    }

    // A must be to the left of B and B can not be nested inside A
    Interval LCA(Interval A, Interval B){
        assert(A.right < B.left);
        int64_t open = bpr->double_enclose(A.left, B.left);
        int64_t close = bpr->find_close(open);
        return Interval(open,close);
    }

};


// sdsl::bit_vector bpr = {1, 1, 1,0,1,0, 0, 1, 1,0, 0, 0};
// Parent_Support PS(bpr);
// for(int64_t i = 0; i < bpr.size(); i++){
//     cout << i << " " << PS.rs_10.rank(i) << endl;
// }

// Prints:

// 0 0
// 1 0
// 2 0
// 3 0
// 4 1
// 5 1
// 6 2
// 7 2
// 8 2
// 9 2
// 10 3
// 11 3

#endif
