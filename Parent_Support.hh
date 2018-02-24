#ifndef PARENT_SUPPORT_HH
#define PARENT_SUPPORT_HH

#include <iostream>
#include "BD_BWT_index/include/BD_BWT_index.hh"
#include "sdsl/rank_support_v.hpp"
#include "sdsl/select_support_mcl.hpp"
#include "sdsl/bp_support_g.hpp"
#include "sdsl/util.hpp"
#include "globals.hh"

// The purpose of this class is to provide functions to map a lexicographic
// range into the lexicographic range of the parent
class Parent_Support{
    
public:
    
    sdsl::select_support_mcl<10,2>* ss_10; // find the i-th leaf in the bpr
    sdsl::rank_support_v<10,2>* rs_10; // inverse of find the i-th leaf
    sdsl::bp_support_g<>* bps; // enclose (= parent)
    
    Parent_Support() {};
    Parent_Support(sdsl::select_support_mcl<10,2>* ss_10, sdsl::rank_support_v<10,2>* rs_10, sdsl::bp_support_g<>* bps) :
        ss_10(ss_10), rs_10(rs_10), bps(bps) {}
        
    // Constructor for the lazy and tests
    Parent_Support(sdsl::bit_vector& bpr){
        
        // TODO: These leak memory. Can't just delete in destructor because of the other constructor.
        this->ss_10 = new sdsl::select_support_mcl<10,2>();
        this->rs_10 = new sdsl::rank_support_v<10,2>();
        this->bps = new sdsl::bp_support_g<>();

        sdsl::util::init_support(*this->ss_10, &bpr);
        sdsl::util::init_support(*this->rs_10, &bpr);
        sdsl::util::init_support(*this->bps, &bpr);
    }
    
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
            int64_t open = bps->enclose(I.left);
            int64_t close = bps->find_close(open);
            return Interval(open,close);
        }
    } 
    
    Interval topology_to_lex(Interval I){ // Takes a topology interval. TODO: DELETE (not the responsibility of this class)
        int64_t left = rs_10->rank(I.left);
        int64_t right = rs_10->rank(I.right + 1) - 1;
        return Interval(left,right);
    }
        
private:
    
    Interval get_leaf_bpr(int64_t leaf_rank){
        int64_t close = ss_10->select(leaf_rank+1); // Indexing starts from 1, hence the +1
        return Interval(close-1,close);
    }

    // A must be to the left of B and B can not be nested inside A
    Interval LCA(Interval A, Interval B){
        assert(A.right < B.left);
        int64_t open = bps->double_enclose(A.left, B.left);
        int64_t close = bps->find_close(open);
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
