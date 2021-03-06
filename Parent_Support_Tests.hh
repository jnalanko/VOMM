#ifndef PARENT_SUPPORT_TESTS
#define PARENT_SUPPORT_TESTS

#include <iostream>
#include "Parent_Support.hh"

using namespace std;

void Parent_Support_Tests(){
    sdsl::bit_vector bpr = {1, 1, 1,0,1,0, 0, 1, 1,0,1,0,1,0, 0, 0};
    std::shared_ptr<Basic_bitvector> bprv = make_shared<Basic_bitvector>(bpr);
    bprv->init_bps_support();
    bprv->init_rank_10_support();
    bprv->init_select_10_support();
    Parent_Support PS(bprv);

    assert(PS.lex_parent(Interval(0,0)) == Interval(0,1));
    assert(PS.lex_parent(PS.lex_parent(Interval(0,0))) == Interval(0,4));
    assert(PS.lex_parent(PS.lex_parent(PS.lex_parent(Interval(0,0)))) == Interval(0,4)); // parent of root should be the root
    
    assert(PS.lex_parent(Interval(1,1)) == Interval(0,1));
    assert(PS.lex_parent(Interval(2,2)) == Interval(2,4));
    assert(PS.lex_parent(Interval(3,3)) == Interval(2,4));
    assert(PS.lex_parent(Interval(4,4)) == Interval(2,4));
    
    cerr << "Parent support tests passed" << endl;
}

#endif
