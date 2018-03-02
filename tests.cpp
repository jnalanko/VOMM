#include "String_Depth_Support_tests.hh"
#include "Maxreps_tester.hh"
#include "MS_enumerator_tests.hh"
#include "LMA_Support_Tests.hh"
#include "Parent_Support_Tests.hh"
#include "context_marking_tests.hh"
#include "score_string_tests.hh"
#include "Topology_tests.hh"

#include <iostream>

using namespace std;

int main(int argc, char** argv){
        
    (void) argc; (void) argv; // Silence unused variable compiler warning

    score_string_tests();
    test_serialization();
    test_RLE();
    test_mark_contexts_entropy_all();
    test_mark_contexts_p_norm_all();
    test_mark_contexts_KL_all();
    test_mark_contexts_formulas_234_all();
    test_brute_bpr_building();
    test_rev_st_bpr_building();
    test_maxrep_rev_st_bpr_building();
    test_maxrep_depth_bounded_rev_st_bpr_building();
    LMA_Support_Tests();
    String_Depth_Support_tests();
    Maxreps_tests();
    MS_Enumerator_tests();
    Parent_Support_Tests();

    cout << "All tests passed" << endl;
    
}



