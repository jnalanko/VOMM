#include "String_Depth_Support_tests.hh"
#include "Maxreps_tester.hh"
#include "MS_enumerator_tests.hh"
#include "LMA_Support_Tests.hh"
#include "Parent_Support_Tests.hh"
#include "context_marking_tests.hh"
#include "score_string_tests.hh"
#include "Topology_tests.hh"
#include "logging.hh"

#include <iostream>

using namespace std;

int main(int argc, char** argv){
        
    (void) argc; (void) argv; // Silence unused variable compiler warning

    disable_logging();

    score_string_tests();
    String_Depth_Support_tests();
    test_precomputed_depths();
    test_serialization();
    test_recursive_scoring();
    test_mark_contexts_entropy_all();
    test_mark_contexts_p_norm_all();
    test_mark_contexts_KL_all();
    test_mark_contexts_formulas_234_all();
    test_brute_bpr_building();
    test_rev_st_bpr_building();
    test_maxrep_rev_st_bpr_building();
    test_maxrep_depth_bounded_rev_st_bpr_building();
    Maxreps_tests();
    test_RLE();
    LMA_Support_Tests();
    MS_Enumerator_tests();
    Parent_Support_Tests();

    cout << "All tests passed" << endl;
    
}



