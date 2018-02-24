#ifndef TOPOLOGY_TESTS_HH
#define TOPOLOGY_TESTS_HH

#include "suffixtree_brute.hh"
#include "brute_tools.hh"
#include "Precalc.hh"
#include <string>
#include <sstream>

using namespace std;

template <typename T>
string vec_to_string(vector<T>& v){
    stringstream ss;
    for(int64_t i = 0; i < v.size(); i++){
        ss << v[i];
        if(i != v.size() - 1)
            ss << ",";
    }
    return ss.str();
}

sdsl::bit_vector to_sdsl(const vector<bool>& v){
    sdsl::bit_vector u(v.size(),0);
    for(int64_t i = 0; i < v.size(); i++)
        u[i] = v[i];
    return u;
}

void test_brute_bpr_building(){

    string text = "abracabra";
    std::reverse(text.begin(), text.end());
    BD_BWT_index<> index((uint8_t*)text.c_str());
    
    vector<Edge> ST = get_suffix_tree(text + '$');
    vector<bool> bpr = ST_to_bpr(ST);
    vector<bool> rev_st_topology_correct = {1,1,0,1,1,0,1,0,1,1,0,1,0,0,0,1,1,0,1,0,0,1,0,1,1,0,1,0,0,0}; // hand checked
    assert(bpr == rev_st_topology_correct);

}

void test_rev_st_bpr_building(){
    srand(525720);
    for(int64_t reps = 0; reps < 100; reps++){
        string S = get_random_string(100, 3);
        BD_BWT_index<> index((uint8_t*)S.c_str());
        sdsl::bit_vector rev_st_bpr = get_rev_st_topology(index);
        
        string S_rev(S.rbegin(), S.rend());
        vector<Edge> ST = get_suffix_tree(S_rev + '$');
        sdsl::bit_vector bpr_brute = to_sdsl(ST_to_bpr(ST));
        
        assert(rev_st_bpr == bpr_brute);
    }
    cout << "Rev ST full topology building OK" << endl;
}

void test_maxrep_rev_st_bpr_building(){
    // Testing maxreps + left extensions of maxreps
    srand(552340);
    for(int64_t reps = 0; reps < 100; reps++){
        string S = get_random_string(100, 3);

        BD_BWT_index<> index((uint8_t*)S.c_str());
        Rev_ST_Maxrep_Iterator iter(&index);
        sdsl::bit_vector rev_st_bpr = get_rev_st_bpr_and_pruning(index,iter).bpr;
        
        string S_rev(S.rbegin(), S.rend());
        vector<Edge> ST = get_suffix_tree(S_rev + '$');
        set<string> maxreps = get_maxreps(S_rev);
        vector<Edge> ST_maxreps;
        for(Edge E : ST){
            if(maxreps.find(E.from) != maxreps.end()){
                // Origin of edge is a maxrep
                ST_maxreps.push_back(E);
            }
        }
        
        sdsl::bit_vector bpr_brute = to_sdsl(ST_to_bpr(ST_maxreps));
        
        assert(rev_st_bpr == bpr_brute);
    }
    cout << "Rev ST maxrep left extension topology building OK" << endl;
}

void test_maxrep_depth_bounded_rev_st_bpr_building(){
    // Testing maxreps + left extensions of maxreps
    srand(12131412);
    for(int64_t reps = 0; reps < 100; reps++){
        string S = get_random_string(100, 3);
        int64_t depth_bound = 3;

        BD_BWT_index<> index((uint8_t*)S.c_str());
        Rev_ST_Depth_Bounded_Maxrep_Iterator iter(&index, depth_bound);
        sdsl::bit_vector rev_st_bpr = get_rev_st_bpr_and_pruning(index,iter).bpr;
        
        string S_rev(S.rbegin(), S.rend());
        vector<Edge> ST = get_suffix_tree(S_rev + '$');
        set<string> maxreps = get_maxreps(S_rev);
        vector<Edge> ST_maxreps;
        for(Edge E : ST){
            if(maxreps.find(E.from) != maxreps.end() && E.from.size() < depth_bound){
                // Origin of edge is a maxrep
                ST_maxreps.push_back(E);
            }
        }
        
        sdsl::bit_vector bpr_brute = to_sdsl(ST_to_bpr(ST_maxreps));        
        assert(rev_st_bpr == bpr_brute);
    }
    cout << "Rev ST depth bounded maxrep left extension topology building OK" << endl;
}

#endif
