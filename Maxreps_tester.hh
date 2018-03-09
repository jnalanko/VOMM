#ifndef MAXREPS_TESTER_HH
#define MAXREPS_TESTER_HH

#include "brute_tools.hh"
#include "suffixtree_brute.hh"
#include "Precalc.hh"
#include "Maxreps.hh"
#include "BD_BWT_index/include/BD_BWT_index.hh"
#include <iostream>
#include <algorithm>

using namespace std;

/*
void depth_bounded_vs_not(string S){
    BD_BWT_index<> bibwt((uint8_t*)S.c_str());
    Rev_ST_Depth_Bounded_Maxrep_Iterator depth_it(&bibwt, 1e18);
    Rev_ST_Maxrep_Iterator notdepth_it(&bibwt);
    
    vector<pair<Interval, Interval>> depth_intervals;
    depth_it.init();
    while(depth_it.next()){
        depth_intervals.push_back({depth_it.get_top().intervals.reverse, depth_it.get_top().intervals.forward});
    }
    
    vector<pair<Interval, Interval> > nondepth_intervals;
    notdepth_it.init();
    while(notdepth_it.next()){
        nondepth_intervals.push_back({notdepth_it.get_top().intervals.reverse, notdepth_it.get_top().intervals.forward});
    }
    
    sort(depth_intervals.begin(), depth_intervals.end());
    sort(nondepth_intervals.begin(), nondepth_intervals.end());
    
    ofstream depth_out("aaa/depth.txt");
    ofstream nondepth_out("aaa/nondepth.txt");
    for(auto I : depth_intervals) depth_out << I.first.toString() << " " << I.second.toString() << endl;
    for(auto I : nondepth_intervals) nondepth_out << I.first.toString() << " " << I.second.toString() << endl;
} */

void randomized_tests(){
    cerr << "Running randomized maxrep tests" << endl;
    srand(4141414);
    for(int64_t rep = 0; rep < 50; rep++){
        string S = get_random_string(200,3);
        //depth_bounded_vs_not(S);
        set<string> maxreps;
        
        BD_BWT_index<> bibwt((uint8_t*)S.c_str());
        Rev_ST_Depth_Bounded_Maxrep_Iterator iter(1e18);
        iter.set_index(&bibwt);
        
        sdsl::bit_vector bpr = get_rev_st_bpr_and_pruning(bibwt, iter).bpr;
        //sdsl::bp_support_g<> bps(bpr);
        
        sdsl::bit_vector pruning = get_rev_st_bpr_and_pruning(bibwt, iter).pruning_marks;
        
        Pruned_Topology_Mapper mapper(make_shared<Basic_bitvector>(bpr), make_shared<Basic_bitvector>(pruning));
        mapper.rev_st_bpr->init_rank_10_support();
        mapper.rev_st_bpr->init_select_10_support();
        mapper.rev_st_bpr->init_bps_support();
        mapper.pruning_marks->init_rank_support();
        mapper.pruning_marks->init_select_support();
        
        sdsl::bit_vector marks = get_rev_st_maximal_marks(bibwt, bpr.size(), iter, mapper);
        
        for(int64_t i = 0; i < marks.size(); i++){
            if(bpr[i] == 1 && marks[i] == 1){
                Interval colex = mapper.node_to_leaves(i);
                string X = colex_range_to_string(S, colex.left, colex.right);
                maxreps.insert(X);
            }
        }
        
        set<string> maxreps_correct = get_maxreps(S); // The correct answer
        /*cout << maxreps.size() << endl;
        ofstream right("aaa/right.txt");
        ofstream wrong("aaa/wrong.txt");
        for(string S : maxreps_correct) right << S << endl;
        for(string S : maxreps) wrong << S << endl;*/
        
        assert(maxreps == maxreps_correct);
    }
}

void Maxreps_tests(){
    randomized_tests();
    
    string s = "mississippi";
    BD_BWT_index<> index((uint8_t*)s.c_str());
    auto maxreps = find_maxreps(index);
    sort(maxreps.begin(), maxreps.end());
    assert(maxreps.size() == 5);
    vector<pair<string, Interval_pair> > answers =
      {make_pair(string(""), Interval_pair(0,11,0,11)),
       make_pair(string("i"), Interval_pair(1,4,1,4)),
       make_pair(string("issi"), Interval_pair(3,4,3,4)),
       make_pair(string("p"), Interval_pair(6,7,6,7)),
       make_pair(string("s"), Interval_pair(8,11,8,11))
    };
    assert(maxreps == answers);
}

#endif
