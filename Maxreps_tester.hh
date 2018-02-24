#ifndef MAXREPS_TESTER_HH
#define MAXREPS_TESTER_HH

#include "Maxreps.hh"
#include "BD_BWT_index/include/BD_BWT_index.hh"
#include <iostream>
#include <algorithm>

using namespace std;

void Maxreps_tests(){
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
    cout << "Maxreps test for mississippi OK" << endl;
}

#endif
