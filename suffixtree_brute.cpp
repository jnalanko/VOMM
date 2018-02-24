#include <iostream>
#include <string>
#include <vector>
#include <set>
#include <map>
#include <sstream>
#include <fstream>

#include "suffixtree_brute.hh"

int main(int argc, char** argv){
    if(argc == 1){
        cerr << "Give file as an argument" << endl;
        return -1;
    }
    ifstream in(argv[1]);
    string S; in >> S; // Assuming no spaces
    if(!in.good()){
        cerr << "Error reading file " << argv[1] << endl;
        return -1;
    }

    string Srev(S.rbegin(), S.rend());
    Srev += '$';
    vector<Edge> edges = get_suffix_tree(Srev);
    set<string> maxreps = get_maxreps(Srev);
    cout << ST_to_dot(edges, maxreps) << endl;
  
}
