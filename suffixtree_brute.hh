#ifndef SUFFIXTREE_BRUTE
#define SUFFIXTREE_BRUTE

#include <iostream>
#include <string>
#include <vector>
#include <set>
#include <map>
#include <sstream>
#include <fstream>

using namespace std;

// There is an edge from string X to string Y iff 
// - X and Y are right-maximal
// - X is a prefix of Y
// - There does not exist a right-maximal prefix Y of Z that is longer than X

struct Edge{
    string from;
    string to;
};

bool is_prefix_of(string X, string Y){
    return Y.substr(0, X.size()) == X;
}

set<string> get_maxreps(string S){
    // Put in sentinels to both ends
    if(S.back() != '$') S += '$';
    if(S[0] != '#') S = '#' + S;
    
    set<char> alphabet;
    for(char c : S) alphabet.insert(c);

    set<string> substrings;
    for(int64_t i = 0; i < S.size(); i++){
        for(int64_t k = 0; k < S.size(); k++){
            substrings.insert(S.substr(i,k));
        }
    }

    set<string> maxreps;
    
    // Add right maximal nodes
    for(string X : substrings){
        int64_t right = 0;
        int64_t left = 0;
        for(char c : alphabet){
            if(substrings.find(X + c) != substrings.end()){
                right++;
            }
            if(substrings.find(c + X) != substrings.end()){
                left++;
            }
        }
        if(right >= 2 && left >= 2) maxreps.insert(X);
    }
    return maxreps;
}

vector<Edge> get_suffix_tree(string S){
    set<char> alphabet;
    for(char c : S) alphabet.insert(c);

    set<string> substrings;
    for(int64_t i = 0; i < S.size(); i++){
        for(int64_t k = 0; k < S.size(); k++){
            substrings.insert(S.substr(i,k));
        }
    }

    set<string> nodes;

    // Add right maximal nodes
    for(string X : substrings){
        int64_t right = 0;
        for(char c : alphabet){
            if(substrings.find(X + c) != substrings.end()){
                right++;
            }
        }
        if(right >= 2) nodes.insert(X);
    }

    // Add leaves
    for(int64_t i = 0; i < S.size(); i++){
        nodes.insert(S.substr(i));
    }
    
    vector<Edge> edges;
    for(string X : nodes){
        for(string Y : nodes){
            if(X == Y) continue;
            if(is_prefix_of(X,Y)){
                bool good = true;
                for(string Z : nodes){
                    if(Z != Y && Z.size() > X.size() && is_prefix_of(Z,Y)) 
                        good = false;
                }
                if(good){
                    edges.push_back({X,Y});
                }
            }
        }
    }
    return edges;
}

vector<bool> ST_to_bpr(vector<Edge> edges){
    set<string> nodes;
    for(Edge E : edges){
        nodes.insert(E.from);
        nodes.insert(E.to);
    }


    vector<string> ancestors;
    ancestors.push_back("");
    vector<bool> bpr = {1}; // Open of root    
    for(string v : nodes){
        if(v == "") continue;

        while(!is_prefix_of(ancestors.back(),v)){
            ancestors.pop_back();
            bpr.push_back(0); // close
        }

        ancestors.push_back(v);
        bpr.push_back(1); // open
    }
    while(!ancestors.empty()){
        ancestors.pop_back();
        bpr.push_back(0); // close
    }

    return bpr;
}

string ST_to_dot(const vector<Edge>& edges, set<string> maxreps){
    stringstream dot;
    dot << "digraph suffixtree {\n";
    set<string> nodes;
    for(Edge E : edges){
        nodes.insert(E.from);
        nodes.insert(E.to);
    }
    map<string, int64_t> node_to_id;
    int64_t counter = 0;
    int64_t leaf_count = 0;
    for(string node : nodes){
        node_to_id[node] = counter;
        if(node.size() > 0 && node.back() == '$'){
            // Leaf: don't draw the label because it is so long
            dot << counter << " [label=\"" << leaf_count << "\"];\n";
            leaf_count++;
        } else{
            // Non-leaf. Draw label. Maxreps are red.
            if(maxreps.find(node) != maxreps.end()){
                dot << counter << " [style=filled, fillcolor=red, label=\"" + node + "\"];\n";
            } else{
                dot << counter << " [label=\"" + node + "\"];\n";
            }
        }
        counter++;
    }

    for(auto E : edges){
        int64_t from = node_to_id[E.from];
        int64_t to = node_to_id[E.to];
        dot << from << " -> " << to << ";\n";
    }
    dot << "}\n";
    return dot.str();
}

#endif

