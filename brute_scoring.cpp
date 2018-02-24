#include <iostream>
#include <string>
#include <vector>
#include <set>
#include <cmath>

using namespace std;


// Computes the right-side empirical entropy of W in T
double entropy(string W, string T){
    T += "$";
    vector<int64_t> right_counts(256);
    int64_t right_extensions = 0;
    for(int64_t i = 0; i < T.size(); i++){
        if(T.substr(i,W.size()) == W){
            right_counts[T[i+W.size()]]++;
            right_extensions++;
        }
    }
    double entropy = 0;
    for(int64_t i = 0; i < 256; i++){
        if(right_counts[i] != 0){
            double prob = ((double)right_counts[i] / right_extensions);
            entropy -= prob * log2(prob);
        }
    }
    return entropy;
}

int64_t count(string W, string T){
    int64_t ans = 0;
    for(int64_t i = 0; i < T.size(); i++){
        if(T.substr(i,W.size()) == W){
            ans++;
        }
    }
    return ans;
}

// W is a context iff f(W)H(W) - \sum_{a \in \Sigma} f(aW) H(aW) >= threshold
// where f is the number of occurrences of W and H(W) is the empirical entropy of W
// to the right, i.e. if P(b) = f(Wb) / f(W), then
// H(W) = - \sum_{b \in \Sigma} P(b) log P(b)
vector<string> get_contexts_entropy(const string& T, double threshold){
    set<char> alphabet;
    for(char c : T) alphabet.insert(c);
    set<string> contexts = {""};
    for(int64_t start = 0; start < T.size(); start++){
        for(int64_t length = 1; length <= T.size()-start; length++){
            string W = T.substr(length);
            double lefthand = count(W,T) * entropy(W,T);
            for(char a : alphabet){
                lefthand -= count(a + W, T) * entropy(a + W, T); 
            }
            if(lefthand >= threshold)
                contexts.insert(W);
        }
    }
    vector<string> ans(contexts.begin(), contexts.end());
    return ans;
}

double score_string_brute(string S, string T){
    return 0;
}

int main(){
    string T = "abababbababbbabababababbabbbabaabbabababbabababaaababab";
    for(string W : get_contexts_entropy(T,0.1)) cout << W << endl;
}
