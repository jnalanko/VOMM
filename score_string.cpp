//
//  score_string.cpp
//  PST
//
//  Created by Niklas Alanko on 24/04/2017.
//  Copyright © 2017 University of Helsinki. All rights reserved.
//

#include <iostream>
#include <string>
#include <sstream>
#include <fstream>
#include <vector>
#include <ctime>
#include <chrono>
#include <streambuf>
#include "String_Depth_Support.hh"
#include "Parent_Support.hh"
#include "LMA_Support.hh"
#include "BPR_tools.hh"
#include "Precalc.hh"
#include "FASTA_parsing.hh"
#include "score_string.hh"
#include "logging.hh"

using namespace std;

string read_raw_file(string filename){

    std::ifstream t(filename.c_str());
    std::string str;

    t.seekg(0, std::ios::end);   
    str.reserve(t.tellg());
    t.seekg(0, std::ios::beg);

    str.assign((std::istreambuf_iterator<char>(t)),
                std::istreambuf_iterator<char>());
    
    return str;
}

class Scoring_Config{
    
private:
    
  Scoring_Config(const Scoring_Config&); // Prevent copy-construction
  Scoring_Config& operator=(const Scoring_Config&);  // Prevent assignment
    
public:
    
    enum Context_Type {UNDEFINED, EQ234, ENTROPY, KL, PNORM};
    
    bool only_maxreps;
    Context_Type context_type;
    double escapeprob;
    string modeldir;
    string filename;
    
    Scoring_Function* scorer = nullptr;
    Loop_Invariant_Updater* updater = nullptr;
    
    Scoring_Config() : only_maxreps(false), context_type(UNDEFINED), escapeprob(-1), scorer(nullptr), updater(nullptr) {}
    
    ~Scoring_Config(){
        delete scorer;
        delete updater;
    }
    
    void assert_all_ok(){
        assert(modeldir != "");
        assert(filename != "");
        assert(context_type != UNDEFINED);
        assert(escapeprob != -1);
        assert(scorer != nullptr);
        assert(updater != nullptr);
    }
    
    void load_info_file(){
        assert(modeldir != "");
        assert(filename != "");
        string path = modeldir + "/" + filename + ".info";
        ifstream file(path);
        string ctype;
        file >> only_maxreps >> ctype;
        if(!file.good()){
            cerr << "Error reading file: " << path << endl;
            exit(-1);
        }
        
        if(ctype == "EQ234") context_type = EQ234;
        else if(ctype == "entropy") context_type = ENTROPY;
        else if(ctype == "KL") context_type = KL;
        else if(ctype == "pnorm") context_type = PNORM;
        else assert(false);
        
    }
    
};

int main(int argc, char** argv){
    if(argc < 4){
        cerr << "Computes the probability of string S against string T" << endl;
        cerr << "Usage: ./score_string [--query-fasta or --query-raw] S.txt <rest of the options>" << endl;
        return -1;
    }
    
    Scoring_Config C;
    

    
    vector<string> queries;
    string reference;
    for(int64_t i = 1; i < argc; i++){
        if(argv[i] == string("--query-raw")){
            i++;
            queries.push_back(read_raw_file(argv[i]));
        } else if(argv[i] == string("--query-fasta")){
            i++;
            auto v = parse_FASTA(argv[i]);
            for(auto pair : v) queries.push_back(pair.first);
        }
        else if(argv[i] == string("--dir")){
            i++;
            C.modeldir = argv[i];
        } else if(argv[i] == string("--file")){
            i++;
            C.filename = argv[i];
        } else if(argv[i] == string("--escapeprob")){
            i++;
            C.escapeprob = stod(argv[i]);
        }
        else{
            cerr << "Invalid argument: " << argv[i] << endl;
            return -1;
        }
    }
    
    C.load_info_file();
    
    if(C.context_type == Scoring_Config::ENTROPY){
        // Type W
        C.scorer = new Basic_Scorer(C.escapeprob, true);
    }
    if(C.context_type != Scoring_Config::ENTROPY){
        // Type aW
        C.scorer = new Basic_Scorer(C.escapeprob, false);
    }
    
    if(C.only_maxreps){
        C.updater = new Maxrep_Depth_Bounded_Updater();
    } else{
        // No pruning at all (model can't have depth pruning + no maxrep pruning, checked by assert_all_ok in model buildnig)
        C.updater = new Basic_Updater();
    }

    C.assert_all_ok();
    
    write_log("Loading the model from " + C.modeldir);
    Global_Data G;
    G.load_all_from_disk(C.modeldir, C.filename);
            
    score_string(reference, G, *C.scorer, *C.updater);
          
    write_log("Starting to process queries");
    for(string query : queries){
        cout << score_string(query, G, *C.scorer, *C.updater) << endl;
    }
    write_log("Done");
    
}
