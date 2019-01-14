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
#include "input_reading.hh"
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
    
    enum class Context_Type {UNDEFINED, EQ234, ENTROPY, KL, PNORM};
    enum class Input_Mode {UNDEFINED, RAW, FASTA};
    
    Input_Mode input_mode;
    string query_filename;
    bool only_maxreps;
    Context_Type context_type;
    double escapeprob;
    string modeldir;
    string reference_filename;
    bool run_length_coding;
    bool recursive_fallback;
    bool lin_scoring;
    int64_t depth_bound;
    
    Scoring_Function* scorer;
    Loop_Invariant_Updater* updater;
    
    Scoring_Config() : input_mode(Input_Mode::UNDEFINED), only_maxreps(false), context_type(Context_Type::UNDEFINED), 
    escapeprob(-1), run_length_coding(false), recursive_fallback(false), lin_scoring(false), depth_bound(-1), scorer(nullptr), updater(nullptr)
     {}
    
    ~Scoring_Config(){
        delete scorer;
        delete updater;
    }
    
    void assert_all_ok(){
        assert(modeldir != "");
        assert(reference_filename != "");
        assert(query_filename != "");
        assert(input_mode != Input_Mode::UNDEFINED);
        assert(context_type != Context_Type::UNDEFINED);
        if(!lin_scoring) assert(escapeprob != -1);
        assert(scorer != nullptr);
        assert(updater != nullptr);
        assert(depth_bound != -1);
    }
    
    void load_info_file(){
        assert(modeldir != "");
        assert(reference_filename != "");
        string path = modeldir + "/" + reference_filename + ".info";
        ifstream file(path);
        string ctype;
        file >> only_maxreps >> ctype >> run_length_coding >> depth_bound;
        if(!file.good()){
            cerr << "Error reading file: " << path << endl;
            exit(-1);
        }
        
        if(ctype == "EQ234") context_type = Context_Type::EQ234;
        else if(ctype == "entropy") context_type = Context_Type::ENTROPY;
        else if(ctype == "KL") context_type = Context_Type::KL;
        else if(ctype == "pnorm") context_type = Context_Type::PNORM;
        else assert(false);
        
    }
    
};

int main(int argc, char** argv){
    if(argc < 4){
        cerr << "Computes the probability of string against a VOMM index" << endl;
        cerr << "Usage: see README.md" << endl;
        return -1;
    }
    
    Scoring_Config C;
    for(int64_t i = 1; i < argc; i++){
        if(argv[i] == string("--query-raw")){
            i++;
            C.query_filename = argv[i];
            C.input_mode = Scoring_Config::Input_Mode::RAW;
        } else if(argv[i] == string("--query-fasta")){
            i++;
            C.query_filename = argv[i];
            C.input_mode = Scoring_Config::Input_Mode::FASTA;
        }
        else if(argv[i] == string("--dir")){
            i++;
            C.modeldir = argv[i];
        } else if(argv[i] == string("--file")){
            i++;
            C.reference_filename = argv[i];
        } else if(argv[i] == string("--escapeprob")){
            i++;
            C.escapeprob = stod(argv[i]);
        } else if(argv[i] == string("--recursive-fallback")){
            C.recursive_fallback = true;
        } else if(argv[i] == string("--lin-scoring")){
            C.lin_scoring = true;
        } else{
            cerr << "Invalid argument: " << argv[i] << endl;
            return -1;
        }
    }
    
    C.load_info_file();
    
    if(C.recursive_fallback){
        C.scorer = new Recursive_Scorer(C.escapeprob, (C.context_type == Scoring_Config::Context_Type::ENTROPY));
    } else{
        C.scorer = new Basic_Scorer(C.escapeprob, (C.context_type == Scoring_Config::Context_Type::ENTROPY));
    }
    
    if(C.only_maxreps){
        C.updater = new Maxrep_Pruned_Updater();
    } else{
        // No pruning at all (model can't have depth pruning + no maxrep pruning, checked by assert_all_ok in model buildnig)
        C.updater = new Basic_Updater();
    }

    C.assert_all_ok();
    
    write_log("Loading the model from " + C.modeldir);
    Global_Data G;
    if(C.lin_scoring)
        G.load_structures_that_lin_scoring_needs(C.modeldir, C.reference_filename);
    else
        G.load_all_from_disk(C.modeldir, C.reference_filename, false);
    write_log("Starting to score ");
        
    if(C.input_mode == Scoring_Config::Input_Mode::RAW){
        Raw_file_stream rfs(C.query_filename);
        if(C.lin_scoring){
            cout << score_string_lin(rfs,G) << endl;
        } else{
            cout << score_string(rfs, G, *C.scorer, *C.updater) << endl;
        }
    }
    
    if(C.input_mode == Scoring_Config::Input_Mode::FASTA){
        FASTA_reader fr(C.query_filename);
        while(!fr.done()){
            Read_stream input = fr.get_next_query_stream();
            if(C.lin_scoring){
                cout << score_string_lin(input,G) << endl;
            } else{
                cout << score_string(input, G, *C.scorer, *C.updater) << endl;
            }
        }
    }
    
    write_log("Done");
    
}
