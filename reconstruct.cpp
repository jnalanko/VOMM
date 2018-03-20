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
#include "score_string.hh"
#include "build_model.hh"
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

class Reconstruction_Config{
    
private:
    
  Reconstruction_Config(const Reconstruction_Config&); // Prevent copy-construction
  Reconstruction_Config& operator=(const Reconstruction_Config&);  // Prevent assignment
    
public:
        
    bool context_stats;
    bool only_maxreps;
    bool run_length_coding;
    
    Context_Callback* cf;

    string modeldir;
    string filename;
    
    Reconstruction_Config() : context_stats(false), only_maxreps(false), run_length_coding(false), cf(nullptr) {}
    
    void assert_all_ok(){
        assert(modeldir != "");
        assert(filename != "");
        assert(cf != nullptr);
    }
    
    void load_info_file(){
        assert(modeldir != "");
        assert(filename != "");
        string path = modeldir + "/" + filename + ".info";
        ifstream file(path);
        string ctype; // Old context type. Unused
        file >> only_maxreps >> ctype >> run_length_coding;
        if(!file.good()){
            cerr << "Error reading file: " << path << endl;
            exit(-1);
        }
    }
    
};

int score_string_main(int argc, char** argv){
    
    if(argc < 4){
        cerr << "Marks contexts again. Usage: see readme." << endl;
        return -1;
    }
    
    Reconstruction_Config C;
    
    for(int64_t i = 1; i < argc; i++){
        if(argv[i] == string("--dir")){
            i++;
            C.modeldir = argv[i];
        } else if(argv[i] == string("--file")){
            i++;
            C.filename = argv[i];
        } else if(argv[i] == string("--entropy")){
            i++;
            double threshold = stod(argv[i]);
            C.cf = new Entropy_Formula(threshold);
        } else if(argv[i] == string("--KL")){
            i++;
            double threshold = stod(argv[i]);
            C.cf = new KL_Formula(threshold);
        } else if(argv[i] == string("--pnorm")){
            i++;
            int64_t p = stoi(argv[i]);
            i++;
            double threshold = stod(argv[i]);
            C.cf = new pnorm_Formula(p, threshold);
        } else if(argv[i] == string("--four-thresholds")){
            double t1,t2,t3,t4;
            i++; t1 = stod(argv[i]);
            i++; t2 = stod(argv[i]);
            i++; t3 = stod(argv[i]);
            i++; t4 = stod(argv[i]);
            C.cf = new EQ234_Formula(t1,t2,t3,t4);
        } else if(argv[i] == string("--context-stats")){
            C.context_stats = true;
        } else{
            cerr << "Invalid argument: " << argv[i] << endl;
            return -1;
        }
    }
    
    C.load_info_file();
    C.assert_all_ok();
    
    write_log("Loading the model from " + C.modeldir);
    Global_Data G;
    G.load_all_from_disk(C.modeldir, C.filename, C.run_length_coding, true);
    write_log("Starting to rebuild contexts");
    
    SLT_Iterator iterator(G.bibwt.get());
    Pruned_Topology_Mapper mapper(G.rev_st_bpr, G.pruning_marks);
    Stats_writer wr;
    if(C.context_stats){
        wr.set_file(C.modeldir + "/stats.depths_and_scores.txt");
    }
    C.cf->init(G.bibwt.get(), G.rev_st_bpr->size(), mapper, &wr);
    iterate_with_callback(iterator, C.cf);
    G.rev_st_context_marks = make_shared<Basic_bitvector>(C.cf->get_result());
    G.rev_st_context_marks->init_rank_support();
    G.rev_st_context_marks->init_select_support();
    
    if(C.context_stats){ 
        write_context_summary(G, C.cf->get_number_of_candidates(), C.modeldir + "/stats.context_summary.txt");
    }
    
    G.store_all_to_disk(C.modeldir, C.filename);
    
    write_log("Done");
    
    return 0;
    
}

int main(int argc, char** argv){
    return score_string_main(argc,argv);
}
