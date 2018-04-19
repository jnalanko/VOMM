#ifndef GLOBALS_HH
#define GLOBALS_HH

#include "sdsl/rank_support_v.hpp"
#include "sdsl/select_support_mcl.hpp"
#include "sdsl/bp_support_g.hpp"
#include "sdsl/io.hpp"
#include <memory>
#include "UniBWT.h"
#include "RLEBWT.hh"
#include "Basic_bitvector.hh"
#include "RLE_bitvector.hh"
#include "All_Ones_Bitvector.hh"
#include <sstream>
#include <string>
#include <vector>
#include <stdexcept>
#include <memory>

using namespace std;


// This class constains all the data needed to represent the markov model in a succinct indexed form.
class Global_Data {

    private:
        
    // Forbid copying
    Global_Data(Global_Data const& other);
    Global_Data& operator=(Global_Data const& other);

    public:
        
    
    std::shared_ptr<Bitvector> slt_bpr; // O(n)
    std::shared_ptr<Bitvector> rev_st_bpr; // O(maxreps)
    std::shared_ptr<Bitvector> rev_st_bpr_context_only; // O(maxreps)
    std::shared_ptr<Bitvector> rev_st_maximal_marks; // O(maxreps)
    std::shared_ptr<Bitvector> slt_maximal_marks; // O(n)
    std::shared_ptr<Bitvector> rev_st_context_marks; // O(maxreps)
    std::shared_ptr<Bitvector> pruning_marks; // O(n)

    std::shared_ptr<BIBWT> bibwt; // Used for construction and reconstruction
    std::shared_ptr<BWT> revbwt; // Constructed during build time, used during scoring time.

    std::shared_ptr<sdsl::int_vector<0>> string_depths; // Built only if used

    Global_Data() {}
    
    bool have_slt(){
        assert(slt_bpr != nullptr);
        return slt_bpr->size() != 0;
    }

    std::string toString() {
        std::stringstream ss;
        ss << "slt_bpr: " << slt_bpr->toString() << std::endl
        << "rev_st_bpr: " << rev_st_bpr->toString() << std::endl
        << "rev_st_bpr_context_only: " << rev_st_bpr_context_only->toString() << std::endl
        << "rev_st_maximal_marks: " << rev_st_maximal_marks->toString() << std::endl
        << "slt_maximal_marks: " << slt_maximal_marks->toString() << std::endl
        << "rev_st_context_marks: " << rev_st_context_marks->toString() << std::endl
        << "pruning_marks: " << pruning_marks->toString() << std::endl;
        return ss.str();
    }

    void store_all_to_disk(string directory, string filename_prefix) {
        revbwt->save_to_disk(directory, filename_prefix + ".rev_bwt");
        bibwt->save_to_disk(directory, filename_prefix + ".bibwt");
        
        slt_bpr->serialize(directory + "/" + filename_prefix + ".slt_bpr");
        rev_st_bpr->serialize(directory + "/" + filename_prefix + ".rev_st_bpr");
        rev_st_bpr_context_only->serialize(directory + "/" + filename_prefix + ".rev_st_bpr_context_only");
        rev_st_maximal_marks->serialize(directory + "/" + filename_prefix + ".rev_st_maximal_marks");
        slt_maximal_marks->serialize(directory + "/" + filename_prefix + ".slt_maximal_marks");
        rev_st_context_marks->serialize(directory + "/" + filename_prefix + ".rev_st_context_marks");
        pruning_marks->serialize(directory + "/" + filename_prefix + ".pruning_marks");
        
        store_to_file(*string_depths, directory + "/" + filename_prefix + ".string_depths");
        
    }
    
    void load_bitvector(std::shared_ptr<Bitvector>& destination, string path){
        string type;
        ifstream info;
        info.exceptions(ifstream::failbit | ifstream::badbit);
        info.open(path + "_info");
        info >> type;
        info.close();
        if(type == "basic"){
            destination = make_shared<Basic_bitvector>();
        } else if(type == "rle"){
            destination = make_shared<RLE_bitvector>();
        } else if(type == "all-ones"){
            destination = make_shared<All_Ones_Bitvector>();
        } else {
            throw(std::runtime_error("Unknown bit vector type: " + type));    
        }
        destination->load(path);
    }
    
    void load_bwt(std::shared_ptr<BWT>& destination, string directory, string filename_prefix){
        string path = directory + "/" +  filename_prefix + ".rev_bwt";
        string type;
        ifstream info;
        info.exceptions(ifstream::failbit | ifstream::badbit);
        info.open(path + "_bwt_info");
        info >> type;
        info.close();
        if(type == "basic_bwt"){
            destination = make_shared<Basic_BWT<>>();
        } else if(type == "rle_bwt"){
            destination = make_shared<RLEBWT<>>();            
        } else{
            throw(std::runtime_error("Unknown BWT type: " + type));    
        }
        destination->load_from_disk(directory, filename_prefix + ".rev_bwt");
    }

    void load_all_from_disk(string directory, string filename_prefix, bool load_bibwt) {

        if(load_bibwt){
            bibwt = make_shared<BD_BWT_index<>>(); // todo: RLE?
            bibwt->load_from_disk(directory, filename_prefix + ".bibwt");
        }
        
        load_bwt(revbwt, directory, filename_prefix);
        
        load_bitvector(slt_bpr, directory + "/" + filename_prefix + ".slt_bpr");
        load_bitvector(rev_st_bpr, directory + "/" + filename_prefix + ".rev_st_bpr");
        load_bitvector(rev_st_bpr_context_only, directory + "/" + filename_prefix + ".rev_st_bpr_context_only");
        load_bitvector(rev_st_maximal_marks, directory + "/" + filename_prefix + ".rev_st_maximal_marks");
        load_bitvector(slt_maximal_marks, directory + "/" + filename_prefix + ".slt_maximal_marks");
        load_bitvector(rev_st_context_marks, directory + "/" + filename_prefix + ".rev_st_context_marks");
        load_bitvector(pruning_marks, directory + "/" + filename_prefix + ".pruning_marks");
        
        string_depths = shared_ptr<sdsl::int_vector<0>>(new sdsl::int_vector<0>());
        load_from_file(*string_depths, directory + "/" + filename_prefix + ".string_depths");
        
    }
    
    void load_structures_that_lin_scoring_needs(string directory, string filename_prefix){
       load_bwt(revbwt, directory, filename_prefix);
       load_bitvector(rev_st_bpr, directory + "/" + filename_prefix + ".rev_st_bpr");
       load_bitvector(pruning_marks, directory + "/" + filename_prefix + ".pruning_marks");
    }

};

#endif
