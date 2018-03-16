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
#include <sstream>
#include <string>
#include <vector>
#include <stdexcept>

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
    std::unique_ptr<BWT> revbwt; // Constructed during build time, used during scoring time.

    std::vector<int64_t> string_depths; // Built only if used

    Global_Data() {}

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

    // DOES NOT STORE string_depths!
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
        
    }

    // DOES NOT STORE string_depths!
    void load_all_from_disk(string directory, string filename_prefix, bool run_length_coding, bool load_bibwt = false) {

        if(load_bibwt){
            bibwt = make_shared<BD_BWT_index<>>(); // todo: RLE?
            bibwt->load_from_disk(directory, filename_prefix + ".bibwt");
        }
        if(run_length_coding){
            std::unique_ptr<BWT> bwt_object(new RLEBWT<>());
            revbwt = std::move(bwt_object);
            
            slt_bpr = shared_ptr<RLE_bitvector>(new RLE_bitvector());
            slt_maximal_marks = shared_ptr<RLE_bitvector>(new RLE_bitvector());
            pruning_marks = shared_ptr<RLE_bitvector>(new RLE_bitvector());
        } else{
            std::unique_ptr<BWT> bwt_object(new Basic_BWT<>());
            revbwt = std::move(bwt_object);
            
            slt_bpr = shared_ptr<Basic_bitvector>(new Basic_bitvector());
            slt_maximal_marks = shared_ptr<Basic_bitvector>(new Basic_bitvector());
            pruning_marks = shared_ptr<Basic_bitvector>(new Basic_bitvector());
        }
        
        rev_st_bpr = shared_ptr<Basic_bitvector>(new Basic_bitvector());
        rev_st_bpr_context_only = shared_ptr<Basic_bitvector>(new Basic_bitvector());
        rev_st_maximal_marks = shared_ptr<Basic_bitvector>(new Basic_bitvector());
        rev_st_context_marks = shared_ptr<Basic_bitvector>(new Basic_bitvector());
        
        revbwt->load_from_disk(directory, filename_prefix + ".rev_bwt");
        
        slt_bpr->load(directory + "/" + filename_prefix + ".slt_bpr");
        rev_st_bpr->load(directory + "/" + filename_prefix + ".rev_st_bpr");
        rev_st_bpr_context_only->load(directory + "/" + filename_prefix + ".rev_st_bpr_context_only");
        rev_st_maximal_marks->load(directory + "/" + filename_prefix + ".rev_st_maximal_marks");
        slt_maximal_marks->load(directory + "/" + filename_prefix + ".slt_maximal_marks");
        rev_st_context_marks->load(directory + "/" + filename_prefix + ".rev_st_context_marks");
        pruning_marks->load(directory + "/" + filename_prefix + ".pruning_marks");

    }

};

#endif
