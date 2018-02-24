#ifndef GLOBALS_HH
#define GLOBALS_HH

#include "sdsl/rank_support_v.hpp"
#include "sdsl/select_support_mcl.hpp"
#include "sdsl/bp_support_g.hpp"
#include "sdsl/io.hpp"
#include <sstream>
#include <string>
#include <vector>
#include <stdexcept>

using namespace std;

class Global_Data {

private:
	
public:
	
    sdsl::bit_vector slt_bpr;
    sdsl::bit_vector rev_st_bpr;
    sdsl::bit_vector rev_st_bpr_context_only;
    sdsl::bit_vector rev_st_maximal_marks;
    sdsl::bit_vector slt_maximal_marks;
    sdsl::bit_vector rev_st_context_marks;
    sdsl::bit_vector pruning_marks;

    sdsl::bp_support_g<> rev_st_bps;
    sdsl::bp_support_g<> slt_bps;
    sdsl::select_support_mcl<10,2> rev_st_ss_10;
    sdsl::rank_support_v<10,2> rev_st_rs_10;
    sdsl::rank_support_v<1> rev_st_maximal_marks_rs;
    sdsl::select_support_mcl<1> slt_maximal_marks_ss;
    sdsl::rank_support_v<1> rev_st_context_marks_rs;
    sdsl::select_support_mcl<1> rev_st_context_marks_ss;
    sdsl::bp_support_g<> rev_st_bpr_context_only_bps;
    sdsl::rank_support_v<1> pruning_marks_rs;
    sdsl::select_support_mcl<1> pruning_marks_ss;

	BD_BWT_index<> bibwt;

    std::vector<int64_t> string_depths; // Built only if used
	
	Global_Data(){}
	
	// Need custom copy functions because need to move the pointers in the support structures
	Global_Data(Global_Data const& other){
		*this = other;
	}
    Global_Data& operator=(Global_Data const& other){
		if(this != &other){
			this->slt_bpr = other.slt_bpr;
			this->rev_st_bpr = other.rev_st_bpr;
			this->rev_st_bpr_context_only = other.rev_st_bpr_context_only;
			this->rev_st_maximal_marks = other.rev_st_maximal_marks;
			this->slt_maximal_marks = other.rev_st_maximal_marks;
			this->rev_st_context_marks = other.rev_st_context_marks;
			this->pruning_marks = other.pruning_marks;
			
			this->rev_st_bps = other.rev_st_bps;
			this->slt_bps = other.slt_bps;
			this->rev_st_ss_10 = other.rev_st_ss_10;
			this->rev_st_rs_10 = other.rev_st_rs_10;
			this->rev_st_maximal_marks_rs = other.rev_st_maximal_marks_rs;
			this->slt_maximal_marks_ss = other.slt_maximal_marks_ss;
			this->rev_st_context_marks_rs = other.rev_st_context_marks_rs;
			this->rev_st_context_marks_ss = other.rev_st_context_marks_ss;
			this->rev_st_bpr_context_only_bps = other.rev_st_bpr_context_only_bps;
			this->pruning_marks_rs = pruning_marks_rs;
			this->pruning_marks_ss = pruning_marks_ss;
			
			rev_st_bps.set_vector(&rev_st_bpr);
			slt_bps.set_vector(&slt_bpr);
			rev_st_ss_10.set_vector(&rev_st_bpr);
			rev_st_rs_10.set_vector(&rev_st_bpr);
			rev_st_maximal_marks_rs.set_vector(&rev_st_bpr);
			slt_maximal_marks_ss.set_vector(&slt_bpr);
			rev_st_context_marks_rs.set_vector(&rev_st_bpr);
			rev_st_context_marks_ss.set_vector(&rev_st_bpr);
			rev_st_bpr_context_only_bps.set_vector(&rev_st_bpr_context_only);
			pruning_marks_rs.set_vector(&pruning_marks);
			pruning_marks_ss.set_vector(&pruning_marks);
		}
		
		return *this;
	}


    

    std::string toString() {
        std::stringstream ss;
        ss << "slt_bpr: " << slt_bpr << std::endl
           << "rev_st_bpr: " << rev_st_bpr << std::endl
           << "rev_st_bpr_context_only: " << rev_st_bpr_context_only << std::endl
           << "rev_st_maximal_marks: " << rev_st_maximal_marks << std::endl
           << "slt_maximal_marks: " << slt_maximal_marks << std::endl
           << "rev_st_context_marks: " << rev_st_context_marks << std::endl
           << "pruning_marks: " << pruning_marks << std::endl;
        return ss.str();
    }

    template<typename T>
    void store_to_disk(T& data, string path) {
        if(!sdsl::store_to_file(data, path)) {
            throw std::runtime_error("Error writing to disk " + path);
        }
    }

    template<typename T>
    void load_from_disk(T& data, string path) {
        //cout << sdsl::has_load<T>::value << " " << path << endl;
        if(!sdsl::load_from_file(data, path)) {
            throw std::runtime_error("Error reading from disk: " + path);
        }
    }


    // DOES NOT STORE string_depths!
    void store_all_to_disk(string directory, string filename_prefix) {

		bibwt.save_to_disk(directory, filename_prefix + "_bibwt");
        store_to_disk(slt_bpr, directory + "/" + filename_prefix + "_slt_bpr.dat");
        store_to_disk(rev_st_bpr, directory + "/" + filename_prefix + "_rev_st_bpr.dat");
        store_to_disk(rev_st_bpr_context_only, directory + "/" + filename_prefix + "_rev_st_bpr_context_only.dat");
        store_to_disk(rev_st_maximal_marks, directory + "/" + filename_prefix + "_rev_st_maximal_marks.dat");
        store_to_disk(slt_maximal_marks, directory + "/" + filename_prefix + "_slt_maximal_marks.dat");
        store_to_disk(rev_st_context_marks, directory + "/" + filename_prefix + "_rev_st_context_marks.dat");
        store_to_disk(pruning_marks, directory + "/" + filename_prefix + "_pruning_marks.dat");
        store_to_disk(rev_st_bps, directory + "/" + filename_prefix + "_rev_st_bps.dat");
        store_to_disk(rev_st_ss_10, directory + "/" + filename_prefix + "_rev_st_ss_10.dat");
        store_to_disk(rev_st_rs_10, directory + "/" + filename_prefix + "_rev_st_rs_10.dat");
        store_to_disk(rev_st_maximal_marks_rs, directory + "/" + filename_prefix + "_rev_st_maximal_marks_rs.dat");
        store_to_disk(rev_st_context_marks_ss, directory + "/" + filename_prefix + "_rev_st_context_marks_ss.dat");
        store_to_disk(slt_maximal_marks_ss, directory + "/" + filename_prefix + "_slt_maximal_marks_ss.dat");
        store_to_disk(rev_st_bpr_context_only_bps, directory + "/" + filename_prefix + "_rev_st_bpr_context_only_bps.dat");
        store_to_disk(pruning_marks_rs, directory + "/" + filename_prefix + "_pruning_marks_rs.dat");
        store_to_disk(pruning_marks_ss, directory + "/" + filename_prefix + "_pruning_marks_ss.dat");
        store_to_disk(slt_bps, directory + "/" + filename_prefix + "_slt_bps.dat");
    }

    // DOES NOT STORE string_depths!
    void load_all_from_disk(string directory, string filename_prefix) {

        bibwt.load_from_disk(directory, filename_prefix + "_bibwt");
        load_from_disk(slt_bpr, directory + "/" + filename_prefix + "_slt_bpr.dat");
        load_from_disk(rev_st_bpr, directory + "/" + filename_prefix + "_rev_st_bpr.dat");
        load_from_disk(rev_st_bpr_context_only, directory + "/" + filename_prefix + "_rev_st_bpr_context_only.dat");
        load_from_disk(rev_st_maximal_marks, directory + "/" + filename_prefix + "_rev_st_maximal_marks.dat");
        load_from_disk(slt_maximal_marks, directory + "/" + filename_prefix + "_slt_maximal_marks.dat");
        load_from_disk(rev_st_context_marks, directory + "/" + filename_prefix + "_rev_st_context_marks.dat");
        load_from_disk(pruning_marks, directory + "/" + filename_prefix + "_pruning_marks.dat");

        /*load_from_disk(rev_st_ss_10, directory + "/" + filename_prefix + "_rev_st_ss_10.dat");
        load_from_disk(rev_st_rs_10, directory + "/" + filename_prefix + "_rev_st_rs_10.dat");
        load_from_disk(rev_st_maximal_marks_rs, directory + "/" + filename_prefix + "_rev_st_maximal_marks_rs.dat");
        load_from_disk(slt_maximal_marks_ss, directory + "/" + filename_prefix + "_slt_maximal_marks_ss.dat");
        load_from_disk(rev_st_context_marks_rs, directory + "/" + filename_prefix + "_rev_st_context_marks_rs.dat");
        load_from_disk(pruning_marks_rs, directory + "/" + filename_prefix + "_pruning_marks_rs.dat");
        load_from_disk(pruning_marks_ss, directory + "/" + filename_prefix + "_pruning_marks_ss.dat"); */

        ifstream in;
        in.exceptions ( ifstream::failbit | ifstream::badbit );
        try {
            in.open(directory + "/" + filename_prefix + "_rev_st_bps.dat");
            rev_st_bps.load(in, &rev_st_bpr);
            in.close();
            
            in.open(directory + "/" + filename_prefix + "_slt_bps.dat");
            slt_bps.load(in, &slt_bpr);
            in.close();
            
            in.open(directory + "/" + filename_prefix + "_rev_st_bpr_context_only_bps.dat");
            rev_st_bpr_context_only_bps.load(in, &rev_st_bpr_context_only);
            in.close();
            
            in.open(directory + "/" + filename_prefix + "_rev_st_ss_10.dat");
            rev_st_ss_10.load(in, &rev_st_bpr);
            in.close();
            
            in.open(directory + "/" + filename_prefix + "_rev_st_rs_10.dat");
            rev_st_rs_10.load(in, &rev_st_bpr);
            in.close();
            
            in.open(directory + "/" + filename_prefix + "_rev_st_maximal_marks_rs.dat");
            rev_st_maximal_marks_rs.load(in, &rev_st_maximal_marks);
            in.close();
            
            in.open(directory + "/" + filename_prefix + "_slt_maximal_marks_ss.dat");
            slt_maximal_marks_ss.load(in, &slt_maximal_marks);
            in.close();
            
            in.open(directory + "/" + filename_prefix + "_rev_st_context_marks_rs.dat");
            rev_st_context_marks_rs.load(in, &rev_st_context_marks);
            in.close();
            
            in.open(directory + "/" + filename_prefix + "_rev_st_context_marks_ss.dat");
            rev_st_context_marks_ss.load(in, &rev_st_context_marks);
            in.close();
            
            in.open(directory + "/" + filename_prefix + "_pruning_marks_rs.dat");
            pruning_marks_rs.load(in, &pruning_marks);
            in.close();
            
            in.open(directory + "/" + filename_prefix + "_pruning_marks_ss.dat");
            pruning_marks_ss.load(in, &pruning_marks);
            in.close();
        } catch(ifstream::failure e) {
            cout << "Error loading data structure from disk" << endl;
            exit(-1);
            // Todo: filename
        }
    }

};

#endif
