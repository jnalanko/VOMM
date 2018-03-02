//
//  UniBWT.h
//  VOMM
//
//  Created by Niklas Alanko on 25/02/2018.
//  Copyright © 2018 University of Helsinki. All rights reserved.
//

#ifndef UniBWT_h
#define UniBWT_h

#include <sdsl/construct.hpp>
#include <vector>
#include <utility>
#include <string>
#include "bwt.hh"
#include "Interval.hh"
#include "sdsl/io.hpp"
#include "Interfaces.hh"
#include <stdexcept>


/*
 * Implements an FM.index for a byte alphabet.
 * All indices and ranks are indexed starting from zero.
 * All intervals are inclusive, i.e. the interval from i to j includes both i and j.
 *
 * Terminology: C-array of an interval is an array with size 256 such that
 * C[c] is the number of characters with lexicographical rank strictly less than c at the given interval.
 */

template<class t_bitvector = sdsl::bit_vector>
class Basic_BWT : public BWT {
    
public:
    
    
    
private:
    
    sdsl::wt_huff<t_bitvector> bwt;
    
    std::vector<int64_t> global_c_array;
    std::vector<uint8_t> alphabet;
    
    std::vector<uint8_t> get_string_alphabet(const uint8_t* s) const;
    int64_t strlen(const uint8_t* str) const;
    void count_smaller_chars(const sdsl::wt_huff<t_bitvector>& bwt, std::vector<int64_t>& counts, Interval I);
    
public:
    
    // The input string must not contain the END byte
    static const uint8_t END = 0x01; // End of string marker.
    Basic_BWT() {}
    Basic_BWT(const uint8_t* input);
    
    uint8_t get_END() const { return END; }
    int64_t size() const { return bwt.size();}
    uint8_t bwt_at(int64_t index) const { return bwt[index]; }
    const std::vector<int64_t>& get_global_c_array() const { return global_c_array; }
    const std::vector<uint8_t>& get_alphabet() const { return alphabet; }
    
    // Results are stored in the provided struct reference
    void compute_interval_data(Interval I, Interval_Data& data){
        if(I.size() == 0)
            data.n_distinct_symbols = 0;
        else
            bwt.interval_symbols(I.left, I.right+1, data.n_distinct_symbols, data.symbols, data.ranks_start, data.ranks_end);
    }
    
   Interval search(Interval I, uint8_t c){
        if(I.size() == 0)
            return Interval(-1,-2);
        
        // Todo: rank at I.left is computed twice. Reuse.
        int64_t num_c_in_interval = bwt.rank(I.right + 1,c) - bwt.rank(I.left,c);
        int64_t start_new = get_global_c_array()[c] + bwt.rank(I.left, c);
        int64_t end_new = start_new + num_c_in_interval - 1;
        
        if(start_new > end_new) return Interval(-1,-2); // num_c_in_interval == 0
        
        return Interval(start_new,end_new);
    }
    
    Interval search_with_precalc(Interval I, uint8_t c, Interval_Data& D){
        (void)I; (void) c; (void) D;
        throw(std::runtime_error("Not implemented error: search_precomputed_data"));
        /*
        if(I.size() == 0)
            return Interval_pair(-1,-2,-1,-2);
        
        // Todo: rank at I.left is computed twice. Reuse.
        int64_t num_c_in_interval = bwt.rank(I.right + 1,c) - bwt.rank(I.left,c);
        int64_t start_new = get_global_c_array()[c] + bwt.rank(I.left, c);
        int64_t end_new = start_new + num_c_in_interval - 1;
        
        if(start_new > end_new) return Interval_pair(-1,-2,-1,-2); // num_c_in_interval == 0
        
        return Interval(start_new,end_new);
        */
    }
    
    void save_to_disk(std::string directory, std::string filename_prefix);
    void load_from_disk(std::string directory, std::string filename_prefix);
    
};

template<class t_bitvector>
const uint8_t Basic_BWT<t_bitvector>::END;

template<class t_bitvector>
void Basic_BWT<t_bitvector>::count_smaller_chars(const sdsl::wt_huff<t_bitvector>& bwt,
                                                    std::vector<int64_t>& counts, Interval I){
    assert(alphabet.size() != 0);
    for(int64_t i = 0; i < (int64_t)alphabet.size(); i++) counts[alphabet[i]] = 0;
    if(I.size() == 0) return;
    
    sdsl::int_vector_size_type nSymbols;
    std::vector<uint8_t> symbols(256);
    std::vector<uint64_t> ranks_i(256);
    std::vector<uint64_t> ranks_j(256);
    bwt.interval_symbols(I.left, I.right+1, nSymbols, symbols, ranks_i, ranks_j);
    
    // Put in the counts
    for(int64_t k = 0; k < nSymbols; k++){
        counts[symbols[k]] = ranks_j[k] - ranks_i[k];
    }
    
    // Do the cumulative sum
    int64_t prev_cumul = 0;
    int64_t prev_count = counts[alphabet[0]];
    int64_t temp = 0;
    for(int64_t k = 0; k < alphabet.size(); k++){
        temp = counts[alphabet[k]];
        if(k == 0) counts[alphabet[k]] = 0;
        else{
            counts[alphabet[k]] = prev_cumul + prev_count;
        }
        prev_count = temp;
        prev_cumul = counts[alphabet[k]];
    }
}

// Returns the alphabet in sorted order
template<class t_bitvector>
std::vector<uint8_t> Basic_BWT<t_bitvector>::get_string_alphabet(const uint8_t* s) const{
    
    std::vector<bool> found(256,false);
    while(*s != 0){
        found[*s] = true;
        s++;
    }
    
    std::vector<uint8_t> alphabet;
    for(int i = 0; i < 256; i++){
        if(found[i]) alphabet.push_back((uint8_t)i);
    }
    
    return alphabet;
}

// strlen(const uint8_t*) is not in the standard library
template<class t_bitvector>
int64_t Basic_BWT<t_bitvector>::strlen(const uint8_t* str) const{
    const uint8_t* start = str;
    while(*str != 0) str++;
    return str - start;
}

template<class t_bitvector>
Basic_BWT<t_bitvector>::Basic_BWT(const uint8_t* input)
: global_c_array(256) {
    if(*input == 0) throw std::runtime_error("Tried to construct BD_BWT_index for an empty string");
    int64_t n = strlen(input);
    
    if(std::find(input, input+n, END) != input + n){
        std::stringstream error;
        error << "Input string contains forbidden byte " << std::hex << END;
        throw std::runtime_error(error.str());
    }
    
    // Build the bwt
    uint8_t* data = (uint8_t*) malloc(sizeof(uint8_t) * (n + 1));
    for(int64_t i = 0; i < n; i++){
        data[i] = input[i];
    }
    data[n] = END;
    
    uint8_t* data_bwt = bwt_dbwt(data,n,END);
    free(data);
    
    // Build wavelet trees
    construct_im(this->bwt, (const char*)data_bwt, 1); // Must cast to signed char* or else breaks. File a bug report to sdsl?
    
    this->alphabet = get_string_alphabet(data_bwt);
    
    free(data_bwt);
    
    // Compute cumulative character counts
    count_smaller_chars(bwt,global_c_array,Interval(0,bwt.size()-1));
}


template<class t_bitvector>
void Basic_BWT<t_bitvector>::save_to_disk(std::string directory, std::string filename_prefix){
    std::string bwt_path = directory + "/" + filename_prefix + "_bwt.dat";
    if(!sdsl::store_to_file(bwt, bwt_path)) {
        throw std::runtime_error("Error writing to disk: " + bwt_path);
    }
    
    // Copy to sdsl bit vector because they have serialization built in
    sdsl::int_vector<64> global_c_array_sdsl(global_c_array.size());
    for(int64_t i = 0; i < global_c_array.size(); i++){
        global_c_array_sdsl[i] = global_c_array[i];
    }
    std::string gca = directory + "/" + filename_prefix + "_gca.dat";
    if(!sdsl::store_to_file(global_c_array_sdsl, gca)) {
        throw std::runtime_error("Error writing to disk: " + gca);
    }
    
    // Copy to sdsl bit vector because they have serialization built in
    sdsl::int_vector<8> alphabet_sdsl(alphabet.size());
    for(int64_t i = 0; i < alphabet.size(); i++){
        alphabet_sdsl[i] = alphabet[i];
    }
    std::string A = directory + "/" + filename_prefix + "_alphabet.dat";
    if(!sdsl::store_to_file(alphabet_sdsl, A)) {
        throw std::runtime_error("Error writing to disk: " + A);
    }
}

template<class t_bitvector>
void Basic_BWT<t_bitvector>::load_from_disk(std::string directory, std::string filename_prefix){
    global_c_array.resize(256);
    
    std::string bwt_path = directory + "/" + filename_prefix + "_bwt.dat";
    if(!sdsl::load_from_file(bwt, bwt_path)) {
        throw std::runtime_error("Error reading from disk: " + bwt_path);
    }
    
    sdsl::int_vector<64> gca_sdsl;
    std::string gca_path = directory + "/" + filename_prefix + "_gca.dat";
    if(!sdsl::load_from_file(gca_sdsl, gca_path)) {
        throw std::runtime_error("Error reading from disk: " + gca_path);
    }
    for(int64_t i = 0; i < 256; i++){
        global_c_array[i] = gca_sdsl[i];
    }
    
    sdsl::int_vector<8> alphabet_sdsl;
    std::string alphabet_path = directory + "/" + filename_prefix + "_alphabet.dat";
    if(!sdsl::load_from_file(alphabet_sdsl, alphabet_path)) {
        throw std::runtime_error("Error reading from disk: " + alphabet_path);
    }
    
    alphabet.resize(alphabet_sdsl.size());
    for(int64_t i = 0; i < alphabet_sdsl.size(); i++){
        alphabet[i] = alphabet_sdsl[i];
    }
}




#endif /* UniBWT_h */
