#ifndef RLEBWT_h
#define RLEBWT_h

#include <sdsl/construct.hpp>
#include <vector>
#include <utility>
#include <string>
#include "bwt.hh"
#include "Interval.hh"
#include "sdsl/io.hpp"
#include "Interfaces.hh"
#include "prezza/rle_string.h"
#include <stdexcept>
#include <iostream>
#include <fstream>


template<class t_bitvector = sdsl::bit_vector>
class RLEBWT : public BWT {
    
public:
    
private:
    
    lzrlbwt::rle_string<> bwt;
    
    std::vector<int64_t> global_c_array;
    std::vector<uint8_t> alphabet;
    
    std::vector<uint8_t> get_string_alphabet(const uint8_t* s) const;
    int64_t strlen(const uint8_t* str) const;
    std::vector<int64_t> compute_c_array(const uint8_t* s, int64_t s_length);
    
public:
    
    // The input string must not contain the END byte
    static const uint8_t END = 0x01; // End of string marker.
    RLEBWT() {}
    //RLEBWT(const uint8_t* input);
    
    virtual void init_from_text(const uint8_t* input);
    virtual void init_from_bwt(const uint8_t* bwt);
    uint8_t get_END() const { return END; }
    int64_t size() const { return bwt.size();}
    //uint8_t bwt_at(int64_t index) const { return bwt[index]; }
    const std::vector<int64_t>& get_global_c_array() const { return global_c_array; }
    const std::vector<uint8_t>& get_alphabet() const { return alphabet; }
    
    
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
    }
    
    void save_to_disk(std::string directory, std::string filename_prefix);
    void load_from_disk(std::string directory, std::string filename_prefix);
    
};

template<class t_bitvector>
const uint8_t RLEBWT<t_bitvector>::END;

template<class t_bitvector>
std::vector<int64_t> RLEBWT<t_bitvector>::compute_c_array(const uint8_t* s, int64_t s_length){
    vector<int64_t> C(256);

    // Put in the counts
    for(int64_t i = 0; i < s_length; i++){
        C[s[i]]++;
    }
    
    // Do a cumulative sum
    for(int64_t i = 1; i < C.size(); i++){
        C[i] += C[i-1];
    }
    
    // Shift to right by one because C[i] is supposed to be the
    // number of characters with lex-rank strictly smaller than i
    for(int64_t i = C.size()-1-1; i >= 0; i--){
        C[i+1] = C[i];
    }
    C[0] = 0;
    return C;
}

// Returns the alphabet in sorted order
template<class t_bitvector>
std::vector<uint8_t> RLEBWT<t_bitvector>::get_string_alphabet(const uint8_t* s) const{
    
    std::vector<bool> found(256,false);
    while(*s != 0){
        found[*s] = true;
        s++;
    }
    
    std::vector<uint8_t> alphabet;
    for(int64_t i = 0; i < 256; i++){
        if(found[i]) alphabet.push_back((uint8_t)i);
    }
    
    return alphabet;
}

// strlen(const uint8_t*) is not in the standard library
template<class t_bitvector>
int64_t RLEBWT<t_bitvector>::strlen(const uint8_t* str) const{
    const uint8_t* start = str;
    while(*str != 0) str++;
    return str - start;
}

template<class t_bitvector>
void RLEBWT<t_bitvector>::init_from_text(const uint8_t* input){

    if(*input == 0) throw std::runtime_error("Tried to construct BD_BWT_index for an empty string");

    global_c_array.resize(256);
    for(int64_t i = 0; i < global_c_array.size(); i++) global_c_array[i] = 0;    
    
    int64_t n = strlen(input);
    
    if(std::find(input, input+n, END) != input + n){
        std::stringstream error;
        error << "Input string contains forbidden byte " << std::hex << END;
        throw std::runtime_error(error.str());
    }
    
    uint8_t* input_copy = (uint8_t*)malloc(sizeof(uint8_t) * (n+1)); // Can't use input directly because it's a const
    for(int64_t i = 0; i < n+1; i++) input_copy[i] = input[i];
    
    uint8_t* data_bwt = build_bwt(input_copy,n,END); // bwt_dbwt reallocs to size n+2, reads first n and appends an END and a null
    std::string cppstring((char*)data_bwt); // For rle_string
    
    // Run length compression
    bwt = lzrlbwt::rle_string<>(cppstring);
    this->alphabet = get_string_alphabet(data_bwt);
    this->global_c_array = compute_c_array(data_bwt, n+1);
    
    free(data_bwt);
    
}

template<class t_bitvector>
void RLEBWT<t_bitvector>::init_from_bwt(const uint8_t* input){
    
    if(*input == 0) throw std::runtime_error("Tried to construct BD_BWT_index for an empty string");

    global_c_array.resize(256);
    for(int64_t i = 0; i < global_c_array.size(); i++) global_c_array[i] = 0;    

    std::string cppstring((char*)input); // For rle_string
    
    // Run length compression
    bwt = lzrlbwt::rle_string<>(cppstring);
    this->alphabet = get_string_alphabet(input);
    this->global_c_array = compute_c_array(input, strlen(input));
    
}


template<class t_bitvector>
void RLEBWT<t_bitvector>::save_to_disk(std::string directory, std::string filename_prefix){
    std::string bwt_path = directory + "/" + filename_prefix + "_bwt.dat";
    ofstream bwt_out(bwt_path);
    bwt.serialize(bwt_out);
    if(!bwt_out.good()){
        cerr << "Error writing to file: " << bwt_path << endl;
        exit(-1);
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
    
    ofstream info(directory + "/" + filename_prefix + "_bwt_info");
    info << "rle_bwt" << endl;
    if(!info.good()){
        cerr << "Error writing to disk: " << directory + "/" + filename_prefix + "_bwt_info" << endl;
        exit(-1);
    }
}

template<class t_bitvector>
void RLEBWT<t_bitvector>::load_from_disk(std::string directory, std::string filename_prefix){
    global_c_array.resize(256);
    
    std::string bwt_path = directory + "/" + filename_prefix + "_bwt.dat";
    ifstream bwt_in(bwt_path);
    if(!bwt_in.good()){
        cerr << "Error opening file: " << bwt_path << endl;
        exit(-1);
    }
    bwt.load(bwt_in);

    
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


#endif
