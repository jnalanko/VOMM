#ifndef RLE_BITVECTOR
#define RLE_BITVECTOR

#include "Interfaces.hh"
#include "prezza/rle_string.h"
#include <stdexcept>

class RLE_bitvector : public Bitvector{
public:
    
    lzrlbwt::rle_string<> bv;
    
    bool have_bps;
    bool have_ss_10;
    bool have_rs_10;
    bool have_rs;
    bool have_ss;
        
    RLE_bitvector() : have_bps(false), have_ss_10(false), have_rs_10(false), have_rs(false), have_ss(false) {}
    RLE_bitvector(const sdsl::bit_vector& B) : have_bps(false), have_ss_10(false), have_rs_10(false), have_rs(true), have_ss(true) {
        string S;
        for(int i = 0; i < B.size(); i++){
            if(B[i] == 0) S += '0';
            else S += '1';
        }
        this->bv = lzrlbwt::rle_string<>(S);
    }
    
    virtual int64_t size(){
        return bv.size();
    }
    
    virtual bool operator[](int64_t i){
        return bv[i] == '0' ? 0 : 1;
    }
    
    virtual bool at(int64_t i){
        return bv[i] == '0' ? 0 : 1;;
    }
    
    virtual void serialize(string path){
        ofstream outfile(path + "_bv_and_support");
        outfile.exceptions(ifstream::failbit | ifstream::badbit);
        bv.serialize(outfile);

        ofstream info(path + "_info");
        info << have_bps << " " << have_ss_10 << " " << have_rs_10 << " " << have_rs << " " << have_ss << endl;
        if(!info.good()){
            cerr << "Error writing to disk: " << path + "_info" << endl;
            exit(-1);
        }
    }
    
    virtual void load(string path){
        ifstream infile(path + "_bv_and_support");
        infile.exceptions(ifstream::failbit | ifstream::badbit);
        bv.load(infile);

        ifstream info;
        info.exceptions(ifstream::failbit | ifstream::badbit);
        try{
            info.open(path + "_info");
            info >> have_bps >> have_ss_10 >> have_rs_10 >> have_rs >> have_ss;
            info.close();
        }  catch(ifstream::failure e) {
            cerr << "Error loading data structure from disk: " << path + "_info" << endl;
            exit(-1);
        }        
    }
    
    virtual int64_t rank(int64_t pos){
        assert(have_rs);
        return bv.rank(pos, '1');
    }
    
    virtual int64_t rank_10(int64_t pos){
        (void) pos;
        throw(std::runtime_error("Rank 10 not implemented for run length coded bit vector"));
    }
    
    virtual int64_t select(int64_t rank){
        assert(have_ss);
        return bv.select(rank, '1');
    }
    
    virtual int64_t select_10(int64_t pos){
        (void) pos;
        throw(std::runtime_error("Select 10 not implemented for run length coded bit vector"));
    }
    
    // Balanced parentheses operations. Only work after bps initialization
    virtual int64_t find_close(int64_t open){
        (void) open;
        throw(std::runtime_error("BPS not implemented for run length coded bit vector"));
    }
    
    virtual int64_t find_open(int64_t close){
        (void) close;
        throw(std::runtime_error("BPS not implemented for run length coded bit vector"));
    }
    
    virtual int64_t enclose(int64_t open){
        (void) open;
        throw(std::runtime_error("BPS not implemented for run length coded bit vector"));
    }
    
    virtual int64_t double_enclose(int64_t open1, int64_t open2){
        (void) open1; (void) open2;
        throw(std::runtime_error("BPS not implemented for run length coded bit vector"));
    }
    
    virtual int64_t excess(int64_t pos){
        (void) pos;
        throw(std::runtime_error("BPS not implemented for run length coded bit vector"));
    }
    
    virtual void init_rank_support(){
        // Already have
    }
    
    virtual void init_select_support(){
        // Already have
    }
    
    virtual void init_rank_10_support(){
        throw(std::runtime_error("Rank 10 not implemented for run length coded bit vector"));
    }
    
    virtual void init_select_10_support(){
        throw(std::runtime_error("Select 10 not implemented for run length coded bit vector"));
    }
    
    virtual void init_bps_support(){
        throw(std::runtime_error("BPS not implemented for run length coded bit vector"));
    }
};

#endif