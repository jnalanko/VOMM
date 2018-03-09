#ifndef BASIC_BITVECTOR_HH
#define BASIC_BITVECTOR_HH

#include "sdsl/bit_vectors.hpp"
#include "Interfaces.hh"

// Wrapper for sdsl
class Basic_bitvector : public Bitvector{
    
public:
     
    sdsl::bit_vector bv;
    sdsl::bp_support_g<> bps;
    sdsl::select_support_mcl<10,2> ss_10;
    sdsl::rank_support_v<10,2> rs_10;
    sdsl::rank_support_v<1> rs;
    sdsl::select_support_mcl<1> ss;
    
    bool have_bps;
    bool have_ss_10;
    bool have_rs_10;
    bool have_rs;
    bool have_ss;
        
    Basic_bitvector() : have_bps(false), have_ss_10(false), have_rs_10(false), have_rs(false), have_ss(false) {}
    Basic_bitvector(sdsl::bit_vector bv) : bv(bv), have_bps(false), have_ss_10(false), have_rs_10(false), have_rs(false), have_ss(false) {}
    
    virtual int64_t size(){
        return bv.size();
    }
    
    virtual bool operator[](int64_t i){
        return bv[i];
    }
    
    virtual bool at(int64_t i){
        return bv[i];
    }
    
    template<typename T>
    void store_check_error(T& data, string path){
        if(!sdsl::store_to_file(data, path)){
            throw std::runtime_error("Error writing to disk: " + path);        
        }
    }
    
    virtual void serialize(string path){
        store_check_error(bv, path + "_bv");
        
        if(have_bps) store_check_error(bps, path + "_bps");
        if(have_ss_10) store_check_error(ss_10, path + "_ss_10");
        if(have_rs_10) store_check_error(rs_10, path + "_rs_10");
        if(have_rs) store_check_error(rs, path + "_rs");
        if(have_ss) store_check_error(ss, path + "_ss");
        
        ofstream info(path + "_info");
        info << have_bps << " " << have_ss_10 << " " << have_rs_10 << " " << have_rs << " " << have_ss << endl;
        if(!info.good()){
            cerr << "Error writing to disk: " << path + "_info" << endl;
            exit(-1);
        }
    }
    
    
    template <typename D, typename S>
    void load_support_check_error(D& data, S& support, string path){
        ifstream in;
        in.exceptions ( ifstream::failbit | ifstream::badbit);
        try{
            in.open(path);
            support.load(in, &data);
            in.close();
        }  catch(ifstream::failure e) {
            cerr << "Error loading data structure from disk: " << path << endl;
            exit(-1);
        }
    }
    
    virtual void load(string path){
        if(!sdsl::load_from_file(bv, path + "_bv")) {
            throw std::runtime_error("Error reading from disk: " + path + "_bv");
        }

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

        if(have_bps) load_support_check_error(bv, bps, path + "_bps");
        if(have_ss_10) load_support_check_error(bv, ss_10, path + "_ss_10");
        if(have_rs_10) load_support_check_error(bv, rs_10, path + "_rs_10");
        if(have_rs) load_support_check_error(bv, rs, path + "_rs");
        if(have_ss) load_support_check_error(bv, ss, path + "_ss");
        

        
    }
    
    virtual int64_t rank(int64_t pos){
        assert(have_rs);
        return rs.rank(pos);
    }
    
    virtual int64_t rank_10(int64_t pos){
        assert(have_rs_10);
        return rs_10.rank(pos);
    }
    
    virtual int64_t select(int64_t rank){
        assert(have_ss);
        return ss.select(rank);
    }
    
    virtual int64_t select_10(int64_t pos){
        assert(have_ss_10);
        return ss_10.select(pos);
    }
    
    // Balanced parentheses operations. Only work after bps initialization
    virtual int64_t find_close(int64_t open){
        assert(have_bps);
        return bps.find_close(open);
    }
    
    virtual int64_t find_open(int64_t close){
        assert(have_bps);
        return bps.find_open(close);
    }
    
    virtual int64_t enclose(int64_t open){
        assert(have_bps);
        return bps.enclose(open);
    }
    
    virtual int64_t double_enclose(int64_t open1, int64_t open2){
        assert(have_bps);
        return bps.double_enclose(open1, open2);
    }
    
    virtual int64_t excess(int64_t pos){
        assert(have_bps);
        return bps.excess(pos);
    }
    
    virtual void init_rank_support(){
        sdsl::util::init_support(rs, &bv);
        have_rs = true;
    }
    
    virtual void init_select_support(){
        sdsl::util::init_support(ss, &bv);
        have_ss = true;
    }
    
    virtual void init_rank_10_support(){
        sdsl::util::init_support(rs_10, &bv);
        have_rs_10 = true;
    }
    
    virtual void init_select_10_support(){
        sdsl::util::init_support(ss_10, &bv);
        have_ss_10 = true;
    }
    
    virtual void init_bps_support(){
        sdsl::util::init_support(bps, &bv);
        have_bps = true;
    }
    
    virtual std::string toString(){
        stringstream ss;
        ss << bv;
        return ss.str();
    }
};



#endif