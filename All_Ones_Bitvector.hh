#ifndef ALL_ONES_BITVECTOR
#define ALL_ONES_BITVECTOR

#include "Interfaces.hh"
#include <stdexcept>
#include <fstream>

class All_Ones_Bitvector : public Bitvector{
    
private:
    
    int64_t length; // Can't be named 'size', because there is a function with the same name
    
public:
    
    All_Ones_Bitvector() : length(0) {}
    All_Ones_Bitvector(int64_t size) : length(size) {}
    
    virtual int64_t size(){
        return length;
    }
    virtual bool operator[](int64_t i){
        (void) i;
        return 1;
    }
    virtual bool at(int64_t i){
        (void) i;
        return 1;
    }
    
    // Rank and select. Only work after initialization
    virtual int64_t rank(int64_t pos){ // Number of ones in [0, pos)
        return pos;
    }
    virtual int64_t rank_10(int64_t pos){
        (void) pos;
        throw(std::runtime_error("Error: All_Ones_Bitvector.rank_10 not implemented."));
        return -1;
    }
    virtual int64_t select(int64_t rank){ // rank \in {1,...,size}
        return rank-1;
    }
    virtual int64_t select_10(int64_t pos){
        (void) pos;
        throw(std::runtime_error("Error: All_Ones_Bitvector.select_10 not implemented."));
        return -1;
    }
    
    virtual int64_t find_close(int64_t open){
        (void) open;
        throw(std::runtime_error("All_Ones_Bitvector.find_close not implemented"));            
    }
    virtual int64_t find_open(int64_t close){
        (void) close;
        throw(std::runtime_error("All_Ones_Bitvector.find_open not implemented"));            
    }
    virtual int64_t enclose(int64_t open){
        (void) open;
        throw(std::runtime_error("All_Ones_Bitvector.enclose not implemented"));
    }
    virtual int64_t double_enclose(int64_t open1, int64_t open2){
        (void) open1; (void) open2;
        throw(std::runtime_error("All_Ones_Bitvector.double_enclose not implemented"));
    }
    virtual int64_t excess(int64_t pos){
        (void) pos;
        throw(std::runtime_error("All_Ones_Bitvector.excess not implemented"));
    }
    
    // Initialization
    virtual void init_rank_support() {}
    virtual void init_select_support() {}
    virtual void init_rank_10_support() {
        throw(std::runtime_error("All_Ones_Bitvector.rank_10 not implemented"));            
    }
    
    virtual void init_select_10_support() {
        throw(std::runtime_error("All_Ones_Bitvector.select_10 not implemented"));            
    }
    
    virtual void init_bps_support(){
        throw(std::runtime_error("All_Ones_Bitvector BPS not implemented"));            
    }
    
    virtual void serialize(std::string path){
        ofstream info(path + "_info");
        info << "all-ones " << this->length << endl;
        if(!info.good()){
            cerr << "Error writing to disk: " << path + "_info" << endl;
            exit(-1);
        }
    }
    virtual void load(std::string path){
        ifstream info;
        info.exceptions ( ifstream::failbit | ifstream::badbit);
        info.open(path + "_info");
        string description;
        info >> description >> this->length;
    }
    
    virtual std::string toString(){
        return "All ones";
    }
    
    virtual ~All_Ones_Bitvector() {} // https://stackoverflow.com/questions/8764353/what-does-has-virtual-method-but-non-virtual-destructor-warning-mean-durin
};

#endif