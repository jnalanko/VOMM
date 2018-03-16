#ifndef COUNTERS_HH
#define COUNTERS_HH
#include "Interfaces.hh"

class Basic_Counters : public Counters {
    
public:
    
    typedef uint32_t counter_type;
    
    vector<counter_type> v;
    const int64_t maxvalue = (((int64_t)1) << 32) - 1; // 2^32 - 1
    
    virtual void init(int64_t size){
        v.clear();
        v.resize(size,0);
    }
    
    virtual void increment(int64_t pos){
        assert(v[pos] < maxvalue);
        v[pos]++;
    }
    
    virtual int64_t get(int64_t pos){
        return v[pos];
    }
    
    virtual int64_t size(){
        return v.size();
    }
    
    virtual void free_memory(){
        vector<counter_type>().swap(v);
        // https://stackoverflow.com/questions/10464992/c-delete-vector-objects-free-memory
    }
};

#endif