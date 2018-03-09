#ifndef COUNTERS_HH
#define COUNTERS_HH
#include "Interfaces.hh"

class Basic_Counters : public Counters {
    
public:
    
    vector<int64_t> v;
    
    virtual void init(int64_t size){
        v.clear();
        v.resize(size,0);
    }
    
    virtual void increment(int64_t pos){
        v[pos]++;
    }
    
    virtual int64_t get(int64_t pos){
        return v[pos];
    }
    
    virtual int64_t size(){
        return v.size();
    }
};

#endif