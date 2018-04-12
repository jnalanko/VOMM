#ifndef COUNTERS_HH
#define COUNTERS_HH
#include "Interfaces.hh"
#include "logging.hh"

class Basic_Counters : public Counters {
    
public:
    
    typedef uint32_t counter_type;
    
    vector<counter_type> v_open;
    vector<counter_type> v_close;
    const int64_t maxvalue = (((int64_t)1) << 32) - 1; // 2^32 - 1
    
    virtual void init(int64_t size){
        v_open.clear();
        v_open.resize(size,0);
        v_close.clear();
        v_close.resize(size,0);
    }
    
    virtual void increment_open(int64_t pos){
        assert(v_open[pos] < maxvalue);
        v_open[pos]++;
    }
    
    virtual void increment_close(int64_t pos){
        assert(v_close[pos] < maxvalue);
        v_close[pos]++;
    }
    
    virtual int64_t get_open(int64_t pos){
        return v_open[pos];
    }
   
    virtual int64_t get_close(int64_t pos){
        return v_close[pos];
    }
    
    virtual int64_t size(){
        return v_open.size();
    }
    
    virtual void free_memory(){
        vector<counter_type>().swap(v_open);
        vector<counter_type>().swap(v_close);
        // https://stackoverflow.com/questions/10464992/c-delete-vector-objects-free-memory
    }
};

class Succinct_Counters : public Counters {
    
public:
    
    static const int64_t LEVEL_1_SIZE = 8; // Number of bits of each counter on the first level
    static const int64_t LEVEL_2_SIZE = 64; // Number of bits of each counter on the second level
    static const int64_t SATURATED = (1 << LEVEL_1_SIZE)-1; // Value indicating overflow to the second level
    
    int64_t overflows = 0;
    
    sdsl::int_vector<LEVEL_1_SIZE> v1; // Level 1
    unordered_map<int64_t, int64_t> v2; // Level 2  
    
    virtual void init(int64_t size){
        v1.resize(size);
        for(int64_t i = 0; i < v1.size(); i++) v1[i] = 0;
        v2.clear();
    }    
    
    // (0,x) = 0 + 4x for x >= 0
    // (1,x) = 1 + 4x for x >= 0;
    // (x,0) = 2 + (4(x-2)) for x >= 2
    // (x,1) = 3 + (4(x-2)) for x >= 2
    
    int64_t encode(int64_t open, int64_t close){
        if(open == 0) return 0 + 4*close;
        if(open == 1) return 1 + 4*close;
        if(close == 0) return 2 + 4*(open-2);
        if(close == 1) return 3 + 4*(open-2);
        
        assert(false);
        return -1; // Should never happen
    }
    
    pair<int64_t, int64_t> decode(int64_t x){
        if(x % 4 == 0) return {0, x/4};
        if(x % 4 == 1) return {1, (x-1)/4};
        if(x % 4 == 2) return {(x+6)/4,0};
        if(x % 4 == 3) return {(x+5)/4,1};
        assert(false); // Should never happen
    }
    
    virtual void increment_open(int64_t pos){
        pair<int64_t, int64_t> openclose = get_openclose(pos);
        int64_t open = openclose.first;
        int64_t close = openclose.second;
        if(v1[pos] == SATURATED){
            v2[pos] = encode(open+1,close);
        } else{
            int64_t new_codeword = encode(open+1,close);
            if(new_codeword >= SATURATED){
                overflows++;
                // Overflow to level 2
                v1[pos] = SATURATED;
                v2[pos] = encode(open+1,close);
            } else v1[pos] = new_codeword;
        }
    }
    
    virtual void increment_close(int64_t pos){
        pair<int64_t, int64_t> openclose = get_openclose(pos);
        int64_t open = openclose.first;
        int64_t close = openclose.second;
        if(v1[pos] == SATURATED){
            v2[pos] = encode(open,close+1);
        } else{
            int64_t new_codeword = encode(open,close+1);
            if(new_codeword >= SATURATED){
                overflows++;
                // Overflow to level 2
                v1[pos] = SATURATED;
                v2[pos] = encode(open,close+1);
            } else v1[pos] = new_codeword;
        }
    }
    
    pair<int64_t,int64_t> get_openclose(int64_t pos){
        if(v1[pos] != SATURATED) return decode(v1[pos]);
        else return decode(v2[pos]);
    }
    
    virtual int64_t get_open(int64_t pos){
        if(v1[pos] != SATURATED) return decode(v1[pos]).first;
        else return decode(v2[pos]).first;
    }
    
    virtual int64_t get_close(int64_t pos){
        if(v1[pos] != SATURATED) return decode(v1[pos]).second;
        else return decode(v2[pos]).second;
    }
    
    virtual int64_t size(){
        return v1.size();
    }
    
    virtual void free_memory(){
        write_log("Saturated counters: " + to_string(overflows));
        sdsl::int_vector<LEVEL_1_SIZE>().swap(v1);
        std::unordered_map<int64_t, int64_t>().swap(v2);
        // https://stackoverflow.com/questions/10464992/c-delete-vector-objects-free-memory
    }
};

#endif