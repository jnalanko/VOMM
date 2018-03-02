#ifndef INTERFACES_HH
#define INTERFACES_HH

#include <utility>
#include <string>
#include "BD_BWT_index/include/Interval.hh"
#include <iostream>

// Abstract base classes (play the role of function pointers, but because they are classes, then can have internal state)

class Topology {
public:
    virtual int64_t rev_st_string_depth(int64_t node) = 0;
    virtual int64_t rev_st_parent(int64_t node) = 0;
    virtual int64_t rev_st_lma(int64_t node) = 0;
    
    // Mapping between colex intervals and topology nodes
    virtual int64_t leaves_to_node(Interval I) = 0;    
    virtual Interval node_to_leaves(int64_t node) = 0;
    
    virtual ~Topology() {} // https://stackoverflow.com/questions/8764353/what-does-has-virtual-method-but-non-virtual-destructor-warning-mean-durin
    
};

class Bitvector{
public:
    virtual int64_t size() = 0;
    virtual bool operator[](int64_t i) = 0;
    virtual bool at(int64_t i) = 0; // The same as operator []
    
    // Rank and select. Only work after initialization
    virtual int64_t rank(int64_t pos) = 0;
    virtual int64_t rank_10(int64_t pos) = 0;
    virtual int64_t select(int64_t rank) = 0;
    virtual int64_t select_10(int64_t pos) = 0;
    
    // Balanced parentheses operations. Only work after bps initialization
    virtual int64_t find_close(int64_t open) = 0;
    virtual int64_t find_open(int64_t close) = 0;
    virtual int64_t enclose(int64_t open) = 0;
    virtual int64_t double_enclose(int64_t open1, int64_t open2) = 0;
    virtual int64_t excess(int64_t pos) = 0;
    
    // Initialization
    virtual void init_rank_support() = 0;
    virtual void init_select_support() = 0;
    virtual void init_rank_10_support() = 0;
    virtual void init_select_10_support() = 0;
    virtual void init_bps_support() = 0;
    
    virtual void serialize(std::string path) = 0;
    virtual void load(std::string path) = 0;
    
    virtual ~Bitvector() {} // https://stackoverflow.com/questions/8764353/what-does-has-virtual-method-but-non-virtual-destructor-warning-mean-durin
};

std::ostream& operator<<(std::ostream& os, const Bitvector& dt);

class BWT{
public:
    struct Interval_Data {
        
        sdsl::int_vector_size_type n_distinct_symbols;
        std::vector<uint8_t> symbols; // Vector of length n_distinct_symbols
        
        // ranks_start[i] = number of symbols smaller than the i-th symbol before the start of the interval
        std::vector<uint64_t> ranks_start;
        
        // ranks_start[i] = number of symbols smaller than the i-th symbol up to the end of the interval
        std::vector<uint64_t> ranks_end;
    };
    
    virtual uint8_t get_END() const = 0;
    virtual int64_t size() const = 0;
    virtual const std::vector<int64_t>& get_global_c_array() const = 0;
    virtual const std::vector<uint8_t>& get_alphabet() const = 0;
    //virtual void compute_interval_data(Interval I, Interval_Data& data) = 0;
    virtual Interval search(Interval I, uint8_t c) = 0;
    virtual Interval search_with_precalc(Interval I, uint8_t c, Interval_Data& D) = 0;
    virtual void save_to_disk(std::string directory, std::string filename_prefix) = 0;
    virtual void load_from_disk(std::string directory, std::string filename_prefix) = 0;
    
    virtual ~BWT() {} // https://stackoverflow.com/questions/8764353/what-does-has-virtual-method-but-non-virtual-destructor-warning-mean-durin
    
};

class BIBWT {
public:
    
    struct Interval_Data {
        
        sdsl::int_vector_size_type n_distinct_symbols;
        std::vector<uint8_t> symbols; // Vector of length n_distinct_symbols
        
        // ranks_start[i] = number of symbols smaller than the i-th symbol before the start of the interval
        std::vector<uint64_t> ranks_start;
        
        // ranks_start[i] = number of symbols smaller than the i-th symbol up to the end of the interval
        std::vector<uint64_t> ranks_end;
    };
    
    virtual uint8_t get_END() const = 0;
    virtual int64_t size() const = 0;
    virtual const std::vector<int64_t>& get_global_c_array() const = 0;
    virtual const std::vector<uint8_t>& get_alphabet() const = 0;
    virtual void compute_local_c_array_forward(Interval& interval, std::vector<int64_t>& c_array) = 0;
    virtual void compute_local_c_array_reverse(Interval& interval, std::vector<int64_t>& c_array) = 0;
    virtual Interval_pair left_extend(Interval_pair intervals, uint8_t c) = 0;
    virtual Interval_pair left_extend(Interval_pair intervals, uint8_t c, const std::vector<int64_t>& local_c_array) = 0;
    virtual Interval_pair right_extend(Interval_pair intervals, uint8_t c) = 0;
    virtual Interval_pair right_extend(Interval_pair intervals, uint8_t c, const std::vector<int64_t>& local_c_array) = 0;
    virtual int64_t backward_step(int64_t lex_rank) const = 0;
    virtual int64_t forward_step(int64_t colex_rank) const = 0;
    virtual bool is_right_maximal(Interval_pair I) = 0;
    virtual bool is_left_maximal(Interval_pair I) = 0;
    virtual void compute_bwt_interval_data(Interval I, Interval_Data& data) = 0;
    virtual void compute_rev_bwt_interval_data(Interval I, Interval_Data& data) = 0;
    virtual Interval_pair left_extend(Interval_pair intervals, Interval_Data& data, int64_t symbol_index) = 0;
    //Interval_pair right_extend(Interval_pair I, const Interval_Data& data, uint8_t c); // TODO;
    virtual void save_to_disk_reverse_only(std::string directory, std::string filename_prefix) = 0;
    //virtual void load_from_disk(std::string directory, std::string filename_prefix) = 0;
    
    virtual ~BIBWT() {} // https://stackoverflow.com/questions/8764353/what-does-has-virtual-method-but-non-virtual-destructor-warning-mean-durin

};

class Scoring_Function{
public:  
    virtual double score(Interval I, int64_t d, char c, Topology& topology, BWT& index) = 0;
    
    virtual ~Scoring_Function() {} // https://stackoverflow.com/questions/8764353/what-does-has-virtual-method-but-non-virtual-destructor-warning-mean-durin

};

class Global_Data;
class Loop_Invariant_Updater{
public:
    virtual std::pair<Interval, int64_t> update(Interval I, int64_t d, char c, Global_Data& data, Topology& topology, BWT& index) = 0;
    
    virtual ~Loop_Invariant_Updater() {} // https://stackoverflow.com/questions/8764353/what-does-has-virtual-method-but-non-virtual-destructor-warning-mean-durin

};


class Topology_Mapper{
    
public:
    
    typedef int64_t node_t;
        
    virtual node_t leaves_to_node(Interval leaves) = 0;
    virtual Interval node_to_leaves(node_t node) = 0;
    virtual node_t find_close(node_t open) = 0;
    
    virtual ~Topology_Mapper() {} // https://stackoverflow.com/questions/8764353/what-does-has-virtual-method-but-non-virtual-destructor-warning-mean-durin

    
};

class Iterator{
    
public:
    
    class Stack_frame{
    public:
        Interval_pair intervals;  // forward interval, reverse interval
        int64_t depth; // depth in the tree
        Stack_frame(Interval_pair intervals, int64_t depth) : intervals(intervals), depth(depth) {}
        Stack_frame(){}
    };
    
    virtual void init() = 0;
    virtual bool next() = 0;
    virtual void set_index(BIBWT* index) = 0;
    virtual Stack_frame get_top() = 0;
    
    virtual ~Iterator() {} // https://stackoverflow.com/questions/8764353/what-does-has-virtual-method-but-non-virtual-destructor-warning-mean-durin
};

class Context_Formula{
    
public:
    virtual sdsl::bit_vector get_rev_st_context_marks(BIBWT* index, int64_t rev_st_bpr_size, Iterator& candidate_iterator, Topology_Mapper& mapper) = 0;

    virtual ~Context_Formula() {} // https://stackoverflow.com/questions/8764353/what-does-has-virtual-method-but-non-virtual-destructor-warning-mean-durin

};




#endif
