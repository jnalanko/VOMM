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
    virtual int64_t rank(int64_t pos) = 0; // Number of ones in [0, pos)
    virtual int64_t rank_10(int64_t pos) = 0;
    virtual int64_t select(int64_t rank) = 0; // rank \in {1,...,size}
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
    
    // Serialize: Writes two files, one to path and one to path + "_info"
    // The info-file contains on line: first an ascii string describing the
    // type of the vector, then a space, and then information related to that
    // particular type.
    virtual void serialize(std::string path) = 0;
    
    // Loads the bit vector at 'path' based on the info file written by serialize.
    virtual void load(std::string path) = 0;
    
    virtual std::string toString() = 0;
    
    virtual ~Bitvector() {} // https://stackoverflow.com/questions/8764353/what-does-has-virtual-method-but-non-virtual-destructor-warning-mean-durin
};


class BWT{
public:
    struct Interval_Data {
        
        sdsl::int_vector_size_type n_distinct_symbols;
        std::vector<uint8_t> symbols; // Vector of length n_distinct_symbols
        
        // ranks_start[i] = number of symbols smaller than the i-th symbol before the start of the interval
        std::vector<uint64_t> ranks_start;
        
        // ranks_start[i] = number of symbols smaller than the i-th symbol up to the end of the interval
        std::vector<uint64_t> ranks_end;
        
        int64_t count(int64_t symbol_index){
            return ranks_end[symbol_index] - ranks_start[symbol_index];
        }
        
    };
    
    
    virtual void init_from_text(const uint8_t* input) = 0;
    virtual void init_from_bwt(const uint8_t* bwt) = 0;
    virtual uint8_t get_END() const = 0;
    virtual int64_t size() const = 0;
    virtual const std::vector<int64_t>& get_global_c_array() const = 0;
    virtual const std::vector<uint8_t>& get_alphabet() const = 0;
    //virtual void compute_interval_data(Interval I, Interval_Data& data) = 0;
    virtual Interval search(Interval I, uint8_t c) = 0;
    virtual Interval search_with_precalc(Interval I, uint8_t c, Interval_Data& D) = 0;
    
    // save_to_disk: also write type information to directory + "/" + filename_prefix + "_bwt_info"
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
        
        int64_t count(int64_t symbol_index){
            return ranks_end[symbol_index] - ranks_start[symbol_index];
        }
        
    };
    
    virtual uint8_t get_END() const = 0;
    virtual int64_t size() const = 0;
    virtual const std::vector<int64_t>& get_global_c_array() const = 0;
    virtual const std::vector<uint8_t>& get_alphabet() const = 0;
    virtual bool is_right_maximal(Interval_pair I) = 0;
    virtual bool is_left_maximal(Interval_pair I) = 0;
    virtual void compute_bwt_interval_data(Interval I, Interval_Data& data) = 0; // data.symbols MUST BE lex-ordered!
    virtual void compute_rev_bwt_interval_data(Interval I, Interval_Data& data) = 0; // data.symbols MUST BE lex-ordered!
    virtual Interval_pair left_extend(Interval_pair intervals, Interval_Data& data, int64_t symbol_index) = 0;
    virtual Interval_pair right_extend(Interval_pair intervals, Interval_Data& data, int64_t symbol_index) = 0;
    virtual void save_to_disk(std::string directory, std::string filename_prefix) = 0;
    virtual void load_from_disk(std::string directory, std::string filename_prefix) = 0;
    
    virtual uint8_t forward_bwt_at(int64_t index) const = 0;
    virtual uint8_t backward_bwt_at(int64_t index) const = 0;
    
    //virtual void save_to_disk_reverse_only(std::string directory, std::string filename_prefix) = 0;
    
    
    
    virtual ~BIBWT() {} // https://stackoverflow.com/questions/8764353/what-does-has-virtual-method-but-non-virtual-destructor-warning-mean-durin

};

class Global_Data;
class Loop_Invariant_Updater{
public:
      /**
       * Updates the matching statistics loop invariant.
       * @param I The colex-interval of the longest match at the current position
       * @param node The deepest node that contains all leaves in interval I. Represented by
       *             the position of the open parenthesis in the BPR of the revese suffix tree.
       *             If the BPR is complete, then this the node with colex-interval I. If the BPR
       *             has been pruned, this is the lowest non-pruned ancestor of the node with
       *             colex-interval I
       * @param d The depth of the longest match at the current position
       * @param c The next character in the text
       * @param data The global data structures.
       * @param topology A class implementing the required topology operations
       * @param index The BWT of the reverse
       * @return Suppose the current match is W. Then we return the colex-interval of the longest
       *         suffix of Wc that is found in the index, and the string length of that.
       */
    virtual std::pair<Interval, int64_t> update(Interval I, int64_t node, int64_t d, char c, Global_Data& data, Topology& topology, BWT& index) = 0;
    
    virtual ~Loop_Invariant_Updater() {} // https://stackoverflow.com/questions/8764353/what-does-has-virtual-method-but-non-virtual-destructor-warning-mean-durin

};

class Scoring_Function{
public:
      /**
       * Computes the log-probability of the current position
       * @param I The colex-interval of the longest match at the current position
       * @param node The deepest node that contains all leaves in interval I. Represented by
       *             the position of the open parenthesis in the BPR of the revese suffix tree.
       *             If the BPR is complete, then this the node with colex-interval I. If the BPR
       *             has been pruned, this is the lowest non-pruned ancestor of the node with
       *             colex-interval I
       * @param d The depth of the longest match at the current position
       * @param c The character for which we want to compute a probability
       * @param topology A class implementing the required topology operations
       * @param index The BWT of the reverse
       * @param data The global data structures.
       * @return Suppose the current longest match is W. Then returns the log-probability of Wc.
       */
    virtual double score(/*Interval I,*/int64_t node, int64_t d, char c, Topology& topology, BWT& index, Global_Data& G) = 0;
    
    virtual ~Scoring_Function() {} // https://stackoverflow.com/questions/8764353/what-does-has-virtual-method-but-non-virtual-destructor-warning-mean-durin

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
        bool is_maxrep;
        Stack_frame(Interval_pair intervals, int64_t depth, bool is_maxrep) : intervals(intervals), depth(depth), is_maxrep(is_maxrep) {}
        Stack_frame(){}
    };
    
    virtual void init() = 0;
    virtual bool next() = 0;
    virtual void set_index(BIBWT* index) = 0;
    virtual Stack_frame get_top() = 0;
    
    virtual ~Iterator() {} // https://stackoverflow.com/questions/8764353/what-does-has-virtual-method-but-non-virtual-destructor-warning-mean-durin
};

class Iterator_Callback{
public:
    virtual void callback(const Iterator::Stack_frame& top) = 0;
    virtual void finish() = 0;
    virtual ~Iterator_Callback() {} // https://stackoverflow.com/questions/8764353/what-does-has-virtual-method-but-non-virtual-destructor-warning-mean-durin

};

class Stats_writer;
class Context_Callback : public Iterator_Callback{
public:
    
    virtual void init(BIBWT* index, int64_t rev_st_bpr_size, Topology_Mapper& mapper, Stats_writer* writer) = 0;
    virtual sdsl::bit_vector get_result() = 0;
    virtual int64_t get_number_of_candidates() = 0;
    virtual ~Context_Callback() {} // https://stackoverflow.com/questions/8764353/what-does-has-virtual-method-but-non-virtual-destructor-warning-mean-durin

};

class Counters{
public:  
    virtual void init(int64_t size) = 0;
    virtual void increment_open(int64_t pos) = 0;
    virtual void increment_close(int64_t pos) = 0;
    virtual int64_t get_open(int64_t pos) = 0;
    virtual int64_t get_close(int64_t pos) = 0;
    virtual int64_t size() = 0;
    virtual void free_memory() = 0;
    
};

class String_Depth_Support{
public:
    virtual int64_t string_depth(int64_t open) = 0;
};

/*
class Context_Formula{
    
public:
    virtual sdsl::bit_vector get_rev_st_context_marks(BIBWT* index, int64_t rev_st_bpr_size, Iterator& candidate_iterator, Topology_Mapper& mapper) = 0;
    virtual ~Context_Formula() {} // https://stackoverflow.com/questions/8764353/what-does-has-virtual-method-but-non-virtual-destructor-warning-mean-durin
}; */




#endif
