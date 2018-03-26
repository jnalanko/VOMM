#ifndef CONTEXT_MARKING_HH
#define CONTEXT_MARKING_HH

#include <iostream>
#include "BD_BWT_index/include/BD_BWT_index.hh"
#include "BD_BWT_index/include/Iterators.hh"
#include "sdsl/rank_support_v.hpp"
#include "sdsl/select_support_mcl.hpp"
#include "sdsl/bp_support_g.hpp"
#include "sdsl/util.hpp"
#include <string>
#include <set>
#include <algorithm>
#include <numeric>
#include <vector>
#include "Interfaces.hh"

using namespace std;

class Stats_writer{
// Wrapper for ofstream that can be disabled
private:

    Global_Data* G;
    bool enabled;
    ofstream outfile;    
    
public:
    
    Stats_writer() : G(nullptr), enabled(false) {
        outfile.exceptions(std::ifstream::failbit | std::ifstream::badbit);
    }
    
    void set_file(string filename){
        outfile.open(filename);
        enabled = true;
    }
    
    void set_data(Global_Data* G){
        this->G = G;
    }

    void write_depths(int64_t stringdepth, int64_t open_paren){
        if(enabled){
        assert(G != nullptr);
            outfile << stringdepth << " " << G->rev_st_bpr->excess(open_paren) - 1 << " ";
        }
    }
    
    template<typename T>
    friend Stats_writer& operator<<(Stats_writer& wr, const T& data); 
    
};

template<typename T>
Stats_writer& operator<<(Stats_writer& wr, const T& data){  
    if(!wr.enabled) return wr;
    else{
        wr.outfile << data;
        return wr;
    }
}  

// Formulas to define which strings are contexts

class Entropy_Formula : public Context_Callback{
    
    
private:
    
    // No copying
    Entropy_Formula(Entropy_Formula const& other);
    Entropy_Formula& operator=(Entropy_Formula const& other);
    
public:
    
    
    double threshold;
    BIBWT* index;
    int64_t rev_st_bpr_size;
    Topology_Mapper* mapper;
    sdsl::bit_vector marks;
    int64_t depth_bound;
    int64_t n_candidates;
    Stats_writer* writer;

    // Reusable space
    BD_BWT_index<>::Interval_Data D_W_forward;
    BD_BWT_index<>::Interval_Data D_W_reverse;
    BD_BWT_index<>::Interval_Data D_aW_reverse;
    
    Entropy_Formula(double threshold) : threshold(threshold), depth_bound(1e18), n_candidates(0), writer(nullptr) {}
    Entropy_Formula(double threshold, double depth_bound) : threshold(threshold), depth_bound(depth_bound), n_candidates(0), writer(nullptr) {}
    
    virtual void init(BIBWT* index, int64_t rev_st_bpr_size, Topology_Mapper& mapper, Stats_writer* writer){
        this->index = index;
        this->rev_st_bpr_size = rev_st_bpr_size;
        this->mapper = &mapper;
        this->writer = writer;
        
        marks = sdsl::bit_vector(rev_st_bpr_size, 0);
        
        // Always mark root
        marks[0] = 1;
        marks[marks.size()-1] = 1;
        
        D_W_forward.symbols.resize(index->get_alphabet().size());
        D_W_forward.ranks_start.resize(index->get_alphabet().size());
        D_W_forward.ranks_end.resize(index->get_alphabet().size());
        D_W_reverse.symbols.resize(index->get_alphabet().size());
        D_W_reverse.ranks_start.resize(index->get_alphabet().size());
        D_W_reverse.ranks_end.resize(index->get_alphabet().size());
        D_aW_reverse.symbols.resize(index->get_alphabet().size());
        D_aW_reverse.ranks_start.resize(index->get_alphabet().size());
        D_aW_reverse.ranks_end.resize(index->get_alphabet().size());
    }
    
    virtual void callback(const Iterator::Stack_frame& top){
        Interval_pair I = top.intervals;
        if(!top.is_maxrep || top.depth > depth_bound) return;
        n_candidates++;
        
        index->compute_bwt_interval_data(I.forward, D_W_forward);
        index->compute_rev_bwt_interval_data(I.reverse, D_W_reverse);
        
        int64_t f_W = I.forward.size(); // Number of occurrences of the current string W
        if(D_W_reverse.symbols[0] == index->get_END()){
            // We are not interested if W occurs in the last position of T, because 
            // there is no information about what follows on the right
            f_W--; // f_W is still guaranteed positive, because initially it was >= 2
        }

        double EQ7 = f_W * get_entropy(*index, f_W, D_W_reverse); // Equation 7 on page 10 of the paper
        
        for(int64_t i = 0; i < D_W_forward.n_distinct_symbols; i++){
            Interval_pair I_aW = index->left_extend(I,D_W_forward,i);
            index->compute_rev_bwt_interval_data(I_aW.reverse, D_aW_reverse);
            
            int64_t f_aW = I_aW.forward.size();
                        
            if(D_aW_reverse.symbols[0] == index->get_END()){
                // We are not interested if aW occurs in the last position of T, because 
                // there is no information about what follows on the right
                f_aW--;
            }
            EQ7 -= f_aW * get_entropy(*index, f_aW, D_aW_reverse);
        }
        
        if(EQ7 >= threshold){
            int64_t open = mapper->leaves_to_node(I.reverse);
            int64_t close = mapper->find_close(open);
            marks[open] = 1;
            marks[close] = 1;
            writer->write_depths(top.depth, open);
            *writer << EQ7 << "\n";
        }
    }
    
    double get_entropy(BIBWT& index, int64_t f_W, BIBWT::Interval_Data& D_reverse){
        double ans = 0;
        //index.compute_rev_bwt_interval_data(I_W.reverse, D);
        for(int64_t i = 0; i < D_reverse.n_distinct_symbols; i++){
            char b = D_reverse.symbols[i];
            if(b == index.get_END()) continue;
            int64_t f_Wb = D_reverse.ranks_end[i] - D_reverse.ranks_start[i];
            if(f_Wb != 0){
                ans -= (f_Wb / (double)f_W) * log2(f_Wb / (double)f_W);
            }
        }
        return ans;
    }
    
    virtual void finish(){}
    
    virtual sdsl::bit_vector get_result(){
        return marks;
    }
    
    virtual int64_t get_number_of_candidates(){
        return n_candidates;
    }
};

class EQ234_Formula : public Context_Callback{
    
private:
    
    // No copying
    EQ234_Formula(EQ234_Formula const& other);
    EQ234_Formula& operator=(EQ234_Formula const& other);

public:
    
    double tau1, tau2, tau3, tau4;
    
    // Reusable space
    BD_BWT_index<>::Interval_Data D_W_forward;
    BD_BWT_index<>::Interval_Data D_aW_reverse;
    BD_BWT_index<>::Interval_Data D_W_reverse;
    sdsl::bit_vector marks;
    Topology_Mapper* mapper;
    BIBWT* index;
    int64_t depth_bound;
    int64_t n_candidates;
    Stats_writer* writer;
    
    EQ234_Formula(double tau1, double tau2, double tau3, double tau4) 
    : tau1(tau1), tau2(tau2), tau3(tau3), tau4(tau4), depth_bound(1e18), n_candidates(0), writer(nullptr) {
        assert(tau1 > 0 && tau2 > 0 && tau3 < 1 && tau3 > 0 && tau4 > 1);
    }
    
    EQ234_Formula(double tau1, double tau2, double tau3, double tau4, double depth_bound) 
    : tau1(tau1), tau2(tau2), tau3(tau3), tau4(tau4), depth_bound(depth_bound), n_candidates(0), writer(nullptr) {
        assert(tau1 > 0 && tau2 > 0 && tau3 < 1 && tau3 > 0 && tau4 > 1);
    }
    
    virtual void init(BIBWT* index, int64_t rev_st_bpr_size, Topology_Mapper& mapper, Stats_writer* writer){
        assert(tau1 > 0 && tau2 > 0 && tau3 < 1 && tau3 > 0 && tau4 > 1);
        marks = sdsl::bit_vector(rev_st_bpr_size,0);
        this->mapper = &mapper;
        this->index = index;
        this->writer = writer;
        
        // Always mark root
        marks[0] = 1;
        marks[marks.size()-1] = 1;
        
        D_W_forward.symbols.resize(index->get_alphabet().size());
        D_W_forward.ranks_start.resize(index->get_alphabet().size());
        D_W_forward.ranks_end.resize(index->get_alphabet().size());
        D_W_reverse.symbols.resize(index->get_alphabet().size());
        D_W_reverse.ranks_start.resize(index->get_alphabet().size());
        D_W_reverse.ranks_end.resize(index->get_alphabet().size());
        D_aW_reverse.symbols.resize(index->get_alphabet().size());
        D_aW_reverse.ranks_start.resize(index->get_alphabet().size());
        D_aW_reverse.ranks_end.resize(index->get_alphabet().size());
    }
    
    virtual void callback(const Iterator::Stack_frame& top){
        Interval_pair I_W = top.intervals;
        if(!top.is_maxrep || top.depth > depth_bound-1) return; // Interval must be maximal
        double f_W = I_W.forward.size();
        if(f_W == index->size()) 
            f_W--; // W is the empry string. It should occur |T| times in T. Subtract the end sentinel.
        index->compute_bwt_interval_data(I_W.forward, D_W_forward);
        index->compute_rev_bwt_interval_data(I_W.reverse, D_W_reverse);
        for(int64_t i = 0; i < D_W_forward.n_distinct_symbols; i++){
            char a = D_W_forward.symbols[i];
            if(a == index->get_END()) continue;
            n_candidates++;
            Interval_pair I_aW = index->left_extend(I_W,D_W_forward,i);
            double f_aW = I_aW.forward.size();
            double T_size = index->size() - 1; // Subtract the final sentinel character
            double eq2 = f_aW / (T_size - (top.depth+1) +1);
            if(eq2 < tau1) continue;
            
            index->compute_rev_bwt_interval_data(I_aW.reverse, D_aW_reverse);
            int64_t index_in_D_W_reverse = 0;
            for(int64_t j = 0; j < D_aW_reverse.n_distinct_symbols; j++){
                char b = D_aW_reverse.symbols[j];
                if(b == index->get_END()) continue;
                while(index_in_D_W_reverse < D_W_reverse.n_distinct_symbols-1 && D_W_reverse.symbols[index_in_D_W_reverse] < b)
                    index_in_D_W_reverse++; // Works because symbol lists are sorted
                double f_aWb = D_aW_reverse.count(j);
                if(f_aWb == 0) continue; // eq3 is zero -> can't be a context
                double f_Wb = D_W_reverse.count(index_in_D_W_reverse);
                if(f_Wb == 0) continue;
                double eq3 = f_aWb / f_aW;
                double eq4_numerator = f_aWb / f_aW;
                double eq4_denominator = f_Wb / f_W;
                double eq4 = eq4_numerator / eq4_denominator;
                if(eq3 >= tau2 && (eq4 <= tau3 || eq4 >= tau4)){
                    // Mark aW, or more precisely the child of W with label starting with a in the reverse st
                    int64_t open = mapper->leaves_to_node(I_aW.reverse);
                    int64_t close = mapper->find_close(open);
                    marks[open] = 1;
                    marks[close] = 1;
                    writer->write_depths(top.depth + 1, open);
                    *writer << eq2 << " " << eq3 << " " << eq4 << "\n";
                    break;
                }
            }
        }
    }
    
    virtual void finish(){}
    
    virtual sdsl::bit_vector get_result(){
        return marks;
    }
    
    virtual int64_t get_number_of_candidates(){
        return n_candidates;
    }
    
};

class pnorm_Formula : public Context_Callback{

private:
    
    // No copying
    pnorm_Formula(pnorm_Formula const& other);
    pnorm_Formula& operator=(pnorm_Formula const& other);
    
public:
    
    int64_t p;
    double threshold;
    
    // Reusable space
    BD_BWT_index<>::Interval_Data D_W_forward;
    BD_BWT_index<>::Interval_Data D_W_reverse;
    BD_BWT_index<>::Interval_Data D_aW_reverse;
    
    sdsl::bit_vector marks;
    Topology_Mapper* mapper;
    BIBWT* index;
    int64_t depth_bound;
    int64_t n_candidates;
    Stats_writer* writer;
    
    pnorm_Formula(int64_t p, double threshold) : p(p), threshold(threshold), depth_bound(1e18), n_candidates(0), writer(nullptr)  {}
    pnorm_Formula(int64_t p, double threshold, double depth_bound) : p(p), threshold(threshold), depth_bound(depth_bound), n_candidates(0), writer(nullptr)  {}
    
    virtual void init(BIBWT* index, int64_t rev_st_bpr_size, Topology_Mapper& mapper, Stats_writer* writer){
        assert(threshold >= 0);
        marks = sdsl::bit_vector(rev_st_bpr_size,0);
        this->mapper = &mapper;
        this->index = index;
        this->writer = writer;
        
        // Always mark root
        marks[0] = 1;
        marks[marks.size()-1] = 1;
            
        D_W_forward.symbols.resize(index->get_alphabet().size());
        D_W_forward.ranks_start.resize(index->get_alphabet().size());
        D_W_forward.ranks_end.resize(index->get_alphabet().size());
        D_W_reverse.symbols.resize(index->get_alphabet().size());
        D_W_reverse.ranks_start.resize(index->get_alphabet().size());
        D_W_reverse.ranks_end.resize(index->get_alphabet().size()); 
        D_aW_reverse.symbols.resize(index->get_alphabet().size());
        D_aW_reverse.ranks_start.resize(index->get_alphabet().size());
        D_aW_reverse.ranks_end.resize(index->get_alphabet().size());
    }
    
    virtual void callback(const Iterator::Stack_frame& top){
        Interval_pair I_W = top.intervals;
        if(!top.is_maxrep || top.depth > depth_bound-1) return; // Interval must be maximal
        double f_W = I_W.forward.size();
        if(f_W == index->size()) 
            f_W--; // W is the empry string. It should occur |T| times in T. Subtract the end sentinel.
        index->compute_bwt_interval_data(I_W.forward, D_W_forward);
        index->compute_rev_bwt_interval_data(I_W.reverse, D_W_reverse);
        for(int64_t i = 0; i < D_W_forward.n_distinct_symbols; i++){
            char a = D_W_forward.symbols[i];
            if(a == index->get_END()) continue;
            Interval_pair I_aW = index->left_extend(I_W,D_W_forward,i);
            double f_aW = I_aW.forward.size();
            double p_norm = 0;
            index->compute_rev_bwt_interval_data(I_aW.reverse, D_aW_reverse);
            int64_t index_in_D_aW_reverse = 0;
            for(int64_t j = 0; j < D_W_reverse.n_distinct_symbols; j++){
                char b = D_W_reverse.symbols[j];
                if(b == index->get_END()) continue;
                n_candidates++;
                while(index_in_D_aW_reverse < D_aW_reverse.n_distinct_symbols-1 && D_aW_reverse.symbols[index_in_D_aW_reverse] < b)
                    index_in_D_aW_reverse++; // Works because symbol lists are sorted
                double f_aWb = 0;
                if(D_aW_reverse.symbols[index_in_D_aW_reverse] == b)
                    f_aWb = D_aW_reverse.count(index_in_D_aW_reverse);
                double f_Wb = D_W_reverse.count(j);
                if(f_aW != 0 && f_W != 0)
                    p_norm += pow(fabs(f_aWb / f_aW - f_Wb / f_W),p);
            }
            p_norm = f_aW * pow(p_norm, 1.0/p);
            if(p_norm >= threshold){
                int64_t open = mapper->leaves_to_node(I_aW.reverse);
                int64_t close = mapper->find_close(open);
                marks[open] = 1;
                marks[close] = 1;
                writer->write_depths(top.depth + 1, open);
                *writer << p_norm << "\n";
            }
        }
    }
    
    virtual void finish(){}
    
    virtual sdsl::bit_vector get_result(){
        return marks;
    }
    
    virtual int64_t get_number_of_candidates(){
        return n_candidates;
    }

};

class KL_Formula : public Context_Callback{
    
private:
    
    // No copying
    KL_Formula(KL_Formula const& other);
    KL_Formula& operator=(KL_Formula const& other);
    
public:
    
    double threshold;
    // Reusable space
    BD_BWT_index<>::Interval_Data D_W_forward;
    BD_BWT_index<>::Interval_Data D_aW_reverse;
    BD_BWT_index<>::Interval_Data D_W_reverse;
    
    sdsl::bit_vector marks;
    Topology_Mapper* mapper;
    BIBWT* index;
    int64_t depth_bound;
    int64_t n_candidates;
    Stats_writer* writer;
    
    KL_Formula(double threshold) : threshold(threshold), depth_bound(1e18), n_candidates(0), writer(nullptr) {}
    KL_Formula(double threshold, double depth_bound) : threshold(threshold), depth_bound(depth_bound), n_candidates(0), writer(nullptr) {}
    
    virtual void init(BIBWT* index, int64_t rev_st_bpr_size, Topology_Mapper& mapper, Stats_writer* writer){
        assert(threshold >= 0);
        marks = sdsl::bit_vector(rev_st_bpr_size,0);
        this->mapper = &mapper;
        this->index = index;
        this->writer = writer;

        // Always mark root
        marks[0] = 1;
        marks[marks.size()-1] = 1;

        D_W_forward.symbols.resize(index->get_alphabet().size());
        D_W_forward.ranks_start.resize(index->get_alphabet().size());
        D_W_forward.ranks_end.resize(index->get_alphabet().size());
        D_W_reverse.symbols.resize(index->get_alphabet().size());
        D_W_reverse.ranks_start.resize(index->get_alphabet().size());
        D_W_reverse.ranks_end.resize(index->get_alphabet().size());
        D_aW_reverse.symbols.resize(index->get_alphabet().size());
        D_aW_reverse.ranks_start.resize(index->get_alphabet().size());
        D_aW_reverse.ranks_end.resize(index->get_alphabet().size());
    }
    
    virtual void callback(const Iterator::Stack_frame& top){
        Interval_pair I_W = top.intervals;
        if(!top.is_maxrep || top.depth > depth_bound-1) return; // Interval must be maximal
        double f_W = I_W.forward.size();
        if(f_W == index->size()) 
            f_W--; // W is the empry string. It should occur |T| times in T. Subtract the end sentinel.
        index->compute_bwt_interval_data(I_W.forward, D_W_forward);
        index->compute_rev_bwt_interval_data(I_W.reverse, D_W_reverse);
        for(int64_t i = 0; i < D_W_forward.n_distinct_symbols; i++){
            char a = D_W_forward.symbols[i];
            if(a == index->get_END()) continue;
            n_candidates++;
            Interval_pair I_aW = index->left_extend(I_W,D_W_forward,i);
            double f_aW = I_aW.forward.size();
            if(f_aW == 0) continue;
            double KL_divergence = 0;
            index->compute_rev_bwt_interval_data(I_aW.reverse, D_aW_reverse);
            int64_t index_in_D_W_reverse = 0;
            for(int64_t j = 0; j < D_aW_reverse.n_distinct_symbols; j++){
                char b = D_aW_reverse.symbols[j];
                if(b == index->get_END()) continue;
                while(index_in_D_W_reverse < D_W_reverse.n_distinct_symbols-1 && D_W_reverse.symbols[index_in_D_W_reverse] < b)
                    index_in_D_W_reverse++; // Works because symbol lists are sorted
                double f_aWb = D_aW_reverse.count(j);
                double f_Wb = D_W_reverse.count(index_in_D_W_reverse);
                if(f_aW != 0 && f_W != 0 && f_aWb != 0 && f_Wb != 0)
                    KL_divergence += f_aWb * log((f_aWb / f_aW) / (f_Wb / f_W));
            }
            if(KL_divergence >= threshold){
                int64_t open = mapper->leaves_to_node(I_aW.reverse);
                int64_t close = mapper->find_close(open);
                marks[open] = 1;
                marks[close] = 1;
                writer->write_depths(top.depth + 1, open);
                *writer << KL_divergence << "\n";
            }
        }
    }
    
    virtual void finish(){}
    
    virtual sdsl::bit_vector get_result(){
        return marks;
    }
    
    virtual int64_t get_number_of_candidates(){
        return n_candidates;
    }
};

void iterate_with_callback(Iterator& iterator, Iterator_Callback* cb){
    
    iterator.init();
    while(iterator.next()){
        cb->callback(iterator.get_top());
    }
    cb->finish();
}

template <typename topology_mapper_t>
sdsl::bit_vector mark_contexts_entropy(BIBWT& index, int64_t rev_st_bpr_length, double threshold,
                                       topology_mapper_t& mapper){
    SLT_Iterator iterator(&index);
    Entropy_Formula F(threshold);
    Stats_writer wr;
    F.init(&index, rev_st_bpr_length, mapper, &wr);
    iterate_with_callback(iterator, &F);
    return F.get_result();
    //return mark_contexts_entropy(index,rev_st_bpr_length,threshold,iterator,mapper);
}

// Based on equations 2,3 and 4 of the paper A Framework for Space-efficient String Kernels.
// Marks both opening and close parenthesis
// Let W be a maximal repeat. The contexts are strings of the form aW such that:
// 2) f(aW) / (|T| - |aW| + 1) > \tau_1
// 3) f(aWb) / f(aW) > \tau_2
// 4) (f(aWb) / f(aW)) / (f(Wb) / f(W)) \in (0..\tau_3] \cup [\tau_4..\infty]
// If aW is a context, marks the shortest reverse suffix tree node, that has aW as a prefix
// More detailed comment about the corner cases in the brute force function for testing
// that should do the same thing as this function
template <typename topology_mapper_t>
sdsl::bit_vector mark_contexts_formulas234(BIBWT& index, int64_t rev_st_bpr_length, double tau1,
                                      double tau2, double tau3, double tau4, topology_mapper_t& mapper){
    SLT_Iterator iterator(&index);
    EQ234_Formula F(tau1,tau2,tau3,tau4);
    Stats_writer wr;
    F.init(&index, rev_st_bpr_length, mapper, &wr);
    iterate_with_callback(iterator, &F);
    return F.get_result();
    //return mark_contexts_formulas234(index,rev_st_bpr_length,tau1,tau2,tau3,tau4,iterator,mapper);                                          
}


template <typename topology_mapper_t>
sdsl::bit_vector mark_contexts_KL(BIBWT& index, int64_t rev_st_bpr_length, double threshold, topology_mapper_t& mapper){
    SLT_Iterator iterator(&index);
    KL_Formula F(threshold);
    Stats_writer wr;
    F.init(&index, rev_st_bpr_length, mapper, &wr);
    iterate_with_callback(iterator, &F);
    return F.get_result();
    //return mark_contexts_KL(index,rev_st_bpr_length,threshold,iterator,mapper);     
}
template <typename topology_mapper_t>
sdsl::bit_vector mark_contexts_p_norm(BIBWT& index, int64_t rev_st_bpr_length, double p, double threshold, topology_mapper_t& mapper){
    SLT_Iterator iterator(&index);
    pnorm_Formula F(p,threshold);
    Stats_writer wr;
    F.init(&index, rev_st_bpr_length, mapper, &wr);
    iterate_with_callback(iterator, &F);
    return F.get_result();
    //return mark_contexts_p_norm(index,rev_st_bpr_length,p,threshold,iterator,mapper);     
}



#endif
