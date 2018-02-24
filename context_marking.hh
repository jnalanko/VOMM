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

double get_entropy(BWT& index, Interval_pair I_W, int64_t f_W, BWT::Interval_Data& D){
    double ans = 0;
    index.compute_rev_bwt_interval_data(I_W.reverse, D);
    for(int64_t i = 0; i < D.n_distinct_symbols; i++){
        char b = D.symbols[i];
        if(b == index.get_END()) continue;
        int64_t f_Wb = D.ranks_end[i] - D.ranks_start[i];
        if(f_Wb != 0){
            ans -= (f_Wb / (double)f_W) * log2(f_Wb / (double)f_W);
        }
    }
    return ans;
}

// Marks both opening and close parenthesis
// Based on equations 7 of the paper A Framework for Space-efficient String Kernels.
template <typename topology_mapper_t>
sdsl::bit_vector mark_contexts_entropy(BWT& index, int64_t rev_st_bpr_length, double threshold,
                                       Iterator& it, topology_mapper_t& mapper){
    
    sdsl::bit_vector marks(rev_st_bpr_length,0);
    // Always mark root
    marks[0] = 1;
    marks[marks.size()-1] = 1;
    
    // Initialize reusable space
    vector<int64_t> local_c_array_forward(256);
    vector<int64_t> local_c_array_reverse(256);
    BD_BWT_index<>::Interval_Data D; // Given to get_entropy
    BD_BWT_index<>::Interval_Data D2; // Used for left-extend
    D.symbols.resize(index.get_alphabet().size());
    D.ranks_start.resize(index.get_alphabet().size());
    D.ranks_end.resize(index.get_alphabet().size());
    D2.symbols.resize(index.get_alphabet().size());
    D2.ranks_start.resize(index.get_alphabet().size());
    D2.ranks_end.resize(index.get_alphabet().size());
    
    
    it.init();
    while(it.next()){
        Interval_pair I = it.get_top().intervals;
        index.compute_local_c_array_reverse(I.reverse, local_c_array_reverse);        
        index.compute_local_c_array_forward(I.forward, local_c_array_forward);   
        
        int64_t f_W = I.forward.size(); // Number of occurrences of the current string W
        if(index.right_extend(I,index.get_END(),local_c_array_reverse).forward.size() > 0){
            // We are not interested if W occurs in the last position of T, because 
            // there is no information about what follows on the right
            f_W--; // f_W is still guaranteed positive, because initially it was >= 2
        }

        double EQ7 = f_W * get_entropy(index, I, f_W, D); // Equation 7 on page 10 of the paper
        index.compute_bwt_interval_data(I.forward, D2);
        
        for(int64_t i = 0; i < D2.n_distinct_symbols; i++){
            char a = D2.symbols[i];
            Interval_pair I_aW = index.left_extend(I,a,local_c_array_forward);
            
            int64_t f_aW = I_aW.forward.size();
                        
            if(index.right_extend(I_aW,index.get_END()).forward.size() > 0){
                // We are not interested if aW occurs in the last position of T, because 
                // there is no information about what follows on the right
                f_aW--;
            }
            EQ7 -= f_aW * get_entropy(index, I_aW, f_aW, D);
        }
        
        if(EQ7 >= threshold){
            int64_t open = mapper.leaves_to_node(I.reverse);
            int64_t close = mapper.find_close(open);
            marks[open] = 1;
            marks[close] = 1;
            //Interval bpr_location = enclose_leaves(I.reverse.left, I.reverse.right, rev_st_ss_10, rev_st_bps);
            //marks[bpr_location.left] = 1;
            //marks[bpr_location.right] = 1;
        }
    }
    
    return marks;
}

template <typename topology_mapper_t>
sdsl::bit_vector mark_contexts_entropy(BWT& index, int64_t rev_st_bpr_length, double threshold,
                                       topology_mapper_t& mapper){
    SLT_Iterator iterator(&index);
    return mark_contexts_entropy(index,rev_st_bpr_length,threshold,iterator,mapper);
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
sdsl::bit_vector mark_contexts_formulas234(BWT& index, int64_t rev_st_bpr_length, double tau1,
                                      double tau2, double tau3, double tau4, Iterator& it, topology_mapper_t& mapper){
    
    assert(tau1 > 0 && tau2 > 0 && tau3 < 1 && tau3 > 0 && tau4 > 1);
    sdsl::bit_vector marks(rev_st_bpr_length,0);
    
    // Always mark root
    marks[0] = 1;
    marks[marks.size()-1] = 1;
    
    vector<int64_t> local_c_array_forward(256);
    // vector<int64_t> local_c_array_reverse(256); // Could optimize this function by using this for Wb and aWb
    it.init();
    while(it.next()){
        Interval_pair I_W = it.get_top().intervals;
        if(!index.is_left_maximal(I_W))
            continue; // Interval must be maximal
        double f_W = I_W.forward.size();
        if(f_W == index.size()) 
            f_W--; // W is the empry string. It should occur |T| times in T. Subtract the end sentinel.
        index.compute_local_c_array_forward(I_W.forward, local_c_array_forward);
        for(char a : index.get_alphabet()){
            if(a == index.get_END()) continue;
            Interval_pair I_aW = index.left_extend(I_W,a,local_c_array_forward);
            double f_aW = I_aW.forward.size();
            double T_size = index.size() - 1; // Subtract the final sentinel character
            
            if(f_aW / (T_size - (it.get_top().depth+1) +1) < tau1) continue;
            for(char b : index.get_alphabet()){
                if(b == index.get_END()) continue;
                Interval_pair I_aWb = index.right_extend(I_aW,b);
                Interval_pair I_Wb = index.right_extend(I_W,b);
                double f_aWb = I_aWb.forward.size();
                double f_Wb = I_Wb.forward.size();
                if(f_Wb == 0) continue;
                double eq3 = f_aWb / f_aW ;
                double eq4_numerator = f_aWb / f_aW;
                double eq4_denominator = f_Wb / f_W;
                double eq4 = eq4_numerator / eq4_denominator;
                if(eq3 >= tau2 && (eq4 <= tau3 || eq4 >= tau4)){
                    // Mark aW, or more precisely the child of W with label starting with a in the reverse st
                    int64_t open = mapper.leaves_to_node(I_aW.reverse);
                    int64_t close = mapper.find_close(open);
                    marks[open] = 1;
                    marks[close] = 1;
                    break;
                }
            }
        }
    }
    
    return marks;
}

template <typename topology_mapper_t>
sdsl::bit_vector mark_contexts_formulas234(BWT& index, int64_t rev_st_bpr_length, double tau1,
                                      double tau2, double tau3, double tau4, topology_mapper_t& mapper){
    SLT_Iterator iterator(&index);
    return mark_contexts_formulas234(index,rev_st_bpr_length,tau1,tau2,tau3,tau4,iterator,mapper);                                          
}

template <typename topology_mapper_t>
sdsl::bit_vector mark_contexts_KL(BWT& index, int64_t rev_st_bpr_length, double threshold, Iterator& it, topology_mapper_t& mapper){
    
    assert(threshold >= 0);
    sdsl::bit_vector marks(rev_st_bpr_length,0);
    
    // Always mark root
    marks[0] = 1;
    marks[marks.size()-1] = 1;
    
    vector<int64_t> local_c_array_forward(256);
    // vector<int64_t> local_c_array_reverse(256); // Could optimize this function by using this for Wb and aWb
    it.init();
    while(it.next()){
        Interval_pair I_W = it.get_top().intervals;
        if(!index.is_left_maximal(I_W))
            continue; // Interval must be maximal
        double f_W = I_W.forward.size();
        if(f_W == index.size()) 
            f_W--; // W is the empry string. It should occur |T| times in T. Subtract the end sentinel.
        index.compute_local_c_array_forward(I_W.forward, local_c_array_forward);
        for(char a : index.get_alphabet()){
            if(a == index.get_END()) continue;
            Interval_pair I_aW = index.left_extend(I_W,a,local_c_array_forward);
            double f_aW = I_aW.forward.size();
            double KL_divergence = 0;
            for(char b : index.get_alphabet()){
                if(b == index.get_END()) continue;
                Interval_pair I_aWb = index.right_extend(I_aW,b);
                Interval_pair I_Wb = index.right_extend(I_W,b);
                double f_aWb = I_aWb.forward.size();
                double f_Wb = I_Wb.forward.size();
                if(f_aW != 0 && f_W != 0 && f_aWb != 0 && f_Wb != 0)
                    KL_divergence += f_aWb * log((f_aWb / f_aW) / (f_Wb / f_W));
            }
            if(KL_divergence >= threshold){
                int64_t open = mapper.leaves_to_node(I_aW.reverse);
                int64_t close = mapper.find_close(open);
                marks[open] = 1;
                marks[close] = 1;
            }
        }
    }
        
    return marks;
}

template <typename topology_mapper_t>
sdsl::bit_vector mark_contexts_KL(BWT& index, int64_t rev_st_bpr_length, double threshold, topology_mapper_t& mapper){
    SLT_Iterator iterator(&index);
    return mark_contexts_KL(index,rev_st_bpr_length,threshold,iterator,mapper);     
}

template <typename topology_mapper_t>
sdsl::bit_vector mark_contexts_p_norm(BWT& index, int64_t rev_st_bpr_length, double p, double threshold, Iterator& it, topology_mapper_t& mapper){
    
    assert(threshold >= 0);
    sdsl::bit_vector marks(rev_st_bpr_length,0);
    
    // Always mark root
    marks[0] = 1;
    marks[marks.size()-1] = 1;
    
    vector<int64_t> local_c_array_forward(256);
    // vector<int64_t> local_c_array_reverse(256); // Could optimize this function by using this for Wb and aWb
    it.init();
    while(it.next()){
        Interval_pair I_W = it.get_top().intervals;
        if(!index.is_left_maximal(I_W))
            continue; // Interval must be maximal
        double f_W = I_W.forward.size();
        if(f_W == index.size()) 
            f_W--; // W is the empry string. It should occur |T| times in T. Subtract the end sentinel.
        index.compute_local_c_array_forward(I_W.forward, local_c_array_forward);
        for(char a : index.get_alphabet()){
            if(a == index.get_END()) continue;
            Interval_pair I_aW = index.left_extend(I_W,a,local_c_array_forward);
            double f_aW = I_aW.forward.size();
            double p_norm = 0;
            for(char b : index.get_alphabet()){
                if(b == index.get_END()) continue;
                Interval_pair I_aWb = index.right_extend(I_aW,b);
                Interval_pair I_Wb = index.right_extend(I_W,b);
                double f_aWb = I_aWb.forward.size();
                double f_Wb = I_Wb.forward.size();
                if(f_aW != 0 && f_W != 0)
                    p_norm += pow(fabs(f_aWb / f_aW - f_Wb / f_W),p);
            }
            p_norm = f_aW * pow(p_norm, 1.0/p);
            if(p_norm >= threshold){
                int64_t open = mapper.leaves_to_node(I_aW.reverse);
                int64_t close = mapper.find_close(open);
                marks[open] = 1;
                marks[close] = 1;
            }
        }
    }
    
    return marks;
}

template <typename topology_mapper_t>
sdsl::bit_vector mark_contexts_p_norm(BWT& index, int64_t rev_st_bpr_length, double p, double threshold, topology_mapper_t& mapper){
    SLT_Iterator iterator(&index);
    return mark_contexts_p_norm(index,rev_st_bpr_length,p,threshold,iterator,mapper);     
}

// Formulas to define which strings are contexts

class Entropy_Formula : public Context_Formula{
    
public:
    
    
    double threshold;
    
    Entropy_Formula(double threshold) : threshold(threshold) {}
    
    virtual sdsl::bit_vector get_rev_st_context_marks(BWT* index, int64_t rev_st_bpr_size, Iterator& candidate_iterator, Topology_Mapper& mapper){
        candidate_iterator.set_index(index);
        return mark_contexts_entropy(*index, rev_st_bpr_size, threshold, candidate_iterator, mapper);
    }
};

class EQ234_Formula : public Context_Formula{
    
public:
    
    double tau1, tau2, tau3, tau4;
    
    EQ234_Formula(double tau1, double tau2, double tau3, double tau4) : tau1(tau1), tau2(tau2), tau3(tau3), tau4(tau4) {}
    
    sdsl::bit_vector get_rev_st_context_marks(BWT* index, int64_t rev_st_bpr_size, Iterator& candidate_iterator, Topology_Mapper& mapper){
        candidate_iterator.set_index(index);
        return mark_contexts_formulas234(*index, rev_st_bpr_size, tau1, tau2, tau3, tau4, candidate_iterator, mapper);
    }
};

class KL_Formula : public Context_Formula{
    
public:
    
    double threshold;
    
    KL_Formula(double threshold) : threshold(threshold) {}
    
    sdsl::bit_vector get_rev_st_context_marks(BWT* index, int64_t rev_st_bpr_size, Iterator& candidate_iterator, Topology_Mapper& mapper){
        candidate_iterator.set_index(index);
        return mark_contexts_KL(*index, rev_st_bpr_size, threshold, candidate_iterator, mapper);
    }
};

class pnorm_Formula : public Context_Formula{
    
public:
    
    int64_t p;
    double threshold;
    
    pnorm_Formula(int64_t p, double threshold) : p(p), threshold(threshold){}
    
    sdsl::bit_vector get_rev_st_context_marks(BWT* index, int64_t rev_st_bpr_size, Iterator& candidate_iterator, Topology_Mapper& mapper){
        candidate_iterator.set_index(index);
        return mark_contexts_p_norm(*index, rev_st_bpr_size, p, threshold, candidate_iterator, mapper);
    }
};

#endif
