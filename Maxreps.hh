#ifndef MAXREPS_HH
#define MAXREPS_HH

#include <iostream>
#include "BD_BWT_index/include/BD_BWT_index.hh"
#include "BPR_tools.hh"
#include <string>
#include <vector>
#include <utility>

// This file is for debugging purposes

// Type of a callback function pointer. Template:
// void callback(BD_BWT_index<>* index, std::string* label, Interval_pair current_intervals, int64_t* preorder_rank, void* callback_data);
// The void pointer is data needed by the user
typedef void (*callback_type)(BD_BWT_index<>*, std::string*, Interval_pair, int64_t*, void*);

void maxreps_with_callbacks_recursion(BD_BWT_index<>& index, Interval_pair I, std::string& label, int64_t& preorder_rank, callback_type callback, void* callback_data){
    
    if(!index.is_right_maximal(I)) return;
    preorder_rank++;
    
    if(index.is_left_maximal(I)){
        // Found a maximal repeat
        callback(&index, &label, I, &preorder_rank, callback_data);
    }
    
    // Recurse
    for(char c : index.get_alphabet()){
        Interval_pair I2 = index.left_extend(I,c);
        label.push_back(c);
        maxreps_with_callbacks_recursion(index,I2,label,preorder_rank,callback,callback_data);
        label.pop_back();
    }
}

void maxreps_with_callbacks(BD_BWT_index<>& index, callback_type callback, void* callback_data){
    std::string label = "";
    int64_t preorder_rank = -1;
    maxreps_with_callbacks_recursion(index, Interval_pair(0,index.size()-1,0,index.size()-1), label, preorder_rank, callback, callback_data);
}

void list_all_callback(BD_BWT_index<>* index, std::string* label, Interval_pair I, int64_t* preorder_rank, void* data){
    (void) index; (void) preorder_rank; // Silence unused variable compiler warning
    typedef std::vector<std::pair<std::string, Interval_pair> > result_type;
    result_type& result = *(reinterpret_cast<result_type*>(data));
    std::string label_rev(label->rbegin(), label->rend());
    result.push_back({label_rev, I});
}

std::vector<std::pair<std::string, Interval_pair> > find_maxreps(BD_BWT_index<>& index){
    std::vector<std::pair<std::string, Interval_pair> > result;
    maxreps_with_callbacks(index, list_all_callback, &result);
    return result;
}

struct mark_contexts_entropy_params{
    double threshold; // Threshold in equation 7 on page 10 of the paper
    sdsl::bit_vector* marks; // Must be initialized to zeroes
    
    // Converting lex range into pair of parenthesis
    sdsl::select_support_mcl<10,2>* rev_st_ss_10; // Find leaf
    sdsl::bp_support_g<>* rev_st_bps; // Enclose
};

// Marks the opening parenthesis of every context
void mark_contexts_entropy_callback(BD_BWT_index<>* index, std::string* label, Interval_pair I, int64_t* preorder_rank, void* data){
    (void) label; (void) preorder_rank; // Silence unused variable compiler warning
    auto get_entropy = [index](Interval_pair I_W) -> double {
        double ans = 0;
        for(char b : index->get_alphabet()){
            Interval_pair I_Wb = index->right_extend(I_W,b);
            if(I_Wb.forward.size() != 0){
                double W = I_W.forward.size();
                double Wb = I_Wb.forward.size();
                ans -= (Wb / W) * log2(Wb / W);
            }
        }
        return ans;
    };
    
    double EQ7 = I.forward.size() * get_entropy(I); // Equation 7 on page 10 of the paper
    for(char a : index->get_alphabet()){
        Interval_pair I_aW = index->left_extend(I,a);
        EQ7 -= I_aW.forward.size() * get_entropy(I_aW);
    }
    
    mark_contexts_entropy_params& params = *(reinterpret_cast<mark_contexts_entropy_params*>(data));
    if(EQ7 >= params.threshold){
        Interval bpr_location = enclose_leaves(I.reverse.left, I.reverse.right, *(params.rev_st_ss_10), *(params.rev_st_bps));
        (*params.marks)[bpr_location.left] = 1;
    }
    
}





#endif
