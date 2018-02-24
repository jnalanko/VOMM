#ifndef LMA_SUPPORT_HH
#define LMA_SUPPORT_HH

#include "sdsl/rank_support_v.hpp"
#include "sdsl/select_support_mcl.hpp"
#include "sdsl/bp_support_g.hpp"
#include "sdsl/util.hpp"
#include "globals.hh"
#include "BD_BWT_index/include/BD_BWT_index.hh"

// Lowest marked ancestor support
class LMA_Support{
public:
    
    sdsl::bit_vector* marks;
    sdsl::rank_support_v<1>* marks_rs;
    sdsl::select_support_mcl<1>* marks_ss;
    sdsl::bit_vector* bpr_marked_only;
    sdsl::bp_support_g<>* bpr_marked_only_bps;
    LMA_Support() {}
    
    // Marks assumes both open and close parentheses are marked
    LMA_Support(sdsl::bit_vector* marks,
                sdsl::rank_support_v<1>* marks_rs, 
                sdsl::select_support_mcl<1>* marks_ss, 
                sdsl::bit_vector* bpr_marked_only,
                sdsl::bp_support_g<>* bpr_marked_only_bps) :
        marks(marks), marks_rs(marks_rs), marks_ss(marks_ss), bpr_marked_only(bpr_marked_only), bpr_marked_only_bps(bpr_marked_only_bps) {}
    
    // Constructor for the lazy and tests
    // Marks assumes both open and close parentheses are marked
    LMA_Support(sdsl::bit_vector& bpr, sdsl::bit_vector& marks) : marks(&marks){
        
        // TODO: These leak memory. Can't just delete in destructor because of the other constructor.
        this->marks_rs = new sdsl::rank_support_v<1>();
        this->marks_ss = new sdsl::select_support_mcl<1>();
        this->bpr_marked_only = new sdsl::bit_vector();
        this->bpr_marked_only_bps = new sdsl::bp_support_g<>();
        
        sdsl::util::init_support(*this->marks_rs, this->marks);
        sdsl::util::init_support(*this->marks_ss, this->marks);
        int64_t nMarked = marks_rs->rank(this->marks->size());
        
        // Build bpr for marked only

        this->bpr_marked_only->resize(nMarked);
        int64_t k = 0;
        for(int64_t i = 0; i < (*this->marks).size(); i++){
            if((*this->marks)[i]){
                (*this->bpr_marked_only)[k] = bpr[i];
                k++;
            }
        }
        sdsl::util::init_support(*(this->bpr_marked_only_bps), this->bpr_marked_only);
    }
    
    // Takes the position of an open parenthesis in the bpr
    // Returns the position of the opening parenthesis of the
    // Nearest marked ancestor. The ancestor can be the node itself.
    // If no such ancestor exists, returns -1
    int64_t LMA(int64_t p){
        //assert((*bpr)[p] == 1); // Only works for open parenthesis
        if((*marks)[p]) return p;
        
        int64_t k = marks_rs->rank(p);
        // The current node is a "virtual" node between positions k-1 and k in the marked only bpr
        
        if(k == 0 || k == bpr_marked_only->size()) return -1;
        
        // Four cases: we are in between (), ((, )) or )(
        int64_t open_in_marked_only;
        if((*bpr_marked_only)[k-1] == 1) // (( or ()
            open_in_marked_only = k-1;
        else if((*bpr_marked_only)[k] == 0) // ))
            open_in_marked_only = bpr_marked_only_bps->find_open(k);        
        else{ // )(
            open_in_marked_only = bpr_marked_only_bps->double_enclose(bpr_marked_only_bps->find_open(k-1),k);
            if(open_in_marked_only == bpr_marked_only_bps->size()){
                // This indicates that the double enclose did not find a common
                // ancestor (see sdsl documentation)
                return -1;
            }
        }
        
        return (*marks_ss)(open_in_marked_only+1); // Map back to the original bpr
    }
    
};


#endif
