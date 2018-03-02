#ifndef LMA_SUPPORT_HH
#define LMA_SUPPORT_HH

#include "sdsl/rank_support_v.hpp"
#include "sdsl/select_support_mcl.hpp"
#include "sdsl/bp_support_g.hpp"
#include "sdsl/util.hpp"
#include "globals.hh"
#include "BD_BWT_index/include/BD_BWT_index.hh"
#include "Interfaces.hh"

// Lowest marked ancestor support
// Does not own any of the data
class LMA_Support{
public:
    
    std::shared_ptr<Bitvector> marks;
    std::shared_ptr<Bitvector> bpr_marked_only;
    
    LMA_Support() {}
    
    // Marks assumes both open and close parentheses are marked
    LMA_Support(std::shared_ptr<Bitvector> marks,
                std::shared_ptr<Bitvector> bpr_marked_only) :
        marks(marks), bpr_marked_only(bpr_marked_only) {}
        
    // Takes the position of an open parenthesis in the bpr
    // Returns the position of the opening parenthesis of the
    // Nearest marked ancestor. The ancestor can be the node itself.
    // If no such ancestor exists, returns -1
    int64_t LMA(int64_t p){
        //assert((*bpr)[p] == 1); // Only works for open parenthesis
        if((*marks)[p]) return p;
        
        int64_t k = marks->rank(p);
        // The current node is a "virtual" node between positions k-1 and k in the marked only bpr
        
        if(k == 0 || k == bpr_marked_only->size()) return -1;
        
        // Four cases: we are in between (), ((, )) or )(
        int64_t open_in_marked_only;
        if((*bpr_marked_only)[k-1] == 1) // (( or ()
            open_in_marked_only = k-1;
        else if((*bpr_marked_only)[k] == 0) // ))
            open_in_marked_only = bpr_marked_only->find_open(k);        
        else{ // )(
            open_in_marked_only = bpr_marked_only->double_enclose(bpr_marked_only->find_open(k-1),k);
            if(open_in_marked_only == bpr_marked_only->size()){
                // This indicates that the double enclose did not find a common
                // ancestor (see sdsl documentation)
                return -1;
            }
        }
        
        return marks->select(open_in_marked_only+1); // Map back to the original bpr
    }
    
};


#endif
