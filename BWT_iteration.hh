#ifndef BWT_ITERATION_HH
#define BWT_ITERATION_HH

#include "Interfaces.hh"

// Iterators that give all nodes that should be included in the topology.

class Rev_ST_Iterator : public Iterator{

    public:
    
    std::stack<Stack_frame> iteration_stack;
    Stack_frame top;
    BIBWT* index;
    typename BIBWT::Interval_Data interval_data;
    
    Rev_ST_Iterator() {}
    Rev_ST_Iterator(BIBWT* index) : index(index) {}
    
    virtual void set_index(BIBWT* index){
        this->index = index;
    }
    
    virtual Iterator::Stack_frame get_top(){
        return top;
    }
    
    virtual void init(){
        // Make space
        interval_data.symbols.resize(index->get_alphabet().size());
        interval_data.ranks_start.resize(index->get_alphabet().size());
        interval_data.ranks_end.resize(index->get_alphabet().size());
        
        // Clear the stack
        while(!iteration_stack.empty()) iteration_stack.pop();
        
        // Push the empty string 
        iteration_stack.push(Stack_frame(Interval_pair(0,index->size()-1,0,index->size()-1),0));
    }
    
    virtual bool next(){
         // 1) Take an interval off the top of a stack
         // 2) If it is a leaf, return
         // 3) Push all right-extensions that are left-maximal to the stack
         //    and all left-extensions that are leaves.
         // WARNING: Depths are not accurate for leaves
         
         if(iteration_stack.empty()) return false;
         
         top = iteration_stack.top();
         iteration_stack.pop();
         
         if(top.intervals.reverse.size() == 1) return true; // Leaf
         
         index->compute_rev_bwt_interval_data(top.intervals.reverse, interval_data);
         std::sort(interval_data.symbols.begin(), interval_data.symbols.begin() + interval_data.n_distinct_symbols);
         
         // Iterate alphabet in reverse lexicographic order, so the smallest is pushed to the
         // stack the last, so the iteration is done in lexicographic order
         for(int64_t i = interval_data.n_distinct_symbols-1; i >= 0; i--){
            char c = interval_data.symbols[i];
            Interval_pair I2 = index->right_extend(top.intervals,c);
            if(I2.forward.size() != 0 && index->is_left_maximal(I2)){
                iteration_stack.push(Stack_frame(I2, top.depth+1));
            }
        }
        
        index->compute_bwt_interval_data(top.intervals.forward, interval_data);
        std::sort(interval_data.symbols.begin(), interval_data.symbols.begin() + interval_data.n_distinct_symbols);
         
        // Push all rev st children that are leaves
         for(int64_t i = interval_data.n_distinct_symbols-1; i >= 0; i--){
            char c = interval_data.symbols[i];
            Interval_pair I2 = index->left_extend(top.intervals,c);
            if(I2.reverse.size() == 1){
                iteration_stack.push(Stack_frame(I2, -1)); // Don't know the depth so just put -1
            }
        }
        
        return true;
    }
};

class SLT_Iterator : public Iterator{
public:
    
    std::stack<Stack_frame> iteration_stack;
    Stack_frame top;
    typename BIBWT::Interval_Data interval_data;
    BIBWT* index;
    
    SLT_Iterator() {} // Need to set index later
    SLT_Iterator(BIBWT* index) : index(index) {}
           
    virtual void set_index(BIBWT* index){
        this->index = index;
    }
    
    virtual Iterator::Stack_frame get_top(){
        return top;
    }
        
    virtual void init(){
        
        // Make space
        interval_data.symbols.resize(index->get_alphabet().size());
        interval_data.ranks_start.resize(index->get_alphabet().size());
        interval_data.ranks_end.resize(index->get_alphabet().size());
        
        // Clear the stack
        while(!iteration_stack.empty()) iteration_stack.pop();
        
        // Push the empty string 
        iteration_stack.push(Stack_frame(Interval_pair(0,index->size()-1,0,index->size()-1),0));
    }
    
    virtual bool next(){
         // 1) Take an interval off the top of a stack
         // 2) Push all right-maximal left extensions to stack
         if(iteration_stack.empty()) return false;
         
         top = iteration_stack.top();
         iteration_stack.pop();
         
         index->compute_bwt_interval_data(top.intervals.forward, interval_data);
         std::sort(interval_data.symbols.begin(), interval_data.symbols.begin() + interval_data.n_distinct_symbols);
         
         // Iterate alphabet in reverse lexicographic order, so the smallest is pushed to the
         // stack the last, so the iteration is done in lexicographic DFS order
         for(int64_t i = interval_data.n_distinct_symbols-1; i >= 0; i--){
            char c = interval_data.symbols[i];
            Interval_pair I2 = index->left_extend(top.intervals,c);
            if(I2.forward.size() != 0 && index->is_right_maximal(I2)){
                iteration_stack.push(Stack_frame(I2, top.depth+1));
            }
        }
        
        return true;
    }    
};

class Depth_Bounded_SLT_Iterator : public Iterator{
public:
    
    std::stack<Stack_frame> iteration_stack;
    Stack_frame top;
    typename BIBWT::Interval_Data interval_data;
    BIBWT* index;
    int64_t depth_bound;
    string label; // debug
    
    Depth_Bounded_SLT_Iterator(int64_t depth_bound) : depth_bound(depth_bound) {}
    Depth_Bounded_SLT_Iterator(BIBWT* index, int64_t depth_bound) : index(index), depth_bound(depth_bound) {}
        
    virtual void set_index(BIBWT* index){
        this->index = index;
    }

    virtual Iterator::Stack_frame get_top(){
        return top;
    }
        
    virtual void init(){
        
        label = ""; // debug
        
        // Make space
        interval_data.symbols.resize(index->get_alphabet().size());
        interval_data.ranks_start.resize(index->get_alphabet().size());
        interval_data.ranks_end.resize(index->get_alphabet().size());
        
        // Clear the stack
        while(!iteration_stack.empty()) iteration_stack.pop();
        
        // Push the empty string 
        iteration_stack.push(Stack_frame(Interval_pair(0,index->size()-1,0,index->size()-1),0));
    }
    
    virtual bool next(){
         // 1) Take an interval off the top of a stack
         // 2) Push all right-maximal left extensions to stack
         if(iteration_stack.empty()) return false;
         
         top = iteration_stack.top();
         iteration_stack.pop();
         
         index->compute_bwt_interval_data(top.intervals.forward, interval_data);
         std::sort(interval_data.symbols.begin(), interval_data.symbols.begin() + interval_data.n_distinct_symbols);
         
         if(top.depth <= depth_bound - 1){
            // Iterate alphabet in reverse lexicographic order, so the smallest is pushed to the
            // stack the last, so the iteration is done in lexicographic DFS order
            for(int64_t i = interval_data.n_distinct_symbols-1; i >= 0; i--){
                char c = interval_data.symbols[i];
                Interval_pair I2 = index->left_extend(top.intervals,c);
                if(I2.forward.size() != 0 && index->is_right_maximal(I2)){
                    iteration_stack.push(Stack_frame(I2, top.depth+1));
                }
            }
         }
                
        return true;
    }    
};

class Rev_ST_Maxrep_Iterator : public Iterator{
    
    private:
        
    std::vector<Stack_frame> iteration_stack;
    std::vector<Stack_frame> maxrep_left_extensions;
    
    public:
    
    Stack_frame top;
    BIBWT* index;
    typename BIBWT::Interval_Data interval_data;
    
    Rev_ST_Maxrep_Iterator() {}
    Rev_ST_Maxrep_Iterator(BIBWT* index) : index(index) {}
    
    virtual void set_index(BIBWT* index){
        this->index = index;
    }
    
    virtual Iterator::Stack_frame get_top(){
        return top;
    }
    
    virtual void init(){
        // Make space
        interval_data.symbols.resize(index->get_alphabet().size());
        interval_data.ranks_start.resize(index->get_alphabet().size());
        interval_data.ranks_end.resize(index->get_alphabet().size());
        
        // Clear the stack
        while(!iteration_stack.empty()) iteration_stack.pop_back();
        while(!maxrep_left_extensions.empty()) maxrep_left_extensions.pop_back();
        
        // Push the empty string 
        iteration_stack.push_back(Stack_frame(Interval_pair(0,index->size()-1,0,index->size()-1),0));
    }
    
    // Iterates only maxreps and all left extensions of maxreps
    virtual bool next(){
        if(iteration_stack.empty()){
            // If all rev st nodes have been considered, unload the maxrep left extension stack
            if(!maxrep_left_extensions.empty()){
                top = maxrep_left_extensions.back();
                maxrep_left_extensions.pop_back();
                return true; // Return control
            }
            return false; // Both stacks empty. Done.
        }
        
        top = iteration_stack.back();
        iteration_stack.pop_back();
        
        bool leftmax = index->is_left_maximal(top.intervals);
        bool rightmax = index->is_right_maximal(top.intervals);
        bool maxrep = leftmax && rightmax;
         
        if(leftmax){
            // Push right-extensions
            index->compute_rev_bwt_interval_data(top.intervals.reverse, interval_data);
            std::sort(interval_data.symbols.begin(), interval_data.symbols.begin() + interval_data.n_distinct_symbols);
            for(int64_t i = interval_data.n_distinct_symbols-1; i >= 0; i--){
                char c = interval_data.symbols[i];
                Interval_pair I2 = index->right_extend(top.intervals,c);
                if(I2.reverse.size() != 0){
                   iteration_stack.push_back(Stack_frame(I2, top.depth+1));
                   // Note: if we are inside an edge, the depth refers to the locus we are
                   // currently, not to the depth of the node with the given colex-interval
                }
            }
        }
        
        if(maxrep){
            // Push left-extensions
            index->compute_bwt_interval_data(top.intervals.forward, interval_data);
            std::sort(interval_data.symbols.begin(), interval_data.symbols.begin() + interval_data.n_distinct_symbols);
            for(int64_t i = interval_data.n_distinct_symbols-1; i >= 0; i--){
                char c = interval_data.symbols[i];
                Interval_pair I2 = index->left_extend(top.intervals,c);
                if(I2.forward.size() != 0){
                    if(!index->is_right_maximal(I2)){ // Left extension of a maxrep
                        maxrep_left_extensions.push_back(Stack_frame(I2, top.depth+1));
                    }
                }
            }
        }
        
        
         
        // Return control if top is a maxrep
        if(maxrep){
            return true;
        }
        else return next(); // Keep iterating
    }
};

class Rev_ST_Depth_Bounded_Maxrep_Iterator : public Iterator{

    private:
    std::stack<Stack_frame> iteration_stack;
    std::stack<Stack_frame> maxrep_left_extensions;
    
    public:
    
    
    Stack_frame top;
    BIBWT* index;
    typename BIBWT::Interval_Data interval_data;
    int64_t depth_bound;
    string label; // debug
    
    Rev_ST_Depth_Bounded_Maxrep_Iterator(int64_t depth_bound) : depth_bound(depth_bound) {}
    Rev_ST_Depth_Bounded_Maxrep_Iterator(BIBWT* index, int64_t depth_bound) : index(index), depth_bound(depth_bound) {}
    
    virtual void set_index(BIBWT* index){
        this->index = index;
    } 
    
    virtual Iterator::Stack_frame get_top(){
        return top;
    }

    virtual void init(){
        // Make space
        interval_data.symbols.resize(index->get_alphabet().size());
        interval_data.ranks_start.resize(index->get_alphabet().size());
        interval_data.ranks_end.resize(index->get_alphabet().size());
        label = "";
        
        // Clear the stack
        while(!iteration_stack.empty()) iteration_stack.pop();
        while(!maxrep_left_extensions.empty()) maxrep_left_extensions.pop();

        
        // Push the empty string 
        iteration_stack.push(Stack_frame(Interval_pair(0,index->size()-1,0,index->size()-1),0));
        maxrep_left_extensions.push(Stack_frame(Interval_pair(0,index->size()-1,0,index->size()-1),0));
        // The empty string is not a left-extension of a maxrep but we need it as a special 
        // case so that our tree has a root
    }
    
    // Iterates only maxreps and all left extensions of maxreps
    virtual bool next(){
        if(!maxrep_left_extensions.empty()){
            top = maxrep_left_extensions.top();
            maxrep_left_extensions.pop();
            return true; // Return control
        } else{
            if(iteration_stack.empty()) return false; // done
            
            Stack_frame now = iteration_stack.top();
            iteration_stack.pop();
            
            bool leftmax = index->is_left_maximal(now.intervals);
            bool rightmax = index->is_right_maximal(now.intervals);
            bool maxrep = leftmax && rightmax;
            
            if(leftmax){
                // Push right-extensions to iteration stack
                index->compute_rev_bwt_interval_data(now.intervals.reverse, interval_data);
                std::sort(interval_data.symbols.begin(), interval_data.symbols.begin() + interval_data.n_distinct_symbols);
                for(int64_t i = interval_data.n_distinct_symbols-1; i >= 0; i--){
                    char c = interval_data.symbols[i];
                    Interval_pair I2 = index->right_extend(now.intervals,c);
                    if(I2.reverse.size() != 0){
                    iteration_stack.push(Stack_frame(I2, now.depth+1));
                    // Note: if we are inside an edge, the depth refers to the locus we are
                    // currently, not to the depth of the node with the given colex-interval
                    }
                }
            }
            
            if(maxrep){
                // Push left-extensions to left-extension stack
                index->compute_bwt_interval_data(now.intervals.forward, interval_data);
                std::sort(interval_data.symbols.begin(), interval_data.symbols.begin() + interval_data.n_distinct_symbols);
                for(int64_t i = interval_data.n_distinct_symbols-1; i >= 0; i--){
                    char c = interval_data.symbols[i];
                    Interval_pair I2 = index->left_extend(now.intervals,c);
                    if(I2.forward.size() != 0){
                        if(now.depth < depth_bound){ // Left extension of a maxrep
                            maxrep_left_extensions.push(Stack_frame(I2, now.depth+1));
                        }
                    }
                }
            }
            
            return next();
        }
    }
};

#endif

