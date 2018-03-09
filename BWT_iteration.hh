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
        iteration_stack.push(Stack_frame(Interval_pair(0,index->size()-1,0,index->size()-1),0,true)); // Suppose empty string is always a maxrep (special case)
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
         
         // Iterate alphabet in reverse lexicographic order, so the smallest is pushed to the
         // stack the last, so the iteration is done in lexicographic order
         for(int64_t i = interval_data.n_distinct_symbols-1; i >= 0; i--){
            Interval_pair I2 = index->right_extend(top.intervals,interval_data,i);
            if(I2.forward.size() != 0 && index->is_left_maximal(I2)){
                iteration_stack.push(Stack_frame(I2, top.depth+1, index->is_right_maximal(I2)));
            }
        }
        
        index->compute_bwt_interval_data(top.intervals.forward, interval_data);
        //std::sort(interval_data.symbols.begin(), interval_data.symbols.begin() + interval_data.n_distinct_symbols);
         
        // Push all rev st children that are leaves
         for(int64_t i = interval_data.n_distinct_symbols-1; i >= 0; i--){
            Interval_pair I2 = index->left_extend(top.intervals,interval_data,i);
            if(I2.reverse.size() == 1){
                iteration_stack.push(Stack_frame(I2, -1, false)); // Don't know the depth so just put -1
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
        
        // Push the empty string (suppose it is a maxrep (holds as long as the text is not an empty string?)
        iteration_stack.push(Stack_frame(Interval_pair(0,index->size()-1,0,index->size()-1),0, true));
    }
    
    virtual bool next(){
         // 1) Take an interval off the top of a stack
         // 2) Push all right-maximal left extensions to stack
         if(iteration_stack.empty()) return false;
         
         top = iteration_stack.top();
         iteration_stack.pop();
         
         index->compute_bwt_interval_data(top.intervals.forward, interval_data);
         //std::sort(interval_data.symbols.begin(), interval_data.symbols.begin() + interval_data.n_distinct_symbols);
         
         if(top.depth <= depth_bound - 1){
            // Iterate alphabet in reverse lexicographic order, so the smallest is pushed to the
            // stack the last, so the iteration is done in lexicographic DFS order
            for(int64_t i = interval_data.n_distinct_symbols-1; i >= 0; i--){
                Interval_pair I2 = index->left_extend(top.intervals,interval_data,i);
                if(I2.forward.size() != 0 && index->is_right_maximal(I2)){
                    iteration_stack.push(Stack_frame(I2, top.depth+1, index->is_left_maximal(I2)));
                }
            }
         }
                
        return true;
    }    
};

class SLT_Iterator : public Depth_Bounded_SLT_Iterator{
public:
        
    SLT_Iterator() : Depth_Bounded_SLT_Iterator(1e18) {}
    SLT_Iterator(BIBWT* index) : Depth_Bounded_SLT_Iterator(index, 1e18) {}
           
};

class Rev_ST_Depth_Bounded_Maxrep_Iterator : public Iterator{

    private:
    std::vector<Stack_frame> iteration_stack;
    
    public:
    
    
    Stack_frame top;
    BIBWT* index;
    typename BIBWT::Interval_Data interval_data;
    int64_t depth_bound;
    
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
        
        // Clear the stack
        while(!iteration_stack.empty()) iteration_stack.pop_back();
        
        // Push the empty string (suppose it is a maxrep)
        iteration_stack.push_back(Stack_frame(Interval_pair(0,index->size()-1,0,index->size()-1),0, true));
    }
    
    bool next(){
         start:
         if(iteration_stack.empty()) return false;
         
         top = iteration_stack.back();
         iteration_stack.pop_back();
         // top has depth at most equal to the depth bound
         
         index->compute_bwt_interval_data(top.intervals.forward, interval_data);
         
         bool leftmax = interval_data.n_distinct_symbols >= 2;
         bool rightmax = index->is_right_maximal(top.intervals);
         top.is_maxrep = leftmax && rightmax;
        
        if(!rightmax){
            // Just have gone outside of the maxrep tree
            return true;
        }
         
        if(rightmax && !leftmax && top.depth == depth_bound){
            // The left-saturation is maxrep, but will not reach it because of
            // the depth bound -> return a node inside the edge
            return true;
        }
        
        if(top.depth <= depth_bound - 1){
            // Iterate alphabet in reverse lexicographic order, so the smallest is pushed to the
            // stack the last, so the iteration is done in lexicographic DFS order
            for(int64_t i = interval_data.n_distinct_symbols-1; i >= 0; i--){
                //char c = interval_data.symbols[i];
                Interval_pair I2 = index->left_extend(top.intervals, interval_data, i);
                if(I2.forward.size() != 0){
                    iteration_stack.push_back(Stack_frame(I2, top.depth+1, false)); // is_maxrep will be computed when the frame is popped
                }
            }
        }
        
        if(leftmax) return true; // Is a maxrep within the depth bound
        else goto start; // Could use recursion but goto is faster and does not increase the size of the stack
    }
};


class Rev_ST_Maxrep_Iterator : public Rev_ST_Depth_Bounded_Maxrep_Iterator{
    
    public:
    
    Stack_frame top;
    BIBWT* index;
    typename BIBWT::Interval_Data interval_data;
    
    Rev_ST_Maxrep_Iterator() : Rev_ST_Depth_Bounded_Maxrep_Iterator(1e18) {}
    Rev_ST_Maxrep_Iterator(BIBWT* index) : Rev_ST_Depth_Bounded_Maxrep_Iterator(index, 1e18) {}
    

};


#endif

