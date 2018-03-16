#ifndef BUILD_MODEL_HH
#define BUILD_MODEL_HH

#include "BPR_tools.hh"
#include "Parent_Support.hh"
#include "LMA_Support.hh"
#include "String_Depth_Support.hh"
#include "globals.hh"
#include "BWT_iteration.hh"
#include "Precalc.hh"
#include "context_marking.hh"
#include "Interfaces.hh"
#include "BPR_Colex_mapping.hh"
#include "UniBWT.h"
#include "RLEBWT.hh"
#include "Basic_bitvector.hh"
#include "RLE_bitvector.hh"
#include "logging.hh"
#include "build_model.hh"
#include <stack>
#include <vector>
#include <memory>

// All components of the model will be stored into G
// slt_it: iterator that gives all nodes that we want in the SLT
// rev_st_it: iterator that gives all nodes that we want in the rev ST.

void build_model(Global_Data& G, string& T, Context_Callback& context_formula,
                        Iterator& slt_it, Iterator& rev_st_it, bool run_length_coding, bool compute_string_depths){
        
    
    write_log("Building the BiBWT");
    G.bibwt = make_shared<BD_BWT_index<>>((uint8_t*)T.c_str());
    
    slt_it.set_index(G.bibwt.get());
    rev_st_it.set_index(G.bibwt.get());
    
    write_log("Building reverse suffix tree BPR and pruning marks");
    Build_REV_ST_BPR_And_Pruning_Callback revstbprcb;
    revstbprcb.init(*G.bibwt);
    iterate_with_callbacks(rev_st_it, &revstbprcb);
    Rev_st_topology RSTT = revstbprcb.get_result();
    G.rev_st_bpr = std::shared_ptr<Bitvector>(new Basic_bitvector(RSTT.bpr));
    
    G.rev_st_bpr->init_rank_10_support();
    G.rev_st_bpr->init_select_10_support();
    G.rev_st_bpr->init_bps_support();
    
    if(run_length_coding){
        write_log("Run length coding the pruning marks vector");
        G.pruning_marks = std::shared_ptr<Bitvector>(new RLE_bitvector(RSTT.pruning_marks));
    } else{
        G.pruning_marks = std::shared_ptr<Bitvector>(new Basic_bitvector(RSTT.pruning_marks));
    }
    
    G.pruning_marks->init_rank_support();
    G.pruning_marks->init_select_support();
    
    Pruned_Topology_Mapper mapper(G.rev_st_bpr, G.pruning_marks);
    
    write_log("Building SLT BPR");
    Build_SLT_BPR_Callback sltbprcb;
    sltbprcb.init(*G.bibwt);
    iterate_with_callbacks(slt_it, &sltbprcb);
    sdsl::bit_vector sdsl_slt_bpr = sltbprcb.get_result();
    if(run_length_coding){
        write_log("Run length coding SLT BPR");
        G.slt_bpr = std::shared_ptr<Bitvector>(new RLE_bitvector(sdsl_slt_bpr));
    } else{
        G.slt_bpr = std::shared_ptr<Bitvector>(new Basic_bitvector(sdsl_slt_bpr));
    }
    G.slt_bpr->init_rank_support();
    
    
    write_log("Marking maximal repeats and contexts");
    Rev_ST_Maximal_Marks_Callback revstmmcb;
    SLT_Maximal_Marks_Callback sltmmcb;
    revstmmcb.init(*G.bibwt, G.rev_st_bpr->size(), mapper);
    sltmmcb.init(*G.bibwt, sdsl_slt_bpr);
    context_formula.init(G.bibwt.get(), G.rev_st_bpr->size(), mapper);
    vector<Iterator_Callback*> marking_callbacks = {&revstmmcb, &sltmmcb, &context_formula};
    iterate_with_callbacks(slt_it, marking_callbacks);
    
    G.rev_st_maximal_marks = std::shared_ptr<Bitvector>(new Basic_bitvector(revstmmcb.get_result()));
    G.rev_st_maximal_marks->init_rank_support();
    
    G.rev_st_context_marks = std::shared_ptr<Bitvector>(new Basic_bitvector(context_formula.get_result()));
    G.rev_st_context_marks->init_rank_support();
    G.rev_st_context_marks->init_select_support();

    sdsl::bit_vector sdsl_slt_maxreps = sltmmcb.get_result();
    if(run_length_coding){
        write_log("Run length coding maximal repeats on SLT BPR");
        G.slt_maximal_marks = std::shared_ptr<Bitvector>(new RLE_bitvector(sdsl_slt_maxreps));
    } else{
        G.slt_maximal_marks = std::shared_ptr<Bitvector>(new Basic_bitvector(sdsl_slt_maxreps));
    }
    G.slt_maximal_marks->init_rank_support();
    G.slt_maximal_marks->init_select_support();
    
    write_log("Building the BPR of contexts only");
    G.rev_st_bpr_context_only = std::shared_ptr<Bitvector>(new Basic_bitvector(get_rev_st_bpr_context_only(&G)));
    G.rev_st_bpr_context_only->init_bps_support();
    
    if(compute_string_depths){
        cerr << "Not implemented error: compute string depths" << endl;
        throw(std::runtime_error("Not implemented error: compute string depths"));
        //G.string_depths = CA.get_rev_st_string_depths(G.bibwt, G.rev_st_bpr.size(), mapper);
    }
    
    // Store reverse BWT for scoring. Todo: reuse already computed bibwt
    
    string T_reverse(T.rbegin(), T.rend());
    if(run_length_coding){
        write_log("Computing run-length coded reverse BWT (todo: reuse the bibwt)");
        BWT* bwt = new RLEBWT<>((uint8_t*)T_reverse.c_str());
        G.revbwt = std::move(std::unique_ptr<BWT>(bwt)); // transfer ownership to G
    } else{
        write_log("Computing reverse BWT (todo: reuse the bibwt)");
        BWT* bwt = new Basic_BWT<>((uint8_t*)T_reverse.c_str());
        G.revbwt = std::move(std::unique_ptr<BWT>(bwt)); // transfer ownership to G
    }
    
    //cout << G.toString() << endl;
    
}

#endif