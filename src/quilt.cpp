//
// This file is part of Velour and is distributed under the University of
// Illinois Open Source License.  See LICENSE.txt for details.
//

//
// quilt.cpp
//
//

#include "types.h"

static size_t hack_allocated_size = 0;
static size_t hack_actual_size = 0;

// Serial version                                                                             
void 
quilt_files(SeqGraph *graph, file_object_vector &work_to_do)
{
    //uintptr_t peak_seq_nodes = graph->node_count;
    //size_t peak_seq_node_mem = (g__SEQNODE_ALLOCATOR->seq_node_memory_ - g__SEQNODE_ALLOCATOR->seq_node_dead_memory_);
    fflush(stdout);
    for(file_object_vector::iterator itr = work_to_do.begin() ; itr != work_to_do.end() ; ++itr) {
        if( itr != work_to_do.begin() ) {
//            graph->resetFlags();
        }
        switch (itr->filetype) {
        case BUCKET:
            sg_load_bucket(graph, itr->filename);
            break;
        case QUILT:
            sg_load_quilt(graph, itr->filename);
            break;
        default:
            fprintf(stderr, "ERROR: Cannot use %s file %s as quilt input.  Exiting...\n",
                FILE_TYPES[itr->filetype], itr->filename);
            exit(EXIT_FAILURE);
        }

        if (!g__COMBINING && itr+1 == work_to_do.end()) {
            g__PSEUDO_NODES_PRESENT = false;
        }
        if( work_to_do.size() == 1 || itr != work_to_do.begin() ) {
            if( g__FULL_STATISTICS ) {
                sg_stat_components(graph, stdout);
            }
            //peak_seq_nodes = max(peak_seq_nodes, graph->node_count);
            //peak_seq_node_mem = max(peak_seq_node_mem, g__SEQNODE_ALLOCATOR->GetSeqLiveMemory());
            //printf("%zu MB current sequence node live memory (flow).\n", g__SEQNODE_ALLOCATOR->GetSeqLiveMemory() / (1024 * 1024));
#ifdef VERIFY
            graph->verify();
#endif
            //sg_remove_tips(graph);
#ifdef VERIFY
            graph->verify();
#endif
            //sg_concatenate(graph);
        }
#ifdef VERIFY
        graph->verify();
#endif

        //printf("Reading %s file %s\t%i subsequences found.\n", 
        //    FILE_TYPES[itr->filetype], itr->filename, num_subsequences);
        fflush(stdout);
    }

    graph->resetFlags();
    sg_remove_tips(graph);
    sg_concatenate(graph);

    printf("%"PRIuPTR" PREGRAPH sequence nodes.\n", graph->node_count);
    printf("%zu MB PREGRAPH sequence node live memory.\n", g__SEQNODE_ALLOCATOR->GetSeqLiveMemory() / (1024*1024));

    sg_pop_bubbles(graph);
    sg_covcutoff(graph);
    sg_pop_bubbles(graph);
    sg_remove_tips(graph);
    sg_pop_bubbles(graph);
    sg_remove_tips(graph);
    sg_concatenate(graph);

    printf("%"PRIuPTR" FINAL sequence nodes.\n", graph->node_count);
    printf("%zu MB FINAL sequence node live memory.\n", g__SEQNODE_ALLOCATOR->GetSeqLiveMemory() / (1024*1024));

    struct lambda1 {
        static void compute_sequence_sizes(SeqNode *node) {
            if (!isNodeDead<SeqNode>(node)) {
                hack_allocated_size += node->sequence.GetAllocatedBytes();
                hack_actual_size += node->sequence.GetLengthInBytes();
            }
        }
    };
    sg_for_each(graph, lambda1::compute_sequence_sizes);

    printf("%zu MB FINAL ALLOC sequence memory.\n", hack_allocated_size / (1024*1024));
    printf("%zu MB FINAL ACTUAL sequence memory.\n", hack_actual_size / (1024*1024));

    //peak_seq_nodes = max(peak_seq_nodes, graph->node_count);
    //peak_seq_node_mem = max(peak_seq_node_mem, g__SEQNODE_ALLOCATOR->GetSeqLiveMemory());
    //printf("%zu MB current sequence node live memory (flow).\n", g__SEQNODE_ALLOCATOR->GetSeqLiveMemory() / (1024 * 1024));

    /*
    for(uintptr_t i = 10000 ; i < graph->buckets ; ++ i) {
        SeqNode *trav;
        if ((trav = graph->head_table[i]) != NULL) {
            if (trav->right_count[0] > 4) {
                emit_scoop_graphviz(graph, trav, 200, "scoop.dot");
                break;
            }
        }
    }
    exit(1);

    graph->verify();
    sg_pop_bubbles(graph);
    graph->verify();
    sg_remove_tips(graph);
    sg_concatenate(graph);
    graph->verify();
    */

	graph->resetFlags();
    //printf("%"PRIuPTR" peak sequence nodes (flow).\n", peak_seq_nodes);
    //printf("%zu MB peak sequence node live memory (flow).\n", peak_seq_node_mem / (1024 * 1024));
}

