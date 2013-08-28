//
// This file is part of Velour and is distributed under the University of
// Illinois Open Source License.  See LICENSE.txt for details.
//

//
// seqGraph.cpp
//

#include "types.h"
#include "seqGraph.h"

static uintptr_t p__SG_NODES_BUILT = 0;
static uintptr_t p__SG_NODES_CONCAT = 0;
static uintptr_t p__SG_NODES_LOADED = 0;
static uintptr_t p__SG_NODES_EMITTED = 0;

static flow_nodelist_t *p__GROWNODELIST = NULL;

//
// helpers
//


static inline SeqNode *
load_seq_node(FILE* file, int fileFormat)
{
    SeqNode_StackAllocated memory;
    SeqNode *stack_node = new (&memory) SeqNode(Sequence::MAX_BASES);
    stack_node->loadFromFile(file, fileFormat);

    if (stack_node->sequence.get_length() == 0) { // for detecting EOF
        stack_node->~SeqNode(); // explicit destruction
        return NULL;
    } else {
        uint16_t alloc_bases = stack_node->sequence.GetAllocatedBases();
#ifdef VELOUR_TBB
        SeqNode *mem = tls_seqnode_allocator->AllocateNodeMemory( stack_node->GetNodeAllocatedBytes() );
#else
        SeqNode *mem = g__SEQNODE_ALLOCATOR->AllocateNodeMemory( stack_node->GetNodeAllocatedBytes() );
#endif // VELOUR_TBB
        SeqNode *sg_node = new (mem) SeqNode(stack_node, alloc_bases);
        stack_node->~SeqNode(); // explicit destruction
        ++ p__SG_NODES_LOADED;
        return sg_node;
    }
}

static inline SeqNode *
load_seq_node_mmap(char *mmap_start, size_t &amount)
{
    SeqNode_StackAllocated memory;
    SeqNode *stack_node = new (&memory) SeqNode(Sequence::MAX_BASES);
    amount = stack_node->loadFromMemory(mmap_start);

    if (stack_node->sequence.get_length() == 0) { // for detecting EOF
        stack_node->~SeqNode(); // explicit destruction
        return NULL;
    } else {
        uint16_t alloc_bases = stack_node->sequence.GetAllocatedBases();
#ifdef VELOUR_TBB
        SeqNode *mem = tls_seqnode_allocator->AllocateNodeMemory( stack_node->GetNodeAllocatedBytes() );
#else
        SeqNode *mem = g__SEQNODE_ALLOCATOR->AllocateNodeMemory( stack_node->GetNodeAllocatedBytes() );
#endif // VELOUR_TBB
        SeqNode *sg_node = new (mem) SeqNode(stack_node, alloc_bases);
        stack_node->~SeqNode(); // explicit destruction
        ++ p__SG_NODES_LOADED;
        return sg_node;
    }
}

static inline void emit_seq_node(SeqNode *sg_node, FILE* file, int fileFormat)
{
    ++ p__SG_NODES_EMITTED;
    sg_node->emitToFile(file, fileFormat, p__SG_NODES_EMITTED);
}

static inline void
copy_counters(KmerNode *kg_node, SeqNode *sg_node, int kg_direction, int sg_direction)
{
	four_counter_t *kg_side = kg_direction == GO_RIGHT ? &kg_node->right_side : &kg_node->left_side;
	four_color_t *kg_cside = kg_direction == GO_RIGHT ? &kg_node->right_side_colors : &kg_node->left_side_colors;
	four_counter_t *sg_side = sg_direction == GO_RIGHT ? &sg_node->right_side : &sg_node->left_side;
	four_color_t *sg_cside = sg_direction == GO_RIGHT ? &sg_node->right_side_colors : &sg_node->left_side_colors;
	assert( *sg_side == 0 );
	assert( *sg_cside == 0 );

	if (kg_direction == sg_direction) {
		*sg_side = *kg_side;
		*sg_cside = *kg_cside;
	} else { // negate counters and complement, as sense is reversed
		counter_t *kg_counters = static_cast<counter_t *>(static_cast<void *>(kg_side));
		counter_t *sg_counters = static_cast<counter_t *>(static_cast<void *>(sg_side));
		Color *kg_color = static_cast<Color *>(static_cast<void *>(kg_cside));
		Color *sg_color = static_cast<Color *>(static_cast<void *>(sg_cside));
		for (int i=0; i < 4; ++i) {
			sg_counters[i] = -kg_counters[COMPLEMENT(i)];
            sg_color[i] = kg_color[COMPLEMENT(i)];
		}
	}

    // TODO: bitconnections?
}

static void
mark_all_neighbors_as_merge_starts(KmerNode *kg_node, int direction)
{
    bool moving_right = (direction == GO_RIGHT);
    counter_t *counters = moving_right ? kg_node->right_count : kg_node->left_count;
 
#ifndef SMALL_NODES
    KmerNode **pointers = moving_right ? kg_node->right : kg_node->left;
#endif

    for (int i = 0 ; i < 4 ; ++i ) {
        counter_t count = counters[i];
        if (count == 0) continue;
        bool sense_changed = (count < 0);
        color_t *colors = moving_right ? kg_node->right_color : kg_node->left_color;
        KmerNode *neighbor = NULL;
        if (colors[i] == 0) {
#ifdef SMALL_NODES
            neighbor = g__KG_HASHTABLE->findNextNode(kg_node, i, moving_right);
#else
            neighbor = pointers[i];
#endif
        } else {
            continue;
        }
        assert( neighbor != NULL );
        //if (isNodeMerged<kg_node_t>(neighbor)) continue; // commented out for performance, FIXME verify benefit
        if (moving_right ^ sense_changed) {
            setNodeMergingRight<KmerNode>(neighbor);
        } else {
            setNodeMergingLeft<KmerNode>(neighbor);
        }
    }
}

// NOTE: assumes sg_node has been allocated MAX_BASES sequence capacity
static void
extend_contig(KmerNode *__arg_kg_node, SeqNode *__arg_sg_node, int __arg_kg_direction, int __arg_sg_direction)
{
    KmerNode *kg_node = __arg_kg_node;
    SeqNode *sg_node = __arg_sg_node;
    int kg_direction = __arg_kg_direction;
    int sg_direction = __arg_sg_direction;

    unsigned iteration_counter = 1;

    while(1) {  // loop instead of recursion

    assert( !isNodeMerged<KmerNode>(kg_node) );
    setNodeMerged<KmerNode>(kg_node);

    if (!sg_node->sequence.CanAddBases(1)) {
        copy_counters(kg_node, sg_node, kg_direction, sg_direction);
        mark_all_neighbors_as_merge_starts(kg_node, kg_direction);
        return;
    }

    bool kg_moving_right = (kg_direction == GO_RIGHT);
    counter_t *kg_counters = kg_moving_right ? kg_node->right_count : kg_node->left_count;

#ifdef SMALL_NODES
    int valid_dir = valid_single_successor(kg_counters);
    if (valid_dir == MULTIPLE_SUCCESSORS) {
        copy_counters(kg_node, sg_node, kg_direction, sg_direction);
        mark_all_neighbors_as_merge_starts(kg_node, kg_direction);
        return;
    }
    if (valid_dir == NO_SUCCESSORS) {
        // no counters to copy to sg_node
        return;
    }

    // single direction
    assert((valid_dir >= 0) && (valid_dir < 4));
    bool sense_changed = kg_counters[valid_dir] < 0;
    color_t *colors = kg_moving_right ? kg_node->right_color : kg_node->left_color;
    KmerNode *next = NULL;
    if (colors[valid_dir] == 0) {
        next = g__KG_HASHTABLE->findNextNode(kg_node, valid_dir, kg_moving_right);
    } else {
        copy_counters(kg_node, sg_node, kg_direction, sg_direction);
        return;
    }
    assert( next != NULL );

#ifdef VELOUR_TBB
# ifdef VERIFY_THREADSAFE
    assert( next->claim_tid == 0 || next->claim_tid == tls_thread_index ); // ensure next node was claimed by this thread
# endif
#endif

#else
    // TODO: need to setup the pointers for sg_node somewhere... 
    unsigned bitconnections = kg_node->bitconnections;
    unsigned connex = kg_moving_right ? RIGHT(bitconnections) : LEFT(bitconnections);
    unsigned connex_count = CONNECTION_COUNT_MAP[connex];
    if (connex_count != 1) { // can't extend
        if (connex_count != 0) {
            copy_counters(kg_node, sg_node, kg_direction, sg_direction);
            mark_all_neighbors_as_merge_starts(kg_node, kg_direction);
        }
        return;
    }

    // single direction
    kg_node_t **pointers = kg_moving_right ? kg_node->right : kg_node->left;
    int valid_dir = valid_single_successor(counters);
    assert((valid_dir >= 0) && (valid_dir < 4));
    kg_node_t *next = pointers[valid_dir];
    bool sense_changed = counters[valid_dir] < 0;
#endif

    if (isNodeMerged(next)) { // can't extend
        copy_counters(kg_node, sg_node, kg_direction, sg_direction);
        return;
    }

    bool kg_next_moving_right = kg_moving_right ^ sense_changed;
    counter_t *kg_next_back_counters = kg_next_moving_right ? next->left_count : next->right_count;
    int next_back_valid = valid_single_successor(kg_next_back_counters);
    if (next_back_valid == MULTIPLE_SUCCESSORS) {
        //|| next_back_valid == NO_SUCCESSORS   ) {  // this can happen with small nodes -- FIXME check this
        copy_counters(kg_node, sg_node, kg_direction, sg_direction);
        return;   // can't extend
    }
      

    // check for: - single copy clipping marking, which inhibits concatenation
    //            - start at base of Y, conscious of possible erosion
    //            - copy count disagreement, for proper single-copy erosion behavior
    //
    // XXX alternative implementation: traverse and check for long _tip_ that shouldn't clip -- only then can ignore/remove clip mark
    if (g__PSEUDO_NODES_PRESENT) {
        assert( isNodeWeight_AllNonSingleCopy(sg_node) || isNodeWeight_AllSingleCopy(sg_node) || isNodeWeight_None(sg_node) );

        // single copy clipping marking
        if ( abs(kg_counters[valid_dir]) == CLIP_SINGLECOPY_COUNTER_VALUE ||
             abs(kg_next_back_counters[next_back_valid]) == CLIP_SINGLECOPY_COUNTER_VALUE ) {
            copy_counters(kg_node, sg_node, kg_direction, sg_direction); // OPT: could just copy/merge-start the one arc instead of all?
            mark_all_neighbors_as_merge_starts(kg_node, kg_direction);
            return;
        }

        counter_t edge_weight = abs(kg_counters[valid_dir]);
        if (edge_weight == CLIP_SINGLECOPY_COUNTER_VALUE) edge_weight = 1;
        assert( edge_weight > 0 );
                
        counter_t *kg_next_forward_counters = kg_next_moving_right ? next->right_count : next->left_count;
        counter_t next_forward_edges_weight = get_counter_sum(kg_next_forward_counters); // XXX: or max?
                
        counter_t *sg_back_counters = (sg_direction == GO_RIGHT) ? sg_node->left_count : sg_node->right_count;
        bool is_sg_back_multiple_successors = (valid_single_successor(sg_back_counters) == MULTIPLE_SUCCESSORS);

        // base of Y: erosion unclear
        if (edge_weight == 1) {
            if (valid_single_successor(kg_next_forward_counters) == MULTIPLE_SUCCESSORS) {
                /*bool kg_next_forward_single_copy = true;
                for (int i=0; i < 4; ++i) {
                    counter_t value = abs(kg_next_forward_counters[i]);
                    if (value != 0 && value != 1 && value != CLIP_SINGLECOPY_COUNTER_VALUE) {
                        kg_next_forward_single_copy = false;
                    }
                }
                if (kg_next_forward_single_copy) { // all single copy edges*/
                    copy_counters(kg_node, sg_node, kg_direction, sg_direction); // OPT: could just copy/merge-start the one arc instead of all
                    mark_all_neighbors_as_merge_starts(kg_node, kg_direction);
                    return;
                //}
            }
            //if (sg_node->min_edge_weight == 0) { // optimization: a non-zero merged node will always be ok
                if (valid_single_successor(sg_back_counters) == MULTIPLE_SUCCESSORS) {
                    /*bool sg_back_single_copy = true;
                    for (int i=0; i < 4; ++i) {
                        counter_t value = abs(sg_back_counters[i]);
                        if (value != 0 && value != 1 && value != CLIP_SINGLECOPY_COUNTER_VALUE) {
                            sg_back_single_copy = false;
                        }
                    }
                    if (sg_back_single_copy) { // all single copy backward edges*/
                        copy_counters(kg_node, sg_node, kg_direction, sg_direction); // OPT: could just copy/merge-start the one arc instead of all
                        mark_all_neighbors_as_merge_starts(kg_node, kg_direction);
                        return;
                    //}
                }
            //}
        }
        counter_t sg_back_edge_weight = get_counter_sum(sg_back_counters); // XXX: or max?

        // known: sg_back is not single copy multiple back successor, thus can forward merge any copy value
        
        // copy count disagreement: single vs multi
        if ((isNodeWeight_None(sg_node) && !is_sg_back_multiple_successors && edge_weight == 1 && (sg_back_edge_weight > 1 || next_forward_edges_weight > 1)) ||
            (isNodeWeight_HasSingleCopy(sg_node) && (edge_weight > 1 || next_forward_edges_weight > 1)) ||
            (isNodeWeight_AllNonSingleCopy(sg_node) && edge_weight == 1) ) {
            copy_counters(kg_node, sg_node, kg_direction, sg_direction); // OPT: could just copy/merge-start the one arc instead of all
            mark_all_neighbors_as_merge_starts(kg_node, kg_direction);
            return;
        }
    }

    // otherwise, we can extend.  Need to add base, and then call recursively.
	if (sg_direction == GO_RIGHT) {
		sg_node->sequence.AppendBase_Unsafe(kg_moving_right ? valid_dir : COMPLEMENT(valid_dir));
#ifdef VERIFY
		Kmer kmer = sg_node->sequence.GetTailKmer(g__FULLKMER_LENGTH);
		//VERIFY_KMER(kmer);
#endif
	} else {
		sg_node->sequence.PrependBase_Unsafe(kg_moving_right ? COMPLEMENT(valid_dir) : valid_dir);
#ifdef VERIFY
		Kmer kmer = sg_node->sequence.GetHeadKmer(g__FULLKMER_LENGTH);
		//VERIFY_KMER(kmer);
#endif
	}
    counter_t edge_weight = abs(kg_counters[valid_dir]);
    if (edge_weight == CLIP_SINGLECOPY_COUNTER_VALUE) edge_weight = 1;
    assert( edge_weight > 0 );
    updateNodeWeight(sg_node, edge_weight);
    sg_node->kmer_occurrences += next->kmer_occurrences;
    if (g__PSEUDO_NODES_PRESENT) {
        assert( isNodeWeight_AllNonSingleCopy(sg_node) || isNodeWeight_AllSingleCopy(sg_node) || isNodeWeight_None(sg_node) );
    }

    // setup variables for next iteration of loop, i.e. recursive function call
    //   WAS: extend_contig(next, sg_node, kg_direction ^ sense_changed, sg_direction);
    kg_node = next;
    //sg_node = sg_node;
    kg_direction = kg_direction ^ sense_changed;
    //sg_direction = sg_direction;

    // track the "recursion depth"
    ++ iteration_counter;
    }
}

static void
build_node(KmerNode *kg_node, int direction, SeqGraph *sg_hashtable, FILE* file, int fileFormat)
{
    SeqNode_StackAllocated memory;
    SeqNode * stack_seqnode = new (&memory) SeqNode(Sequence::MAX_BASES);

    // mark the other direction for collection on the next pass
    mark_all_neighbors_as_merge_starts(kg_node, !direction);

    // copy counters
    copy_counters(kg_node, stack_seqnode, !direction, !direction);

	// set initial sequence
    stack_seqnode->sequence.InitializeWithKmer_Unsafe(kg_node->kmer, g__FULLKMER_LENGTH);

    // set initial kmer occurrences count
    stack_seqnode->kmer_occurrences = kg_node->kmer_occurrences;

	// attempt to extend sequence
    extend_contig(kg_node, stack_seqnode, direction, direction);

	// emit the sequence node
    if (fileFormat != NOFORMAT ) {
        emit_seq_node(stack_seqnode, file, fileFormat);
    }

	// insert node in sequence graph
    if (sg_hashtable != NULL) {
        uint16_t alloc_bases = stack_seqnode->sequence.GetAllocatedBases();
#ifdef VELOUR_TBB
        SeqNode *mem = tls_seqnode_allocator->AllocateNodeMemory( stack_seqnode->GetNodeAllocatedBytes() );
#else
        SeqNode *mem = g__SEQNODE_ALLOCATOR->AllocateNodeMemory( stack_seqnode->GetNodeAllocatedBytes() );
#endif // VELOUR_TBB
        SeqNode *sg_node = new (mem) SeqNode(stack_seqnode, alloc_bases);
        sg_hashtable->insertNode(sg_node);
    }

    stack_seqnode->~SeqNode(); // explicit destruction

    ++ p__SG_NODES_BUILT;
}

// XXX: either propagate kmer hashtable pointer into functions, or remove the variable
struct KmerConcatPass1 {
    KmerConcatPass1(KmerGraph *hashtable, SeqGraph *sg_hashtable, FILE* file, int fileFormat) :
        hashtable(hashtable), sg_hashtable(sg_hashtable), file(file), fileFormat(fileFormat) {}
    void operator()(KmerNode *node) {
        if (isNodeMerged<KmerNode>(node))
            return; // already merged this node
#ifdef SMALL_NODES
        unsigned left_count = 0, right_count = 0;
        node->getNeighborCounts(&left_count, &right_count);
#else
        unsigned bitconnections = trav->bitconnections;
        unsigned left_count = CONNECTION_COUNT_MAP[LEFT(bitconnections)];
        unsigned right_count = CONNECTION_COUNT_MAP[RIGHT(bitconnections)];
#endif
        if (left_count != 1) {
            build_node(node, GO_RIGHT, sg_hashtable, file, fileFormat);
        } else if (right_count != 1) {
            build_node(node, GO_LEFT, sg_hashtable, file, fileFormat);
        }
    }
  private:
    KmerGraph *hashtable;
    SeqGraph *sg_hashtable;
    FILE *file;
    int fileFormat;
};

struct KmerConcatPass2 {
    KmerConcatPass2(KmerGraph *hashtable, SeqGraph *sg_hashtable, FILE* file, int fileFormat) :
        hashtable(hashtable), sg_hashtable(sg_hashtable), file(file), fileFormat(fileFormat) {}
    void operator()(KmerNode *node) {
        if (isNodeMerged<KmerNode>(node))
            return;  // already merged this node
        if (isNodeMergingRight<KmerNode>(node))
            build_node(node, GO_RIGHT, sg_hashtable, file, fileFormat);
        else if (isNodeMergingLeft<KmerNode>(node))
            build_node(node, GO_LEFT, sg_hashtable, file, fileFormat);
     }
  private:
    KmerGraph *hashtable;
    SeqGraph *sg_hashtable;
    FILE *file;
    int fileFormat;
};

struct KmerConcatPass3 {
    KmerConcatPass3(KmerGraph *hashtable, SeqGraph *sg_hashtable, FILE* file, int fileFormat) :
        hashtable(hashtable), sg_hashtable(sg_hashtable), file(file), fileFormat(fileFormat) {}
    void operator()(KmerNode *node) {
        if (!isNodeMerged<KmerNode>(node)) {
            build_node(node, GO_RIGHT, sg_hashtable, file, fileFormat); // XXX: just create a seq node for this kmer
        }
     }
  private:
    KmerGraph *hashtable;
    SeqGraph *sg_hashtable;
    FILE *file;
    int fileFormat;
};


// generate sequence nodes by concatenating kmer nodes
static void
make_seqgraph(KmerGraph *kg_hashtable, SeqGraph *sg_hashtable, FILE* file, int fileFormat)
{
    p__SG_NODES_BUILT = 0;

    printf("======== CONCATENATION: KMER --> SEQUENCE GRAPH CONSTRUCTION PASS 1 ========\n");
    {
        KmerConcatPass1 f1(kg_hashtable, sg_hashtable, file, fileFormat);
        kg_for_each(kg_hashtable, f1);
    }

    printf("======== CONCATENATION: KMER --> SEQUENCE GRAPH CONSTRUCTION PASS 2 ========\n");
    {
        KmerConcatPass2 f2(kg_hashtable, sg_hashtable, file, fileFormat);
        kg_for_each(kg_hashtable, f2);
    }

	// corner case: catch nodes with exactly one arc in each direction, both being boundary arcs
	if( g__PSEUDO_NODES_PRESENT ) {
		printf("======== CONCATENATION: KMER --> SEQUENCE GRAPH CONSTRUCTION PASS 3 ========\n");
        {
            KmerConcatPass3 f3(kg_hashtable, sg_hashtable, file, fileFormat);
            kg_for_each(kg_hashtable, f3);
        }
	}

    /*XXX printf("======== CONCATENATION: KMER --> SEQUENCE GRAPH CONSTRUCTION DEBUG LEFTOVERS ========\n");
    for (int i = 0 ; i < HASH_BUCKETS ; ++ i) {
        for (KmerNode *trav = kg_hashtable[i] ; trav != NULL ; trav = trav->next) {
			if (!isNodeMerged<KmerNode>(trav)) {
				kmer_print_node(trav->kmer, g__FULLKMER_LENGTH);
			}
        }
    }*/

/*
    // TODO: (dna can be circular) concatenate cycles ...
#ifdef VERIFY
    printf("======== CONCATENATION: KMER --> SEQUENCE GRAPH CONSTRUCTION VERIFY ========\n");
    for (int i = 0 ; i < HASH_BUCKETS ; ++ i) {
        for (kg_node_t *trav = kg_hashtable[i] ; trav != NULL ; trav = trav->next) {
//            assert( isNodeMerged<kg_node_t>(trav) || isCycle<kg_node_t>(trav) ); // FIXME TODO
        }
    }
#endif // VERIFY
*/
/*
#ifdef VERIFY
    // SHOULD HAVE GRABBED ALL OF THE NODES!!  (except for small cycles in weird circumstances)
    // FIXME: does this pass actually do anything!?
    printf("==================== SEQ GRAPH CONSTRUCTION PASS 3 ====================\n");
    for (int i = 0 ; i < HASH_BUCKETS ; ++ i) {
        for (kg_node_t *trav = kg_hashtable[i] ; trav != NULL ; trav = trav->next) {
            if (trav->dummy_value != MERGED_VALUE) {
                bool moving_right = true;  // look for the cycle to the right.
                kg_node_t *trav2 = trav;
                do {
                    counter_t *counters = moving_right ? trav2->right_count : trav2->left_count;
                    int valid_dir = valid_single_successor(counters);
                    assert((valid_dir != MULTIPLE_SUCCESSORS) && (valid_dir != NO_SUCCESSORS));
                    bool sense_changed = counters[valid_dir] < 0;
                    assert(trav2 != NULL);
#ifdef SMALL_NODES
                    trav2 = getnext_node(HASHTABLE, trav2->kmer, valid_dir, moving_right, g__FULLKMER_LENGTH);
                    assert( g__PSEUDO_NODES_PRESENT || trav2 != NULL );
#else
                    trav2 = (moving_right ? trav2->right : trav2->left)[valid_dir];
#endif
                    moving_right ^= sense_changed;
                } while (trav2 != trav);
            }
        }
    }
#endif
*/
}

// generate sequence component by concatenating a kmer component
void sg_build_component(KmerGraph *kg_hashtable, std::deque<KmerNode*> &nodes, SeqGraph *sg_hashtable)
{
    // FIXME: ugly flags hack
    for(std::deque<KmerNode*>::iterator it=nodes.begin(); it != nodes.end(); ++it) {
        resetNodeMergingFlags<KmerNode>(*it);
    }

    // NOTE: this code should mirror the code in 'make_seqgraph()'
    //printf("======== COMPONENT CONCATENATION: KMER --> SEQUENCE GRAPH CONSTRUCTION PASS 1 ========\n");
    {
        KmerConcatPass1 f1(kg_hashtable, sg_hashtable, NULL, NOFORMAT);
        for(std::deque<KmerNode*>::iterator it=nodes.begin(); it != nodes.end(); ++it) {
            f1(*it);
        }
    }

    //printf("======== COMPONENT CONCATENATION: KMER --> SEQUENCE GRAPH CONSTRUCTION PASS 2 ========\n");
    {
        KmerConcatPass2 f2(kg_hashtable, sg_hashtable, NULL, NOFORMAT);
        for(std::deque<KmerNode*>::iterator it=nodes.begin(); it != nodes.end(); ++it) {
            f2(*it);
        }
    }

	// corner case: catch nodes with exactly one arc in each direction, both being boundary arcs
	if( g__PSEUDO_NODES_PRESENT ) {
		//printf("======== COMPONENT CONCATENATION: KMER --> SEQUENCE GRAPH CONSTRUCTION PASS 3 ========\n");
        {
            KmerConcatPass3 f3(kg_hashtable, sg_hashtable, NULL, NOFORMAT);
            for(std::deque<KmerNode*>::iterator it=nodes.begin(); it != nodes.end(); ++it) {
                f3(*it);
            }
        }
	}

#ifdef VELOUR_TBB
    size_t currentLiveMemory = g__NODE_ALLOCATORS->GetLiveMemory();
    size_t oldPeak = ATOMIC_LOAD(g__peakLiveMemory);
    if (oldPeak < currentLiveMemory) { // doesn't spin -- if C&S failed then assume new value is more correct
        ATOMIZE(g__peakLiveMemory).compare_and_swap(currentLiveMemory, oldPeak);
    }
#else
    g__peakLiveMemory = max(g__peakLiveMemory, g__NODE_ALLOCATORS->GetLiveMemory());
#endif

    // "delete" the kmer nodes
    for(std::deque<KmerNode*>::iterator it=nodes.begin(); it != nodes.end(); ++it) {
#ifdef VELOUR_TBB
        atomic_setNodeDead<KmerNode>(*it); // explicitly atomically set node dead, despite deallocation sets dead too
#else
        setNodeDead<KmerNode>(*it);
#endif // VELOUR_TBB

        //kg_hashtable->removeNode(*it); // TODO OPT: need to fix kgraph parallel iterator if enabled

        // XXX XXX: node contents must not be broken by the deallocation since we didn't remove the node!!!

        // XXX XXX: deallocating the last allocated kmernode will decrement allocation pointer, should be okay
        //            since we don't allocate anymore kmer nodes
 
#ifdef VELOUR_TBB
        tls_kmernode_allocator->DeallocateNodeMemory(*it); // FIXME: can't deallocate if not removed above?
#else
        g__KMERNODE_ALLOCATOR->DeallocateNodeMemory(*it); // FIXME: can't deallocate if not removed above?
#endif
    }
}

static void
sg_mark_all_neighbors_as_merge_starts(SeqGraph *graph, SeqNode *sg_node, int direction)
{
    bool moving_right = (direction == GO_RIGHT);
    counter_t *counters = moving_right ? sg_node->right_count : sg_node->left_count;
 
#ifndef SMALL_NODES
    SeqNode **pointers = moving_right ? sg_node->right : sg_node->left;
#endif

    for (int i = 0 ; i < 4 ; ++i ) {
        counter_t count = counters[i];
        if (count == 0) continue;
        bool sense_changed; // = (count < 0);
        color_t *colors = moving_right ? sg_node->right_color : sg_node->left_color;
        SeqNode *neighbor = NULL;
        if (colors[i] == 0) {
#ifdef SMALL_NODES
            neighbor = graph->findNextNode(sg_node, i, moving_right, &sense_changed);
#else
            neighbor = pointers[i];
#endif
        } else {
            continue;
        }
        assert( neighbor != NULL );
        //if (isNodeMerged<SeqNode>(neighbor)) continue; // commented out for performance, FIXME verify benefit
        if (moving_right ^ sense_changed) {
            setNodeMergingRight<SeqNode>(neighbor);
        } else {
            setNodeMergingLeft<SeqNode>(neighbor);
        }
    }
}

static inline void
sg_copy_counters(SeqNode *next_node, SeqNode *root_sg_node, int direction, int root_direction)
{
	assert(next_node != root_sg_node);

	four_counter_t *sg_side = direction == GO_RIGHT ? &next_node->right_side : &next_node->left_side;
	four_color_t *sg_cside = direction == GO_RIGHT ? &next_node->right_side_colors : &next_node->left_side_colors;
	four_counter_t *root_side = root_direction == GO_RIGHT ? &root_sg_node->right_side : &root_sg_node->left_side;
	four_color_t *root_cside = root_direction == GO_RIGHT ? &root_sg_node->right_side_colors : &root_sg_node->left_side_colors;
	//assert( *root_side == 0 );

	if (direction == root_direction) {
		*root_side = *sg_side;
		*root_cside = *sg_cside;
	} else { // negate counters and complement, as sense is reversed
		counter_t *sg_counters = static_cast<counter_t *>(static_cast<void *>(sg_side));
		counter_t *root_counters = static_cast<counter_t *>(static_cast<void *>(root_side));
		Color *sg_colors = static_cast<Color *>(static_cast<void *>(sg_cside));
		Color *root_colors = static_cast<Color *>(static_cast<void *>(root_cside));
		for (int i=0; i < 4; ++i) {
			root_counters[i] = -sg_counters[COMPLEMENT(i)];
            root_colors[i] = sg_colors[COMPLEMENT(i)];
		}
	}

    // TODO: bitconnections?
}

static SeqNode *
sg_extend_contig(SeqGraph *__arg_graph, SeqNode *__arg_root_sg_node, int __arg_root_direction)
{
    SeqGraph *graph = __arg_graph;
    SeqNode *root_sg_node = __arg_root_sg_node;
    int root_direction = __arg_root_direction;

    unsigned iteration_counter = 1;

    while(1) {  // loop instead of recursion

    bool moving_right = (root_direction == GO_RIGHT);
    counter_t *counters = moving_right ? root_sg_node->right_count : root_sg_node->left_count;

#ifdef SMALL_NODES
    int valid_dir = valid_single_successor(counters);
    if (valid_dir == MULTIPLE_SUCCESSORS) {
        sg_mark_all_neighbors_as_merge_starts(graph, root_sg_node, root_direction);
        return root_sg_node;
    }
    if (valid_dir == NO_SUCCESSORS) {
        return root_sg_node;
    }
    assert( abs(counters[valid_dir]) > 0 );

    // single direction
    assert((valid_dir >= 0) && (valid_dir < 4));
    bool sense_changed;
    color_t *colors = moving_right ? root_sg_node->right_color : root_sg_node->left_color;
    SeqNode *next = NULL;
    if (colors[valid_dir] == 0) {
        next = graph->findNextNode(root_sg_node, valid_dir, moving_right, &sense_changed);
    } else {
        return root_sg_node;
    }
    assert( next != NULL );
#else
	/*
    // TODO: need to setup the pointers for sg_node somewhere... 
    unsigned bitconnections = sg_node->bitconnections;
    unsigned connex = moving_right ? RIGHT(bitconnections) : LEFT(bitconnections);
    unsigned connex_count = CONNECTION_COUNT_MAP[connex];
    if (connex_count != 1) { // can't extend
        if (connex_count != 0) {
            sg_copy_counters(sg_node, root_sg_node, direction);
            sg_mark_all_neighbors_as_merge_starts(sg_node, direction);
        }
        return root_sg_node;
    }

    // single direction
    SeqNode **pointers = moving_right ? sg_node->right : sg_node->left;
    int valid_dir = valid_single_successor(counters);
    assert((valid_dir >= 0) && (valid_dir < 4));
    SeqNode *next = pointers[valid_dir];
    bool sense_changed = counters[valid_dir] < 0;
	*/
#endif

    if (isNodeMerged(next)) { // can't extend
        return root_sg_node;
    }

    bool next_moving_right = moving_right ^ sense_changed;
    counter_t *next_back_counters = next_moving_right ? next->left_count : next->right_count;
    int next_back_valid = valid_single_successor(next_back_counters);
    if (next_back_valid == MULTIPLE_SUCCESSORS) {
        return root_sg_node;   // can't extend
    }

    // check for: - single copy clipping marking, which inhibits concatenation
    //            - start at base of Y, conscious of possible erosion
    //            - copy count disagreement, for proper single-copy erosion behavior
    //
    // XXX alternative implementation: traverse and check for long _tip_ that shouldn't clip -- only then can ignore/remove clip mark
    if (g__PSEUDO_NODES_PRESENT) {
        assert( isNodeWeight_AllNonSingleCopy(root_sg_node) || isNodeWeight_AllSingleCopy(root_sg_node) || isNodeWeight_None(root_sg_node) );

        // single copy clipping marking
        if ( abs(counters[valid_dir]) == CLIP_SINGLECOPY_COUNTER_VALUE ||
             abs(next_back_counters[next_back_valid]) == CLIP_SINGLECOPY_COUNTER_VALUE ) {
            sg_mark_all_neighbors_as_merge_starts(graph, root_sg_node, root_direction); // OPT: just merge-start the one arc?
            return root_sg_node;
        }

        counter_t edge_weight = abs(counters[valid_dir]);
        if (edge_weight == CLIP_SINGLECOPY_COUNTER_VALUE) edge_weight = 1;
        assert( edge_weight > 0 );

        counter_t *next_forward_counters = next_moving_right ? next->right_count : next->left_count;
        counter_t next_forward_edges_weight = get_counter_sum(next_forward_counters); // XXX: or max?

        counter_t *root_back_counters = (root_direction == GO_RIGHT) ? root_sg_node->left_count : root_sg_node->right_count;
        bool is_root_back_multiple_successors = (valid_single_successor(root_back_counters) == MULTIPLE_SUCCESSORS);

        // base of Y: erosion unclear
        if (edge_weight == 1) {
            if (valid_single_successor(next_forward_counters) == MULTIPLE_SUCCESSORS) {
                /*bool next_forward_single_copy = true;
                for (int i=0; i < 4; ++i) {
                    counter_t value = abs(next_forward_counters[i]);
                    if (value != 0 && value != 1 && value != CLIP_SINGLECOPY_COUNTER_VALUE) {
                        next_forward_single_copy = false;
                    }
                }
                if (next_forward_single_copy) { // all single copy edges*/
                    sg_mark_all_neighbors_as_merge_starts(graph, root_sg_node, root_direction); // OPT: just merge-start the one arc?
                    return root_sg_node;
                //}
            }
            //if (root_sg_node->min_edge_weight == 0) { // optimization: a non-zero merged node will always be ok
                if (valid_single_successor(root_back_counters) == MULTIPLE_SUCCESSORS) {
                    /*bool root_back_single_copy = true;
                    for (int i=0; i < 4; ++i) {
                        counter_t value = abs(root_back_counters[i]);
                        if (value != 0 && value != 1 && value != CLIP_SINGLECOPY_COUNTER_VALUE) {
                            root_back_single_copy = false;
                        }
                    }
                    if (root_back_single_copy) { // all single copy backward edges*/
                        sg_mark_all_neighbors_as_merge_starts(graph, root_sg_node, root_direction); // OPT: just merge-start the one arc?
                        return root_sg_node;
                    //}
                }
            //}
        }
        counter_t sg_back_edge_weight = get_counter_sum(root_back_counters); // XXX: or max?

        // known: sg_back is not single copy multiple back successor, thus can forward merge any copy value

        // copy count disagreement: single vs multi
        if ((isNodeWeight_None(root_sg_node) && !is_root_back_multiple_successors && edge_weight == 1 && (sg_back_edge_weight > 1 || next_forward_edges_weight > 1)) ||
            (isNodeWeight_HasSingleCopy(root_sg_node) && (edge_weight > 1 || next_forward_edges_weight > 1)) ||
            (isNodeWeight_AllNonSingleCopy(root_sg_node) && edge_weight == 1) ) {
            sg_mark_all_neighbors_as_merge_starts(graph, root_sg_node, root_direction); // OPT: just merge-start the one arc?
            return root_sg_node;
        }
    }

    if (!root_sg_node->sequence.CanAddBases(next->sequence.get_length() - g__FULLKMER_LENGTH + 1)) { // impossible to concatenate the next node
        sg_mark_all_neighbors_as_merge_starts(graph, root_sg_node, root_direction); // OPT: just merge-start the one arc?
        return root_sg_node;
    }

    // otherwise, we can extend.  Need to add next sequence minus overlap, and then call recursively.
	//   1. remove hashtable head OR tail entry for root node's defunct side
	//   2. concatenate next node
	//   3. remove hashtable head AND tail entries for next node
	//   4. insert hashtable head OR tail entry for root node's new side

	//if (debugprint) { printf("\nEXTENDING!\n"); }
	if (root_direction == GO_RIGHT) {
		graph->removeNodeTail(root_sg_node);
	} else {
		graph->removeNodeHead(root_sg_node);
	}

    // NOTE: we update internal edge weight _before_ copying the new counters!
    counter_t edge_weight = abs(counters[valid_dir]);
    if (edge_weight == CLIP_SINGLECOPY_COUNTER_VALUE) edge_weight = 1;
    assert( edge_weight > 0 );
    updateNodeWeight(root_sg_node, edge_weight);

    // update kmer occurrences count
    root_sg_node->kmer_occurrences += next->kmer_occurrences;

    SeqNode *oldnode = root_sg_node;
    root_sg_node->sequence.Concatenate(graph, &root_sg_node, &next->sequence, root_direction, sense_changed, g__FULLKMER_LENGTH); // NOTE: can cause node growth
    if (root_sg_node != oldnode) { // node was grown
        if (p__GROWNODELIST != NULL) {
            p__GROWNODELIST->push_back(root_sg_node);
        }
    }

    sg_copy_counters(next, root_sg_node, next_moving_right, moving_right);

    graph->removeNodeUnsafe(next); // the graph iterator is never pointing at 'next'

	if (root_direction == GO_RIGHT) {
		graph->insertNodeTail(root_sg_node);
	} else {
		graph->insertNodeHead(root_sg_node);
	}

	++ p__SG_NODES_CONCAT;

    if (g__PSEUDO_NODES_PRESENT) {
        assert( isNodeWeight_AllNonSingleCopy(root_sg_node) || isNodeWeight_AllSingleCopy(root_sg_node) || isNodeWeight_None(root_sg_node) );
    }

	// free "next" node (already disconnected from seq graph)
#ifdef VELOUR_TBB
    tls_seqnode_allocator->DeallocateNodeMemory(next);
#else
    g__SEQNODE_ALLOCATOR->DeallocateNodeMemory(next);
#endif // VELOUR_TBB

	// extend contig to next node
    
    // setup variables for next iteration of loop, i.e. recursive function call
    //   WAS: sg_extend_contig(graph, root_sg_node, root_direction);
    // graph = graph;
    // root_sg_node = root_sg_node;
    // root_direction = root_direction;
    
    // track the "recursion depth"
    ++ iteration_counter;
    }
    assert( false );
}

static SeqNode *
concat_node(SeqGraph *graph, SeqNode *sg_node, int direction)
{
    assert( !isNodeMerged<SeqNode>(sg_node) );
    setNodeMerged<SeqNode>(sg_node);

    // mark the other direction for collection on the next pass
    sg_mark_all_neighbors_as_merge_starts(graph, sg_node, !direction);

	// try to extend this node
    return sg_extend_contig(graph, sg_node, direction);
}

struct SeqConcatPass1 {
    SeqConcatPass1(SeqGraph *graph) : graph(graph) {}
    SeqNode * operator()(SeqNode *node) {
        if (isNodeMerged<SeqNode>(node))
            return node; // already merged this node
#ifdef SMALL_NODES
        unsigned left_count = 0, right_count = 0;
        node->getNeighborCounts(&left_count, &right_count);
#else
        unsigned bitconnections = node->bitconnections;
        unsigned left_count = CONNECTION_COUNT_MAP[LEFT(bitconnections)];
        unsigned right_count = CONNECTION_COUNT_MAP[RIGHT(bitconnections)];
#endif
        if (left_count != 1) {
            return concat_node(graph, node, GO_RIGHT);
        } else if (right_count != 1) {
            return concat_node(graph, node, GO_LEFT);
        }
        return node;
    }
  private:
    SeqGraph *graph;
};

struct SeqConcatPass2 {
    SeqConcatPass2(SeqGraph *graph) : graph(graph) {}
    SeqNode * operator()(SeqNode *node) {
        if (isNodeMerged<SeqNode>(node))
            return node;  // already merged this node
        if (isNodeMergingRight<SeqNode>(node))
            return concat_node(graph, node, GO_RIGHT);
        else if (isNodeMergingLeft<SeqNode>(node))
            return concat_node(graph, node, GO_LEFT);
        else
            return node;
    }
  private:
    SeqGraph *graph;
};


static void
concat_seqgraph(SeqGraph *graph, bool silent)
{
    p__SG_NODES_CONCAT = 0;

    if (!silent) {
        printf("======== CONCATENATION: SEQUENCE GRAPH PASS 1 ========\n");
    }

    {
        SeqConcatPass1 f1(graph);
        sg_for_each_grow(graph, f1);
    }

    if (!silent) {
        printf("======== CONCATENATION: SEQUENCE GRAPH PASS 2 ========\n");
    }

    {
        SeqConcatPass2 f2(graph);
        sg_for_each_grow(graph, f2);
    }

    // XXX: leftovers?

	//sg_verify_graph();
/*
    // TODO: (dna can be circular) concatenate cycles ...
#ifdef VERIFY
    printf("======== CONCATENATION: SEQUENCE GRAPH VERIFY ========\n");
    for (int i = 0 ; i < HASH_BUCKETS ; ++ i) {
        for (SeqNode *trav = sg_hashtable[i] ; trav != NULL ; trav = trav->head_next) {
            assert( isNodeMerged<SeqNode>(trav) || isCycle<SeqNode>(trav) );
        }
    }
#endif // VERIFY
*/
}
    
void sg_concatenate(SeqGraph *graph, bool silent)
{
    graph->resetFlags();
	concat_seqgraph(graph,silent);

    if (!silent) {
        printf("%"PRIuPTR" nodes concatenated.\n", p__SG_NODES_CONCAT);
    }
}

void sg_nodelist_concatenate(SeqGraph *graph, bool silent, flow_nodelist_t *nodelist)
{
    for (flow_nodelist_t::iterator it = nodelist->begin(); it != nodelist->end(); ++it) {
        resetNodeMergingFlags<SeqNode>(*it);
    }

    p__SG_NODES_CONCAT = 0;

    if (!silent) {
        printf("======== CONCATENATION: SEQUENCE GRAPH PASS 1 ========\n");
    }

    {
        flow_nodelist_t grown_nodes;
        p__GROWNODELIST = &grown_nodes;

        SeqConcatPass1 f1(graph);
        for (flow_nodelist_t::iterator it = nodelist->begin(); it != nodelist->end(); ++it) {
            SeqNode *node = *it;
            if (!isNodeDead(node)) { f1(node); }
        }

        // then add the newly grown nodes
#ifdef VELOUR_TBB
        for (flow_nodelist_t::iterator it = grown_nodes.begin(); it != grown_nodes.end(); ++it) {
            nodelist->push_back(*it);
        }
#else
        nodelist->insert(nodelist->end(), grown_nodes.begin(), grown_nodes.end());
#endif
        p__GROWNODELIST = NULL; // important
    }

    if (!silent) {
        printf("======== CONCATENATION: SEQUENCE GRAPH PASS 2 ========\n");
    }

    {
        flow_nodelist_t grown_nodes;
        p__GROWNODELIST = &grown_nodes;

        SeqConcatPass2 f2(graph);
        for (flow_nodelist_t::iterator it = nodelist->begin(); it != nodelist->end(); ++it) {
            SeqNode *node = *it;
            if (!isNodeDead(node)) { f2(node); }
        }

        // then add the newly grown nodes
#ifdef VELOUR_TBB
        for (flow_nodelist_t::iterator it = grown_nodes.begin(); it != grown_nodes.end(); ++it) {
            nodelist->push_back(*it);
        }
#else
        nodelist->insert(nodelist->end(), grown_nodes.begin(), grown_nodes.end());
#endif
        p__GROWNODELIST = NULL; // important
    }
    
    if (!silent) {
        printf("%"PRIuPTR" nodes concatenated.\n", p__SG_NODES_CONCAT);
    }
}


//
// exported functions
//

void sg_build(SeqGraph *sg_hashtable, KmerGraph *kg_hashtable)
{
    make_seqgraph(kg_hashtable, sg_hashtable, NULL, NOFORMAT);
    printf("%"PRIuPTR" sequence nodes built from kmer graph\n", p__SG_NODES_BUILT);
}

void sg_dump_quilt_from_kmergraph(KmerGraph *kg_hashtable, const char *filename)
{
	size_t retval;
    int filedes = open(filename, O_CREAT | O_TRUNC | O_WRONLY, 0666);
    FILE* file = fdopen(filedes, "w");
    // FIXME: performance -- setbuffer() a larger buffer?

    // QUILT File Header: node count, kmer length
    uintptr_t dummy = 0;
    retval = fwrite(&dummy, sizeof(uintptr_t), 1, file);
	assert(retval == 1 && "fwrite() error." );
    retval = fwrite(&g__FULLKMER_LENGTH, sizeof(unsigned), 1, file);
	assert(retval == 1 && "fwrite() error." );

    p__SG_NODES_EMITTED = 0;
    make_seqgraph(kg_hashtable, NULL, file, QUILT);
	assert( p__SG_NODES_EMITTED == p__SG_NODES_BUILT );
    printf("%"PRIuPTR" sequence nodes built from kmer graph\n", p__SG_NODES_BUILT);
    fclose(file);

    // update file header's node count
    file = fopen(filename,"r+");
    retval = fwrite(&p__SG_NODES_EMITTED, sizeof(p__SG_NODES_EMITTED), 1, file);
	assert(retval == 1 && "fwrite() error." );
    fclose(file);
}

void sg_dump_pregraph_from_kmergraph(KmerGraph *kg_hashtable, const char *filename)
{
    int filedes = open(filename, O_CREAT | O_TRUNC | O_WRONLY, 0666);
    FILE* file = fdopen(filedes, "w");

    // PREGRAPH File Header: node count, read count, boundary count, kmer length, double stranded
    //   note: using hex so can fixup node count later
    fprintf(file, "0x%08lx\t%" PRIu64 "\t%d\t%d\t%hd\n", 0UL, g__READ_COUNT, 0, g__FULLKMER_LENGTH, 1);

    p__SG_NODES_EMITTED = 0;
    make_seqgraph(kg_hashtable, NULL, file, PREGRAPH);
	assert( p__SG_NODES_EMITTED == p__SG_NODES_BUILT );
    printf("%"PRIuPTR" sequence nodes built from kmer graph\n", p__SG_NODES_BUILT);
    fclose(file);

    // update file header's node count
    file = fopen(filename,"r+");
    fprintf(file, "0x%08"PRIxPTR"", p__SG_NODES_EMITTED);
    fclose(file);
}

void sg_load_bucket(SeqGraph *sg_hashtable, const char *filename)
{
    int filedes = open(filename, O_RDONLY);
    FILE* file = fdopen(filedes, "r");
    // FIXME: performance -- setbuffer() a larger buffer?

    printf("%" PRIuPTR " sequence nodes present before bucket file %s\n", sg_hashtable->node_count, filename);

    p__SG_NODES_LOADED = 0;

    size_t component_size;
    while (fread(&component_size, sizeof(component_size), 1, file) == 1) {
        uintptr_t total_nodes;
        if (fread(&total_nodes, sizeof(total_nodes), 1, file) != 1) {
            fprintf(stderr, "ERROR: failed to fread()\n");
            exit(EXIT_FAILURE);
        }
        for (uintptr_t component_node_count=0; component_node_count < total_nodes; ++component_node_count) {
            SeqNode *node = NULL;
            node = load_seq_node(file, BUCKET);
            assert( node != NULL && "Component node count must be wrong?" );
            sg_hashtable->insertNodeAndUpdateColors(node);
        }
    }

    assert( feof(file) && "fread() failed but not end-of-file!?" );
    fclose(file);

    printf("%"PRIuPTR" sequence nodes loaded from bucket file %s\n", p__SG_NODES_LOADED, filename);
    printf("%"PRIuPTR" sequence nodes present after bucket file %s\n", sg_hashtable->node_count, filename);
}

uintptr_t sg_load_stream_bucket(SeqGraph *sgraph, const uintptr_t node_limit, FILE *file)
{
    uintptr_t nodes_loaded = 0;

    size_t component_size;
    while (nodes_loaded < node_limit && fread(&component_size, sizeof(component_size), 1, file) == 1) {
        uintptr_t total_nodes;
        if (fread(&total_nodes, sizeof(total_nodes), 1, file) != 1) {
            fprintf(stderr, "ERROR: failed to fread()\n");
            exit(EXIT_FAILURE);
        }
        uintptr_t component_node_count;
        for (component_node_count=0; component_node_count < total_nodes; ++component_node_count) {
            SeqNode *node = load_seq_node(file, BUCKET);
            assert( node != NULL && "Component node count must be wrong?" );
            sgraph->insertNodeAndUpdateColors(node);
        }
        nodes_loaded += component_node_count;
    }
    printf("%"PRIuPTR" nodes loaded from stream bucket.\n", nodes_loaded);

    return nodes_loaded;
}

uintptr_t sg_load_stream_component(SeqGraph *sgraph, FILE *file)
{
    uintptr_t nodes_loaded = 0;

    size_t component_size;
    if (fread(&component_size, sizeof(component_size), 1, file) == 1) {
        uintptr_t total_nodes;
        if (fread(&total_nodes, sizeof(total_nodes), 1, file) != 1) {
            fprintf(stderr, "ERROR: failed to fread()\n");
            exit(EXIT_FAILURE);
        }
        for (nodes_loaded=0; nodes_loaded < total_nodes; ++nodes_loaded) {
            SeqNode *node = load_seq_node(file, BUCKET);
            assert( node != NULL && "Component node count must be wrong?" );
            sgraph->insertNode(node);
        }
    }
    //printf("%"PRIuPTR" component nodes loaded from stream bucket.\n", nodes_loaded);

    return nodes_loaded;
}

uintptr_t sg_load_mmap_stream_component(char *mmap_start, size_t *amount, SeqNode **chain)
{
    uintptr_t nodes_loaded = 0;
    off_t offset = 0;

    size_t component_size;
    memcpy(&component_size, (mmap_start + offset), sizeof(component_size));
    offset += sizeof(component_size);

    uintptr_t total_nodes;
    memcpy(&total_nodes, (mmap_start + offset), sizeof(total_nodes));
    offset += sizeof(total_nodes);

    for (nodes_loaded=0; nodes_loaded < total_nodes; ++nodes_loaded) {
        size_t nodesize;
        SeqNode *node = load_seq_node_mmap((mmap_start+offset), nodesize);
        offset += nodesize;
        assert( node != NULL && "Component node count must be wrong?" );
        node->head_next = *chain;
        *chain = node;
    }

    assert( offset == component_size );

    *amount = offset;

    return nodes_loaded;
}

void sg_load_quilt(SeqGraph *sg_hashtable, const char *filename)
{
	size_t retval;
    int filedes = open(filename, O_RDONLY);
    FILE* file = fdopen(filedes, "r");
    // FIXME: performance -- setbuffer() a larger buffer?
    
    printf("%" PRIuPTR " sequence nodes present before quilt file %s\n", sg_hashtable->node_count, filename);

    // QUILT File Header: node count, kmer length
    uintptr_t file_nodecount;
    unsigned file_kmerlength;
    retval = fread(&file_nodecount, sizeof(sg_hashtable->node_count), 1, file);
	assert(retval == 1 && "fread() error." );
    retval = fread(&file_kmerlength, sizeof(unsigned), 1, file);
	assert(retval == 1 && "fread() error." );
    assert( file_kmerlength == g__FULLKMER_LENGTH );

    p__SG_NODES_LOADED = 0;
	SeqNode *node = NULL;
    while ((node = load_seq_node(file, QUILT)) != NULL) {
		sg_hashtable->insertNodeAndUpdateColors(node);
    }
    fclose(file);
    assert( p__SG_NODES_LOADED == file_nodecount );

    printf("%"PRIuPTR" sequence nodes loaded from quilt file %s\n", p__SG_NODES_LOADED, filename);
    printf("%" PRIuPTR " sequence nodes present after quilt file %s\n", sg_hashtable->node_count, filename);
}

struct lambda_dump_quilt {
    lambda_dump_quilt(SeqGraph *graph, FILE* file) : graph(graph), file(file) {}
    void operator()(SeqNode *node) { emit_seq_node(node, file, QUILT); }
  private:
    SeqGraph *graph;
    FILE *file;
};

void sg_dump_quilt(SeqGraph *sg_hashtable, const char *filename)
{
	size_t retval;
    int filedes = open(filename, O_CREAT | O_TRUNC | O_WRONLY, 0666);
    FILE* file = fdopen(filedes, "w");
    // FIXME: performance -- setbuffer() a larger buffer?

    // QUILT File Header: node count, kmer length
    retval = fwrite(&sg_hashtable->node_count, sizeof(sg_hashtable->node_count), 1, file); // make this a dummy write
	assert(retval == 1 && "fwrite() error." );
    retval = fwrite(&g__FULLKMER_LENGTH, sizeof(unsigned), 1, file);
	assert(retval == 1 && "fwrite() error." );

    p__SG_NODES_EMITTED = 0;
    {
        lambda_dump_quilt f(sg_hashtable, file);
        sg_for_each(sg_hashtable, f);
    }
    fclose(file);
    printf("%"PRIuPTR" sequence nodes dumped to quilt file %s\n", p__SG_NODES_EMITTED, filename);

	assert( p__SG_NODES_EMITTED == sg_hashtable->node_count );

/*	// update file header's node count
    file = fopen(filename,"r+");
    retval = fwrite(&p__SG_NODES_EMITTED, sizeof(unsigned), 1, file);
	assert(retval == 1 && "fwrite() error." );
    fclose(file); */
}

struct lambda_dump_pregraph {
    lambda_dump_pregraph(SeqGraph *graph, FILE* file) : graph(graph), file(file) {}
    void operator()(SeqNode *node) { emit_seq_node(node, file, PREGRAPH); }
  private:
    SeqGraph *graph;
    FILE *file;
};

void sg_dump_pregraph(SeqGraph *sg_hashtable, const char *filename)
{
    int filedes = open(filename, O_CREAT | O_TRUNC | O_WRONLY, 0666);
    FILE* file = fdopen(filedes, "w");
    // FIXME: performance -- setbuffer() a larger buffer?

    // PREGRAPH File Header: node count, read count, kmer length
    //   note: using hex so can fixup node count later
    //fprintf(file, "%" PRIuPTR "\t%" PRIu64 "\t%d\t%d\t%hd\n", sg_hashtable->node_count, g__READ_COUNT, 0, g__FULLKMER_LENGTH, 1);
    fprintf(file, "%" PRIuPTR "\t%" PRIu64 "\t%d\t%hd\n", sg_hashtable->node_count, g__READ_COUNT, g__FULLKMER_LENGTH, 1);

    p__SG_NODES_EMITTED = 0;
    {
        lambda_dump_pregraph f(sg_hashtable, file);
        sg_for_each(sg_hashtable, f);
    }
    fclose(file);
    printf("%"PRIuPTR" sequence nodes dumped to pregraph file %s\n", p__SG_NODES_EMITTED, filename);
	assert( p__SG_NODES_EMITTED == sg_hashtable->node_count);
}
