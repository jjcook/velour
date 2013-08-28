//
// This file is part of Velour and is distributed under the University of
// Illinois Open Source License.  See LICENSE.txt for details.
//

//
// sg_bubblePopping.cpp
//
//   sequence graph bubble popping -- only pops single copy bubbles
//
// TODO OPT: avoid loops when computing consensus!
//

#include "types.h"

//
// CONFIGURABLE PARAMETERS
//
//   TODO: set these based on k-mer length?
static const unsigned MAX_BUBBLE_LENGTH = 200;
static const unsigned MAX_BUBBLE_GAPS = 3;
static   const double MAX_BUBBLE_DIVERGENCE = 0.2;
//static   const double MAX_BUBBLE_COVERAGE = 1;


// internal statistics variables // TODO: need to reset some of these between iterations
static uintptr_t p__BUBBLES_POPPED = 0;
static uintptr_t p__BUBBLES_CONSIDERED = 0;
static uintptr_t p__BUBBLE_NODES_DEALLOCATED = 0;

//static uintptr_t p__PROFILE_NONE = 0;
//static uintptr_t p__PROFILE_SINGLECOPY = 0;

// sequence comparison
static unsigned matrix[MAX_BUBBLE_LENGTH+1][MAX_BUBBLE_LENGTH+1];
static const unsigned cost_indel = 0;
static const unsigned cost_map[4][4] = {
    {1, 0, 0, 0},
    {0, 1, 0, 0},
    {0, 0, 1, 0},
    {0, 0, 0, 1}
};

//
// decide if sequences are sufficiently similar
//
// TODO: use a banded variant for efficiency
//
// TODO: early exit if lengths are too dissimilar?
//
static bool compareSequences(Sequence *pop_sequence, Sequence *solid_sequence)
{
    assert( pop_sequence != NULL );
    assert( solid_sequence != NULL );

    uint16_t pop_length = pop_sequence->get_length();
    uint16_t solid_length = solid_sequence->get_length();

    if (pop_length > MAX_BUBBLE_LENGTH || solid_length > MAX_BUBBLE_LENGTH) return false;

    for (unsigned i = 0; i <= pop_length; i++)
        matrix[i][0] = 0;
    for (unsigned j = 0; j <= solid_length; j++)
        matrix[0][j] = 0;

    for (unsigned i = 1; i <= pop_length; i++) {
        for (unsigned j = 1; j <= solid_length; j++) {
            unsigned x1, x2, x3;
            x1 = matrix[i - 1][j - 1] + cost_map[(unsigned) pop_sequence->GetBase(i - 1)][(unsigned) solid_sequence->GetBase(j - 1)];
            x2 = matrix[i - 1][j] + cost_indel;
            x3 = matrix[i][j - 1] + cost_indel;
            matrix[i][j] = max(x1, max(x2, x3));
        }
    }

    unsigned maxScore = matrix[pop_length][solid_length];

    uint16_t maxLength = max(pop_length, solid_length);

    if (maxScore < maxLength - MAX_BUBBLE_GAPS)
        return false;

    if ((1 - ((double)maxScore / maxLength)) > MAX_BUBBLE_DIVERGENCE)
        return false;

    return true;
}

/*
// isolated bubbles: check one side of a node for potential bubbles to pop, and pop them
//
//                                                          //
//     _o_                                                  //
//    /   \                                                 //
// --o--o--o--   ===>   --o--                               //
//                                                          //
//                                                          //
//     _o_                  _o_                             //
//    /   \                /   \                            //
// --o--o--o--   ===>   --o--o--o--   ===>   --o--          //
//    \   /                                                 //
//     -o-                                                  //
//
//  XXX NOTE: this code will not detect bubbles comprised of non-concatenated nodes
//
static bool check_simple_pop_direction(SeqGraph *graph, SeqNode *node, bool direction)
{
    //
    // discover bubble-compliant simple paths: the poppable and solid paths
    //
    SeqNode *poppable_segment[4] = {NULL};
    SeqNode *solid_segment[4] = {NULL};

    unsigned poppable_count = 0;

    SeqNode *join[4] = {NULL};
    counter_t *join_back_counters[4] = {NULL}; // XXX OPT: store pointer to single counter, obviating join_back?
    Nucleotide join_back[4];

    counter_t *node_counters = direction ? node->right_count : node->left_count;
    color_t *node_colors = direction ? node->right_color : node->left_color;

    for (unsigned i=0 ; i < 4 ; ++i) {
        if (node_counters[i] == 0 || node_colors[i] != 0) continue;

        bool next_sense_changed;
        SeqNode *next = graph->findNextNode(node, i, direction, &next_sense_changed);
        assert( next != NULL );

        bool next_direction = direction ^ next_sense_changed;

        if (next->sequence.get_length() > MAX_BUBBLE_LENGTH) continue; // early exit: node is too long to be a bubble

        // check the backward direction for the given segment
        counter_t *next_back_counters = next_direction ? next->left_count : next->right_count;
        int next_back_arc = valid_single_successor(next_back_counters);
        assert( next_back_arc != NO_SUCCESSORS );
        if (next_back_arc == MULTIPLE_SUCCESSORS) continue;

        // check the forward direction and memoize the join node
        counter_t *next_counters = next_direction ? next->right_count : next->left_count;
        color_t *next_colors = next_direction ? next->right_color : next->left_color;

        int next_arc = valid_single_successor(next_counters);
        if (next_arc != NO_SUCCESSORS && next_arc != MULTIPLE_SUCCESSORS) {

            if (next_colors[next_arc] != 0) continue;

            bool next_join_sense_changed;
            SeqNode *next_join = graph->findNextNode(next, next_arc, next_direction, &next_join_sense_changed);
            assert( next_join != NULL );

            bool next_join_direction = next_direction ^ next_join_sense_changed;

            // compute this once we know about next_join
            Nucleotide next_head_rightmost_base = next->sequence.GetHeadKmerRightmostBase(g__FULLKMER_LENGTH);
            Nucleotide next_tail_leftmost_base = next->sequence.GetTailKmerLeftmostBase(g__FULLKMER_LENGTH);
            Nucleotide next_join_back_valid = next_direction ? next_tail_leftmost_base : next_head_rightmost_base;
            if (next_join_sense_changed) { next_join_back_valid = COMPLEMENT(next_join_back_valid); }

            //
            // record this simple path
            //
            join[i] = next_join;
            join_back_counters[i] = next_join_direction ? next_join->left_count : next_join->right_count; // the back counters
            join_back[i] = next_join_back_valid;

            if ((isNodeWeight_AllSingleCopy(next) || isNodeWeight_None(next)) && // TODO: _None is probably useless?
                    cnorm(node_counters[i]) <= MAX_BUBBLE_COVERAGE &&
                    cnorm((join_back_counters[i])[join_back[i]]) <= MAX_BUBBLE_COVERAGE) {
                poppable_segment[i] = next;
                ++ poppable_count;
            } else {
                solid_segment[i] = next;
            }
        }
    }

    if (poppable_count == 0) { // nothing to pop
        return false;
    }

    // at this point, we've identified linear paths that are short enough to be poppable

    bool modified = false;

    // NOTE: only compare segments that end at the same join node
    // NOTE: only compare poppable segments to solid segments

    // TODO: use poppable_count to exit early?
    for (unsigned i=0 ; i < 4 ; ++i) {
        for (unsigned j=(i+1) ; j < 4 ; ++j) {
            if (join[i] == join[j]) {
                unsigned pop_index;
                unsigned solid_index;
                SeqNode *the_pop_segment = NULL;
                SeqNode *the_solid_segment = NULL;
                if (poppable_segment[i] != NULL && solid_segment[j] != NULL) {
                    pop_index = i; the_pop_segment = poppable_segment[i];
                    solid_index = j; the_solid_segment = solid_segment[j];
                } else if (poppable_segment[j] != NULL && solid_segment[i] != NULL) {
                    pop_index = j; the_pop_segment = poppable_segment[j];
                    solid_index = i; the_solid_segment = solid_segment[i];
                }
                if (the_pop_segment != NULL && the_solid_segment != NULL) {
                    ++ p__BUBBLES_CONSIDERED;

                    // check if the bubble can be popped
                    bool pop = compareSequences(&the_pop_segment->sequence, &the_solid_segment->sequence);
                    if (!pop) {
                        // TODO OPT: don't RC unless the two sequences are potentially similar, e.g. length
                        Sequence_StackAllocated sequence_memory;
                        Sequence *revcomp = new (&sequence_memory) Sequence(the_pop_segment->sequence, Sequence::MAX_BASES);
                        revcomp->ReverseComplement();
                        pop = compareSequences(revcomp, &the_solid_segment->sequence);
                    }

                    if (pop) {

#ifdef VERIFY
                        graph->verify_node(node);
                        graph->verify_node(the_pop_segment);
                        graph->verify_node(join[pop_index]);
#endif

                        // TODO: TBB version
                        // TODO: TBB: make sure I'm removing / deallocating the_pop_segment node safely
                        //assert( cnorm(node_counters[pop_index]) > 0 && cnorm(node_counters[pop_index]) <= MAX_BUBBLE_COVERAGE );
                        //assert( graph->findNextNode(node, pop_index, direction) == the_pop_segment );
                        //assert( node_counters[pop_index] != 0 );
                        node_counters[pop_index] = 0;

                        //assert( cnorm((join_back_counters[pop_index])[ join_back[pop_index] ]) > 0 &&
                        //       cnorm((join_back_counters[pop_index])[ join_back[pop_index] ]) <= MAX_BUBBLE_COVERAGE );
                        //assert( graph->findNextNode(join[pop_index], join_back[pop_index], (join_back_counters[pop_index] != &join[pop_index]->left_count[0])) == the_pop_segment );
                        //assert( (join_back_counters[pop_index])[ join_back[pop_index] ] != 0 );
                        (join_back_counters[pop_index])[ join_back[pop_index] ] = 0;

                        graph->removeNode(the_pop_segment); // remove node -- TODO TBB / non-nodelist graph iterator safety?
#ifdef VELOUR_TBB
                        tls_seqnode_allocator->DeallocateNodeMemory(the_pop_segment);
#else
                        g__SEQNODE_ALLOCATOR->DeallocateNodeMemory(the_pop_segment);
#endif // VELOUR_TBB
                        poppable_segment[pop_index] = NULL; // nullify pointer so segment isn't referenced again for comparison
                        ++ p__BUBBLES_POPPED; // update statistic

#ifdef VERIFY
                        graph->verify_node(node);
                        graph->verify_node(join[pop_index]);
#endif

                        modified = true;
                    }
                }
            }
        }
    }

    return modified;
}
*/

//
// PURPOSE: pop isolated and overlapping bubbles, on one side of the node
//
//     _o_                                                  //
//    /   \                                                 //
// --o--o--o--   ===>   --o--o--o--                         //
//
//                                                          //
//                                                          //
//     _o_                  _o_                             //
//    /   \                /   \                            //
// --o--o--o--   ===>   --o--o--o--   ===>   --o--o--o--    //
//    \   /                                                 //
//     -o-                                                  //
//                                                          
//     _o_                     _o_                          //
//    /   \                   /   \                         //
// --o--o--o--o--   ===>   --o--o--o--   ===>  --o--o--o--  //
//       \   /                                              //
//        -o-                                               //
//
// method:
//   1. determine trunk (coverage majority) vector from root node
//   2. discover all non-trunk poppable paths from root node,
//        memoize join node for each
//   3. follow trunk path, looking for the memoized join nodes
//
//  XXX NOTE: this code will not detect poppable paths comprised of non-concatenated nodes
//
//  XXX NOTE: ?? not well designed for use on a partitioned graph
//
//  NOTE: NOT TBB PARALLEL SAFE
//
//  return TRUE if modified the graph
//
static bool check_complex_pop_direction(SeqGraph *graph, SeqNode *node, bool direction)
{
    // PRE-CONDITION: do not process node if its subgraph is not present
    if (direction == GO_RIGHT && node->right_side_colors != 0) return false;
    if (direction == GO_LEFT && node->left_side_colors != 0) return false;

    SeqNode *next_arr[4] = {NULL};
    bool next_sense_changed_arr[4];

    counter_t *node_counters = direction ? node->right_count : node->left_count;

    //
    // determine the trunk vector (index of node's majority neighbor)
    //
    unsigned node_majority_idx;
    double next_majority_kmer_coverage = 0.0;
    for (unsigned i=0; i < 4 ; ++i) {
        if (node_counters[i] != 0) {
            next_arr[i] = graph->findNextNode(node, i, direction, &next_sense_changed_arr[i]);
            SeqNode *next = next_arr[i];
            assert( next != NULL );

            double next_kmer_coverage = next->getNodeKmerCoverage(g__FULLKMER_LENGTH);
            if (next_kmer_coverage > next_majority_kmer_coverage) {
                node_majority_idx = i;
                next_majority_kmer_coverage = next_kmer_coverage;
            }
        }
    }
    assert(next_majority_kmer_coverage != 0.0);
            
    //
    // discover poppable segments
    //
    SeqNode *poppable_segment[4] = {NULL};
    unsigned poppable_count = 0;

    SeqNode *join[4] = {NULL};
    counter_t *join_back_counters[4] = {NULL}; // XXX OPT: store pointer to single counter, obviating join_back?
    Nucleotide join_back[4];

    for (unsigned i=0 ; i < 4 ; ++i) {
        if (i == node_majority_idx) continue; // skip the trunk
        if (node_counters[i] == 0) continue;

        bool next_sense_changed = next_sense_changed_arr[i];
        SeqNode *next = next_arr[i];
        assert( next != NULL );
        
        if (next->sequence.get_length() > MAX_BUBBLE_LENGTH) continue; // early exit: node is too long to be a bubble
        // TODO: if (next->getNodeKmerCoverage(g__FULLKMER_LENGTH) > MAX_BUBBLE_COVERAGE) continue; // early exit: node is not poppable

        bool next_direction = direction ^ next_sense_changed;

        // check the backward direction for the given segment
        counter_t *next_back_counters = next_direction ? next->left_count : next->right_count;
        int next_back_arc = valid_single_successor(next_back_counters);
        assert( next_back_arc != NO_SUCCESSORS );
        if (next_back_arc == MULTIPLE_SUCCESSORS) continue;

        // check the forward direction and memoize the join node
        counter_t *next_counters = next_direction ? next->right_count : next->left_count;
        color_t *next_colors = next_direction ? next->right_color : next->left_color;

        int next_arc = valid_single_successor(next_counters);
        if (next_arc != NO_SUCCESSORS && next_arc != MULTIPLE_SUCCESSORS) {

            if (next_colors[next_arc] != 0) continue;  // ok, bubble popping is unordered

            bool next_join_sense_changed;
            SeqNode *next_join = graph->findNextNode(next, next_arc, next_direction, &next_join_sense_changed);
            assert( next_join != NULL );

            bool next_join_direction = next_direction ^ next_join_sense_changed;

            // compute this once we know about next_join
            Nucleotide next_head_rightmost_base = next->sequence.GetHeadKmerRightmostBase(g__FULLKMER_LENGTH);
            Nucleotide next_tail_leftmost_base = next->sequence.GetTailKmerLeftmostBase(g__FULLKMER_LENGTH);
            Nucleotide next_join_back_valid = next_direction ? next_tail_leftmost_base : next_head_rightmost_base;
            if (next_join_sense_changed) { next_join_back_valid = COMPLEMENT(next_join_back_valid); }

            // record as poppable
            join[i] = next_join;
            join_back_counters[i] = next_join_direction ? next_join->left_count : next_join->right_count; // the back counters
            join_back[i] = next_join_back_valid;

            poppable_segment[i] = next;
            ++ poppable_count;
        }
    }

    if (poppable_count == 0) { // nothing to pop
        return false;
    }

    // at this point, we've identified the join nodes and linear paths that are short enough to be poppable

    bool modified = false;

    //
    // look for join nodes connected to the trunk path
    //
    //   XXX NOTE: significant lost opportunity, since doesn't check all paths
    //
    SeqNode *segment_root = next_arr[node_majority_idx];
    bool root_direction = direction ^ next_sense_changed_arr[node_majority_idx];

    // initialize consensus sequence
    Sequence_StackAllocated sequence_memory;
    Sequence *consensus = new (&sequence_memory) Sequence(segment_root->sequence, Sequence::MAX_BASES);
    bool consensus_sense_changed = false;

    SeqNode *from = segment_root;
    bool from_direction = root_direction;
    //while (1) {
        // look for join nodes, and if found, try to pop bubbles
        counter_t *from_counters = from_direction ? from->right_count : from->left_count;
        for (unsigned i=0; i < 4; ++i) { // 'i' is the idx from the latest trunk node 'from'
            if (from_counters[i] == 0) continue;
            bool peek_sense_changed;
            SeqNode *peek = graph->findNextNode(from, i, from_direction, &peek_sense_changed);
            assert( peek != NULL );
            //bool peek_direction = from_direction ^ to_sense_changed;

            for (unsigned k=0; k < 4; ++k) { // 'k' is the idx from the original branching node
                if (peek == join[k]) {
                    //
                    // found a bubble join, try to pop the bubble
                    //
                    SeqNode *the_pop_segment = poppable_segment[k];

                    ++ p__BUBBLES_CONSIDERED;

                    // check if the bubble can be popped
                    bool pop = compareSequences(&the_pop_segment->sequence, consensus);
                    if (!pop) {
                        // TODO OPT: don't RC unless the two sequences are potentially similar, e.g. length
                        Sequence_StackAllocated sequence_memory;
                        Sequence *revcomp = new (&sequence_memory) Sequence(the_pop_segment->sequence, Sequence::MAX_BASES);
                        revcomp->ReverseComplement();
                        pop = compareSequences(revcomp, consensus);
                    }

                    if (pop) {
                        // FIXME: this is only correct if trunk segment is a single node!!!
                        segment_root->kmer_occurrences += the_pop_segment->kmer_occurrences;
                        // END FIXME

                        assert( node_counters[k] != 0 );
                        node_counters[k] = 0;

                        assert( (join_back_counters[k])[ join_back[k] ] != 0 );
                        (join_back_counters[k])[ join_back[k] ] = 0;

                        setNodeDead<SeqNode>(the_pop_segment);
                        the_pop_segment->connections = 0;
/*
                        graph->removeNode(the_pop_segment); // remove node
#ifdef VELOUR_TBB
                        tls_seqnode_allocator->DeallocateNodeMemory(the_pop_segment);
#else
                        g__SEQNODE_ALLOCATOR->DeallocateNodeMemory(the_pop_segment);
#endif // VELOUR_TBB
*/

                        join[k] = NULL; // nullify pointer so segment isn't popped again
                        ++ p__BUBBLES_POPPED; // update statistic

/*
#ifdef VERIFY
                        graph->verify_node(node);
                        graph->verify_node(join[k]);
#endif
*/

                        -- poppable_count;
                        modified = true;
                    }
                }
                if (poppable_count == 0) break;
            }
            if (poppable_count == 0) break;
        }
        //if (poppable_count == 0) break;

        //consensus_sense_changed ^= to_sense_changed;

        //bool to_sense_changed;
        //SeqNode *to = graph->findNextNode(from, i, from_direction, &to_sense_changed);
        //assert( to != NULL );
        //bool to_direction = from_direction ^ to_sense_changed;

        // TODO TODO TODO: keep adding to the trunk (e.g. below), but have to spread out the coverage update
    //}

    return modified;

    // OLD OLD OLD
/*
    for (unsigned i=0 ; i < 4 ; ++i) {
        if (node_counters[i] == 0 || node_colors[i] != 0) continue;
        if (poppable_segment[i] != NULL) continue; // don't check direction that is poppable
        if (poppable_count == 0) break;

        bool next_sense_changed;
        SeqNode *next = graph->findNextNode(node, i, direction, &next_sense_changed);
        assert( next != NULL );

        bool next_direction = direction ^ next_sense_changed;

        if (next->sequence.get_length() > MAX_BUBBLE_LENGTH) continue; // early exit: too long to be a bubble consensus

        // don't check the backward direction

        // NOTE: this node can't be a join node, since no segment in between

        counter_t *next_counters = next_direction ? next->right_count : next->left_count;
        color_t *next_colors = next_direction ? next->right_color : next->left_color;

        // determine trunk forward direction (non-single copy) // TODO TODO TODO
        int trunk_arc = -1;
        for (unsigned k=0 ; k < 4 ; ++k) {   // find first trunk direction and take it
            if (next_counters[k] > MAX_BUBBLE_COVERAGE) {
                trunk_arc = k;
                break;
            }
        }
        if (trunk_arc == -1) continue;    // no trunk to follow

        if (next_colors[trunk_arc] != 0) continue;

        // initialize consensus sequence
        Sequence_StackAllocated sequence_memory;
        Sequence *consensus = new (&sequence_memory) Sequence(next->sequence, Sequence::MAX_BASES);

        bool consensus_sense_changed = false;

        //int depth = 1;
        //printf("consensus(%02d): ", depth); consensus->Print(); puts("\n");

        // okay, lets go looking for join nodes!
        SeqNode *from = next;
        bool from_direction = next_direction;
        while (trunk_arc >= 0) {

            counter_t *from_counters = from_direction ? from->right_count : from->left_count;
            color_t *from_colors = from_direction ? from->right_color : from->left_color;

            assert( from_counters[trunk_arc] != 0 );
            if (from_colors[trunk_arc] != 0) break;

            bool to_sense_changed;
            SeqNode *to = graph->findNextNode(from, trunk_arc, from_direction, &to_sense_changed);
            assert( to != NULL );

            bool to_direction = from_direction ^ to_sense_changed;

            consensus_sense_changed ^= to_sense_changed;

            // check if "to" is a join node
            for (unsigned k=0 ; k < 4 ; ++k) { // k = the node arc the poppable segment is sourced from
                if (to == join[k]) { // found a join node!
                    SeqNode *the_pop_segment = poppable_segment[k];

                    ++ p__BUBBLES_CONSIDERED;

                    // check if the bubble can be popped
                    bool pop = compareSequences(&the_pop_segment->sequence, consensus);
                    if (!pop) {
                        // TODO OPT: don't RC unless the two sequences are potentially similar, e.g. length
                        Sequence_StackAllocated sequence_memory;
                        Sequence *revcomp = new (&sequence_memory) Sequence(the_pop_segment->sequence, Sequence::MAX_BASES);
                        revcomp->ReverseComplement();
                        pop = compareSequences(revcomp, consensus);
                    }

                    if (pop) {

#ifdef VERIFY
                        graph->verify_node(node);
                        graph->verify_node(the_pop_segment);
                        graph->verify_node(join[k]);
#endif

                        //assert( cnorm(node_counters[pop_index]) > 0 && cnorm(node_counters[pop_index]) <= MAX_BUBBLE_COVERAGE );
                        //assert( graph->findNextNode(node, pop_index, direction) == the_pop_segment );
                        //assert( node_counters[pop_index] != 0 );
                        node_counters[k] = 0;

                        //assert( cnorm((join_back_counters[pop_index])[ join_back[pop_index] ]) > 0 &&
                        //       cnorm((join_back_counters[pop_index])[ join_back[pop_index] ]) <= MAX_BUBBLE_COVERAGE );
                        //assert( graph->findNextNode(join[pop_index], join_back[pop_index], (join_back_counters[pop_index] != &join[pop_index]->left_count[0])) == the_pop_segment );
                        //assert( (join_back_counters[pop_index])[ join_back[pop_index] ] != 0 );
                        (join_back_counters[k])[ join_back[k] ] = 0;

                        graph->removeNode(the_pop_segment); // remove node
#ifdef VELOUR_TBB
                        tls_seqnode_allocator->DeallocateNodeMemory(the_pop_segment);
#else
                        g__SEQNODE_ALLOCATOR->DeallocateNodeMemory(the_pop_segment);
#endif // VELOUR_TBB

                        join[k] = NULL; // nullify pointer so segment isn't popped again
                        ++ p__BUBBLES_POPPED; // update statistic

#ifdef VERIFY
                        graph->verify_node(node);
                        graph->verify_node(join[k]);
#endif

                        -- poppable_count;
                        modified = true;
                    }
                }
            }

            if (poppable_count == 0) break;

            // early exit if too long to be a bubble consensus
            if (consensus->get_length() + to->sequence.get_length() > MAX_BUBBLE_LENGTH) break;

            //printf("   to        : ", depth); to->sequence.Print(); puts("\n");

            consensus->Concatenate_Unsafe(&to->sequence, next_direction, consensus_sense_changed, g__FULLKMER_LENGTH);

            //++ depth;
            //printf("consensus(%02d): ", depth); consensus->Print(); puts("\n");

            from = to;
            from_direction = to_direction;

            // determine next trunk direction
            trunk_arc = -1;
            for (unsigned k=0 ; k < 4 ; ++k) {
                counter_t *counters = from_direction ? from->right_count : from->left_count;
                if (counters[k] > MAX_BUBBLE_COVERAGE) {
                    trunk_arc = k;
                    break;
                }
            }
        }
    }
    return modified;
*/
}

namespace {

struct serial_pop_bubbles_functor {
    serial_pop_bubbles_functor(SeqGraph *graph, bool *modified) : graph(graph), modified(modified) {}
    void operator()(SeqNode *node) {
        assert( node != NULL );
        if (isNodeDead(node)) return;

        //if (isNodeWeight_AllSingleCopy(node)) { ++ p__PROFILE_SINGLECOPY; }
        //if (isNodeWeight_None(node)) { ++ p__PROFILE_NONE; }

#ifdef SMALL_NODES
        unsigned left_count = 0, right_count = 0;
        node->getNeighborCounts(&left_count, &right_count);
#else
        // ?
#endif

        if (left_count >= 2) {
            *modified |= check_complex_pop_direction(graph, node, GO_LEFT);
        }

        if (right_count >= 2) {
            *modified |= check_complex_pop_direction(graph, node, GO_RIGHT);
        }
    }
  private:
    SeqGraph *graph;
    bool *modified;
};

} // namespace: anonymous

// serial
void sg_pop_bubbles(SeqGraph *graph, bool silent)
{
  if (!g__BUBBLE_POPPING) return;

	bool modified = true;
	unsigned pass = 1;

	while (modified) {
        uintptr_t start_node_count = graph->node_count;
        modified = false;

        if (!silent) {
            printf("==================== SG: BUBBLE POPPING PASS %d ====================\n", pass);
            fflush(stdout);
        }

        serial_pop_bubbles_functor f(graph, &modified);
        sg_for_each_mutable(graph, f);

        // TODO: do a cleaning pass to delete dead nodes to update graph->node_count
        uintptr_t nodes_deallocated_this_time = (start_node_count - graph->node_count);
        p__BUBBLE_NODES_DEALLOCATED += nodes_deallocated_this_time;

        if (!silent) {
            printf("sg_pop_bubbles: %"PRIuPTR" of %"PRIuPTR" bubbles popped, %"PRIuPTR" nodes deallocated\n",
                    p__BUBBLES_POPPED, p__BUBBLES_CONSIDERED, p__BUBBLE_NODES_DEALLOCATED);
            fflush(stdout);
        }

        ++pass;

        // XXX: concatenate if don't support multi-linear node traversal
        if (modified) {
            struct lambda {
                static void id(SeqNode *) {}
            };
            sg_for_each_mutable(graph, lambda::id); // make sure things are deleted, TODO OPT delete this

            sg_concatenate(graph, silent);
        }
    }
    //sg_concatenate(graph, silent);

    //printf("profiling bubble potential: %"PRIuPTR" none, %"PRIuPTR" single copy, %"PRIuPTR" total\n",
    //        p__PROFILE_NONE, p__PROFILE_SINGLECOPY, graph->node_count);
    //fflush(stdout);
}

/*
// serial
void sg_nodelist_pop_bubbles(SeqGraph *graph, bool silent, flow_nodelist_t *nodelist)
{
  if (!g__BUBBLE_POPPING) return;

	bool modified = true;
	unsigned pass = 1;

	while (modified) {
        uintptr_t start_node_count = graph->node_count;
        modified = false;

        if (!silent) {
            printf("==================== SG: BUBBLE POPPING PASS %d ====================\n", pass);
            fflush(stdout);
        }

        serial_pop_bubbles_functor f(graph, &modified);
        for (flow_nodelist_t::iterator it = nodelist->begin(); it != nodelist->end(); ++it) {
            SeqNode *node = *it;
            if (!isNodeDead(node)) { f(node); }
        }

        uintptr_t nodes_deallocated_this_time = (start_node_count - graph->node_count);
        p__BUBBLE_NODES_DEALLOCATED += nodes_deallocated_this_time;

        if (!silent) {
            printf("sg_pop_bubbles: %"PRIuPTR" of %"PRIuPTR" bubbles popped, %"PRIuPTR" nodes deallocated\n",
                    p__BUBBLES_POPPED, p__BUBBLES_CONSIDERED, p__BUBBLE_NODES_DEALLOCATED);
            fflush(stdout);
        }

        ++pass;

        // XXX: concatenate if don't support multi-linear node traversal
        if (modified) {
            sg_nodelist_concatenate(graph, silent, nodelist);
        }
    }
    //sg_nodelist_concatenate(graph, silent, nodelist);
    printf("profiling bubble potential: %"PRIuPTR" none, %"PRIuPTR" single copy, %"PRIuPTR" total\n",
            p__PROFILE_NONE, p__PROFILE_SINGLECOPY, graph->node_count);
    fflush(stdout);
}
*/
