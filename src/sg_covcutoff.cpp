//
// This file is part of Velour and is distributed under the University of
// Illinois Open Source License.  See LICENSE.txt for details.
//

//
// sg_covcutoff.cpp
//
//   sequence graph coverage cutoff -- delete nodes with insufficient kmer coverage
//
// XXX NOTE: do not apply to partitioned graphs!!!

#include "types.h"

//
// CONFIGURABLE PARAMETERS
//
// NONE?

// internal statistics variables // TODO: need to reset some of these between iterations
static uintptr_t p__MINCOV_DELETED = 0;
static uintptr_t p__MAXCOV_DELETED = 0;

namespace {

struct serial_covcutoff_functor {
    serial_covcutoff_functor(SeqGraph *graph) : graph(graph) {}
    void operator()(SeqNode *node) {
        assert( node != NULL );
        if (isNodeDead(node)) return;

        double node_kmer_coverage = node->getNodeKmerCoverage(g__FULLKMER_LENGTH);

        bool remove_min = g__COVCUTOFF_MIN > 0.0 ? (node_kmer_coverage < g__COVCUTOFF_MIN) : false;
        bool remove_max = g__COVCUTOFF_MAX > 0.0 ? (node_kmer_coverage > g__COVCUTOFF_MAX) : false;

        if (remove_min || remove_max) {
            if (remove_min) ++ p__MINCOV_DELETED;
            if (remove_max) ++ p__MAXCOV_DELETED;

            // disconnect node from other nodes
            Nucleotide head_rightmost_base = node->sequence.GetHeadKmerRightmostBase(g__FULLKMER_LENGTH);
            Nucleotide tail_leftmost_base = node->sequence.GetTailKmerLeftmostBase(g__FULLKMER_LENGTH);
            for (int i=0; i < 4; ++i) {
                if (node->left_count[i] != 0) {
                    bool sense_changed;
                    SeqNode *next = graph->findNextNode(node, i, GO_LEFT, &sense_changed);
                    assert(next != NULL);

                    int next_back = sense_changed ? COMPLEMENT(head_rightmost_base) : head_rightmost_base;
                    counter_t *next_count = sense_changed ? next->left_count : next->right_count;
                    assert( next_count[next_back] != 0 );
                    next_count[next_back] = 0;
                }
                if (node->right_count[i] != 0) {
                    bool sense_changed;
                    SeqNode *next = graph->findNextNode(node, i, GO_RIGHT, &sense_changed);
                    assert(next != NULL);

                    int next_back = sense_changed ? COMPLEMENT(tail_leftmost_base) : tail_leftmost_base;
                    counter_t *next_count = sense_changed ? next->right_count : next->left_count;
                    assert( next_count[next_back] != 0 );
                    next_count[next_back] = 0;
                }
            }

            // delete node
            setNodeDead<SeqNode>(node);
            node->connections = 0;
/*            graph->removeNode(node);
#ifdef VELOUR_TBB
            tls_seqnode_allocator->DeallocateNodeMemory(node);
#else
            g__SEQNODE_ALLOCATOR->DeallocateNodeMemory(node);
#endif // VELOUR_TBB
*/
        }
    }
  private:
    SeqGraph *graph;
};

} // namespace: anonymous

// serial
void sg_covcutoff(SeqGraph *graph, bool silent)
{
    if (g__COVCUTOFF_MIN == 0.0 && g__COVCUTOFF_MAX == 0.0) return;

    if (!silent) {
        printf("==================== SG: COVERAGE CUTOFF: MIN %g MAX %g ====================\n",
                g__COVCUTOFF_MIN, g__COVCUTOFF_MAX);
        fflush(stdout);
    }

    serial_covcutoff_functor f(graph);
    sg_for_each_mutable(graph, f);

    if (!silent) {
        printf("sg_covcutoff: deallocated %"PRIuPTR" nodes below, %"PRIuPTR" nodes above\n",
                p__MINCOV_DELETED, p__MAXCOV_DELETED);
        fflush(stdout);
    }
}
