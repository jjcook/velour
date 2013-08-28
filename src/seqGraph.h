//
// This file is part of Velour and is distributed under the University of
// Illinois Open Source License.  See LICENSE.txt for details.
//

//
// seqGraph.h
//
//   sequence graph: variable sized nodes
//

#ifndef __SEQUENCE_GRAPH_H
#define __SEQUENCE_GRAPH_H

#include "seqNode.h"

class SeqGraph {
  public:
    SeqGraph(uintptr_t buckets);
    ~SeqGraph();

    void insertNode(SeqNode *);
    void insertNodeAndUpdateColors(SeqNode *);

    void insertNodeHead(SeqNode *);
    void insertNodeTail(SeqNode *);

    SeqNode *findNextNode(const SeqNode *, Nucleotide next_base, bool moving_right, bool *sense_changed_ptr=NULL) const;

    void removeNodeUnsafe(SeqNode *node);
    void removeNode(SeqNode *);

    void removeNodeHead(SeqNode *);
    void removeNodeTail(SeqNode *);

#ifdef VELOUR_TBB
    void atomic_insertNode(SeqNode *);

    void atomic_insertNodeHead(SeqNode *);
    void atomic_insertNodeTail(SeqNode *);

    SeqNode *atomic_findNextNode(const SeqNode *, Nucleotide next_base, bool moving_right, bool *sense_changed_ptr=NULL) const;

    void atomic_removeNodeUnsafe(SeqNode *node);
    void atomic_removeNode(SeqNode *);
    void atomic_removeNodeHead(SeqNode *);
    void atomic_removeNodeTail(SeqNode *);
#endif // VELOUR_TBB

    template<typename T> friend void sg_for_each(SeqGraph*, T& unary_functor);
    template<typename T> friend void sg_for_each_grow(SeqGraph*, T& unary_functor);
    template<typename T> friend void sg_for_each_mutable(SeqGraph*, T& unary_functor);

    void verify(bool silent=false); // XXX const
    void verify_node(SeqNode *node); // XXX const

    //bool isMarked(void);
    //void setMarked(void);

    void resetFlags(void);
    void resetSliceFlagsAndSliceColor(void);

    void getNeighborCounts(SeqNode *node, unsigned *l_count, unsigned *r_count) const;

    void bulkMoveAllNodes(SeqGraph *target);

    void PrepareForMajorGarbageCollect(void);

  private:
    typedef uintptr_t hash_index_t;
    hash_index_t hash(const Kmer) const;
 
    SeqNode* findNode(const Kmer, bool moving_right, bool* sense_changed) const;
#ifdef VELOUR_TBB
    SeqNode* atomic_findNode(const Kmer, bool moving_right, bool* sense_changed) const;
#endif // VELOUR_TBB

  private:
    const uintptr_t buckets;
    const uintptr_t hashmask;
    SeqNode **head_table;
    SeqNode **tail_table;

    //bool is_marked; // non-empty flags

  public: // XXX: privatize 
    // pseudo_edge_count;
    uintptr_t node_count;

#ifdef VERIFY
#  ifdef VELOUR_TBB
    tbb::atomic<unsigned> iterator_lock;
#  else
    unsigned iterator_lock;
#  endif
#endif
};

inline void SeqGraph::insertNode(SeqNode *node)
{
#ifdef VERIFY
    // XXX: disabled: only happens in import related; the insertions aren't interesting to the iterator
    //assert( this->iterator_lock == 0 && "Can't insert/remove in a graph iterator.");
#endif // VERIFY
    insertNodeHead(node);
    insertNodeTail(node);
    ++node_count;
}

inline void SeqGraph::insertNodeHead(SeqNode *node)
{
	Kmer head_kmer = node->sequence.GetHeadKmer(g__FULLKMER_LENGTH);
	Kmer head_canon_kmer = canonicalKmer(head_kmer, g__FULLKMER_LENGTH);
    hash_index_t hash_index = hash(head_canon_kmer);
    SeqNode **head_bucket = &head_table[hash_index];

#ifdef VERIFY
	// check that another node in this bucket with this head kmer doesn't exist
    for (SeqNode *trav = *head_bucket; trav != NULL ; trav = trav->head_next) {
        assert( trav != node );
        Kmer trav_head_kmer = trav->sequence.GetHeadKmer(g__FULLKMER_LENGTH);
        assert( trav_head_kmer != head_kmer );
    }
#endif

    node->head_next = *head_bucket;
    *head_bucket = node;
}

inline void SeqGraph::insertNodeTail(SeqNode *node)
{
	Kmer tail_kmer = node->sequence.GetTailKmer(g__FULLKMER_LENGTH);
	Kmer tail_canon_kmer = canonicalKmer(tail_kmer, g__FULLKMER_LENGTH);
    hash_index_t hash_index = hash(tail_canon_kmer);
    SeqNode **tail_bucket = &tail_table[hash_index];

#ifdef VERIFY
	// check that another node in this bucket with this tail kmer doesn't exist
    for (SeqNode *trav = *tail_bucket; trav != NULL ; trav = trav->tail_next) {
        assert( trav != node );
        Kmer trav_tail_kmer = trav->sequence.GetTailKmer(g__FULLKMER_LENGTH);
        assert( trav_tail_kmer != tail_kmer );
    }
#endif

    node->tail_next = *tail_bucket;
    *tail_bucket = node;
}

#ifdef VELOUR_TBB
inline void SeqGraph::atomic_insertNode(SeqNode *node)
{
#ifdef VERIFY
    // XXX: disabled: only happens in import related; the insertions aren't interesting to the iterator
    //assert( this->iterator_lock == 0 && "Can't insert/remove in a graph iterator.");
#endif // VERIFY
    atomic_insertNodeHead(node);
    atomic_insertNodeTail(node);
    ATOMIZE(node_count).fetch_and_increment();
}

inline void SeqGraph::atomic_insertNodeHead(SeqNode *node)
{
	Kmer head_kmer = node->sequence.GetHeadKmer(g__FULLKMER_LENGTH);
	Kmer head_canon_kmer = canonicalKmer(head_kmer, g__FULLKMER_LENGTH);
    hash_index_t hash_index = hash(head_canon_kmer);
    SeqNode **head_bucket = &head_table[hash_index];

/*#ifdef VERIFY
	// check that another node in this bucket with this head kmer doesn't exist
    for (SeqNode *trav = *head_bucket; trav != NULL ; trav = trav->head_next) {
        assert( trav != node );
        Kmer trav_head_kmer = trav->sequence.GetHeadKmer(g__FULLKMER_LENGTH);
        assert( trav_head_kmer != head_kmer );
    }
#endif*/

    while(1) {
        SeqNode *first_node;
        while ((first_node = ATOMIC_LOAD(*head_bucket)) == reinterpret_cast<SeqNode*>(1)) {} // spin waiting for unlock
        node->head_next = first_node;
        if (ATOMIZE(*head_bucket).compare_and_swap(node, first_node) == first_node) {
            break;
        }
    }
}

inline void SeqGraph::atomic_insertNodeTail(SeqNode *node)
{
	Kmer tail_kmer = node->sequence.GetTailKmer(g__FULLKMER_LENGTH);
	Kmer tail_canon_kmer = canonicalKmer(tail_kmer, g__FULLKMER_LENGTH);
    hash_index_t hash_index = hash(tail_canon_kmer);
    SeqNode **tail_bucket = &tail_table[hash_index];

/*#ifdef VERIFY
	// check that another node in this bucket with this tail kmer doesn't exist
    for (SeqNode *trav = *tail_bucket; trav != NULL ; trav = trav->tail_next) {
        assert( trav != node );
        Kmer trav_tail_kmer = trav->sequence.GetHeadKmer(g__FULLKMER_LENGTH);
        assert( trav_tail_kmer != tail_kmer );
    }
#endif*/

    while(1) {
        SeqNode *first_node;
        while ((first_node = ATOMIC_LOAD(*tail_bucket)) == reinterpret_cast<SeqNode*>(1)) {} // spin waiting for unlock
        node->tail_next = first_node;
        if (ATOMIZE(*tail_bucket).compare_and_swap(node, first_node) == first_node) {
            break;
        }
    }
}
#endif // VELOUR_TBB

template<typename T> inline void sg_for_each(SeqGraph* graph, T& unary_functor)
{
    if (graph->node_count == 0) { return; }
#ifdef VERIFY
    ++ graph->iterator_lock;
#endif // VERIFY
    for(uintptr_t i = 0 ; i < graph->buckets ; ++ i) {
        for (SeqNode *trav = graph->head_table[i] ; trav != NULL ; trav = trav->head_next) {
            if (isNodeDead<SeqNode>(trav)) continue;
            unary_functor(trav);
        }
    }
#ifdef VERIFY
    assert( graph->iterator_lock > 0 );
    -- graph->iterator_lock;
#endif // VERIFY
}

template<typename T> inline void sg_for_each_grow(SeqGraph* graph, T& unary_functor)
{
    if (graph->node_count == 0) { return; }
#ifdef VERIFY
    ++ graph->iterator_lock;
#endif // VERIFY
    for(uintptr_t i = 0 ; i < graph->buckets ; ++ i) {
        SeqNode **bucket = &graph->head_table[i];
        SeqNode **prev = bucket;
        for (SeqNode *trav = *prev ; trav != NULL ; trav = *prev) {
            if (isNodeDead<SeqNode>(trav)) continue;
            SeqNode *retnode = unary_functor(trav);
            if (retnode != trav) { // node was reallocated, re-check this bucket // TODO OPT somehow don't reprocess previos nodes?
                prev = bucket;
                if (*prev == retnode) { // skip the node we just processed
                    prev = &( (*prev)->head_next );
                }
            } else {
                prev = &(trav->head_next);
            }
        }
    }
#ifdef VERIFY
    assert( graph->iterator_lock > 0 );
    -- graph->iterator_lock;
#endif // VERIFY
}

template<typename T> inline void sg_for_each_mutable(SeqGraph* graph, T& unary_functor)
{
#ifdef VERIFY
    assert( graph->iterator_lock == 0 && "Can't mutate under another graph iterator.");
#endif // VERIFY
    if (graph->node_count == 0) { return; }
#ifdef VERIFY
    ++ graph->iterator_lock;
#endif // VERIFY
    for(uintptr_t i = 0 ; i < graph->buckets ; ++ i) {
        SeqNode **prev = &graph->head_table[i];
        for (SeqNode *trav = *prev ; trav != NULL ; trav = *prev) {
            bool remove_now = isNodeDead<SeqNode>(trav);
            if (remove_now) {
                *prev = trav->head_next;   // remove node from 'head' hashtable
                graph->removeNodeTail(trav); // remove node from 'tail' hashtable
#ifdef VELOUR_TBB
                tls_seqnode_allocator->DeallocateNodeMemory(trav);
#else
                g__SEQNODE_ALLOCATOR->DeallocateNodeMemory(trav);
#endif // VELOUR_TBB
                --(graph->node_count);
                continue;
            }
            if (!isNodeDead<SeqNode>(trav)) {
                unary_functor(trav);
            }
            remove_now = isNodeDead<SeqNode>(trav);
            if (remove_now) {
                *prev = trav->head_next;   // remove node from 'head' hashtable
                graph->removeNodeTail(trav); // remove node from 'tail' hashtable
#ifdef VELOUR_TBB
                tls_seqnode_allocator->DeallocateNodeMemory(trav);
#else
                g__SEQNODE_ALLOCATOR->DeallocateNodeMemory(trav);
#endif // VELOUR_TBB
                --(graph->node_count);
                continue;
            }
            prev = &(trav->head_next);
        }
    }
#ifdef VERIFY
    assert( graph->iterator_lock > 0 );
    -- graph->iterator_lock;
#endif // VERIFY
}

inline void SeqGraph::removeNodeUnsafe(SeqNode *node)
{
    removeNodeHead(node);
    removeNodeTail(node);
    assert( node_count > 0 );
    --node_count;
}

inline void SeqGraph::removeNode(SeqNode *node)
{
#ifdef VERIFY
    assert( this->iterator_lock == 0 && "Can't remove in a graph iterator.");
#endif // VERIFY
    removeNodeUnsafe(node);
}

#ifdef VELOUR_TBB
inline void SeqGraph::atomic_removeNodeUnsafe(SeqNode *node)
{
    atomic_removeNodeHead(node);
    atomic_removeNodeTail(node);
    assert( node_count > 0 );
    ATOMIZE(node_count).fetch_and_decrement();
}

inline void SeqGraph::atomic_removeNode(SeqNode *node)
{
#ifdef VERIFY
    assert( this->iterator_lock == 0 && "Can't remove in a graph iterator.");
#endif // VERIFY
    atomic_removeNodeUnsafe(node);
}
#endif // VELOUR_TBB

inline SeqGraph::hash_index_t SeqGraph::hash(const Kmer kmer) const {
    return ((kmer ^ (kmer >> 17) ^ (kmer >> 37) ^ (kmer << 11)) & hashmask);
}

/*    
inline bool SeqGraph::isMarked(void)
{
    return this->is_marked;
}

inline void SeqGraph::setMarked(void)
{
    this->is_marked = true;
}
*/

inline void SeqGraph::resetFlags(void)
{
    struct lambda {
        static void f(SeqNode *node) { resetNodeMergingFlags<SeqNode>(node); }
    };
    sg_for_each(this, lambda::f);
}

inline void SeqGraph::resetSliceFlagsAndSliceColor(void)
{
    struct lambda {
        static void f(SeqNode *node) { resetNodeSliceFlags(node); node->slice_color = 0; }
    };
    sg_for_each(this, lambda::f);
}

inline void SeqGraph::PrepareForMajorGarbageCollect(void)
{
    bzero(head_table, buckets * sizeof(SeqNode *));
    bzero(tail_table, buckets * sizeof(SeqNode *));
    node_count = 0;
}

#ifdef VELOUR_TBB
// hack: for parallel component
static inline SeqNode * findNextNode_helper(SeqGraph *graph, SeqNode *node, int i, int direction)
{
    return graph->atomic_findNextNode(node, i, direction);
}
#endif // VELOUR_TBB

//
// construction, loading, and emittal
//

// FIXME: move these?
    
void sg_build(SeqGraph *sg_hashtable, KmerGraph *kg_hashtable);
void sg_build_component(KmerGraph *kg_hashtable, std::deque<KmerNode*> &nodes, SeqGraph *sg_hashtable);

void sg_concatenate(SeqGraph *sg_hashtable, bool silent=false);
void sg_nodelist_concatenate(SeqGraph *sg_hashtable, bool silent, flow_nodelist_t *);

void sg_load_bucket(SeqGraph *hashtable, const char *filename);
uintptr_t sg_load_stream_bucket(SeqGraph *hashtable, const uintptr_t node_limit, FILE *file);
uintptr_t sg_load_stream_component(SeqGraph *sgraph, FILE *file);
uintptr_t sg_load_mmap_stream_component(char *mmap_start, size_t *amount, SeqNode **chain);
void sg_load_quilt(SeqGraph *hashtable, const char *filename);
void sg_dump_quilt(SeqGraph *hashtable, const char *filename);
void sg_dump_pregraph(SeqGraph *hashtable, const char *filename);

void sg_dump_quilt_from_kmergraph(KmerGraph *kg_hashtable, const char *filename);
void sg_dump_pregraph_from_kmergraph(KmerGraph *kg_hashtable, const char *filename);

struct PrivateSeqGraphSet
{
    static const unsigned MIN_PRIVATE_GRAPH_LOG2 = 8;
    static const unsigned MAX_PRIVATE_GRAPH_LOG2 = 22;

    PrivateSeqGraphSet()
    {
        for (unsigned log2_index = 0; log2_index < MIN_PRIVATE_GRAPH_LOG2; ++ log2_index) {
            psg_list[log2_index] = NULL;
        }
        for (unsigned log2_index = MIN_PRIVATE_GRAPH_LOG2; log2_index <= MAX_PRIVATE_GRAPH_LOG2; ++ log2_index) {
            psg_list[log2_index] = new SeqGraph(1 << log2_index);
        }
    }

    ~PrivateSeqGraphSet()
    {
        for (unsigned log2_index = MIN_PRIVATE_GRAPH_LOG2; log2_index <= MAX_PRIVATE_GRAPH_LOG2; ++ log2_index) {
            delete psg_list[log2_index];
        }
    }

    SeqGraph * GetPrivateGraph(uintptr_t nodes)
    {
        SeqGraph *private_sgraph = NULL;
        uintptr_t desired_buckets = round_up_power_of_2(nodes >> 2); // avg 4 nodes per bucket
        unsigned log2_index = log2_of_power_of_2(desired_buckets);
        if (log2_index > MAX_PRIVATE_GRAPH_LOG2) {
            printf("(performance) unexpected private sgraph of %"PRIuPTR" buckets wanted.  using %"PRIuPTR" instead.\n",
                    desired_buckets, (static_cast<uintptr_t>(1) << MAX_PRIVATE_GRAPH_LOG2));
            fflush(stdout);
            private_sgraph = psg_list[MAX_PRIVATE_GRAPH_LOG2];
        } else if (log2_index < MIN_PRIVATE_GRAPH_LOG2) {
            private_sgraph = psg_list[MIN_PRIVATE_GRAPH_LOG2];
        } else {
            private_sgraph = psg_list[log2_index];
        }
        return private_sgraph;
    }

    private:
        SeqGraph *psg_list[MAX_PRIVATE_GRAPH_LOG2+1];
};

#endif // __SEQUENCE_GRAPH_H
