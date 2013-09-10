//
// This file is part of Velour and is distributed under the University of
// Illinois Open Source License.  See LICENSE.txt for details.
//

//
// kmerGraph.h
//
//   kmer graph: hashtable of kmer nodes
//
    
// TODO: statistics?

#ifndef __KMER_GRAPH_H
#define __KMER_GRAPH_H

class KmerGraph {
  public:
    KmerGraph(uintptr_t buckets);
    ~KmerGraph();

    KmerNode* findNodeOrAllocate(const Kmer);
    KmerNode* findNode(const Kmer) const;

    KmerNode *findNextNode(const KmerNode *, Nucleotide base, bool moving_right) const;

    void insertNode(KmerNode *node); // for garbage collection to rebuild hashtable
    //void removeNode(const Kmer);

    template<typename T> friend void kg_for_each(KmerGraph*, T& unary_functor, bool enable_garbage_collect=false);
    template<typename T> friend void kg_for_each_mutable(KmerGraph*, T& unary_functor);

#ifdef VELOUR_TBB
    template<typename T> friend void kg_for_each_parallel(KmerGraph *, T& unary_functor); // TODO: garbage collect for kmer->seq
#endif

    void resetFlags(void);

    void PrepareForMajorGarbageCollect(void);

    typedef uintptr_t hash_index_t;
  private:
    hash_index_t hash(const Kmer) const;

  public:
  //private:
    const uintptr_t buckets;
    const uintptr_t hashmask;
    KmerNode **table;
  
  public: // XXX: privatize
    // pseudo_edge_count;
    uintptr_t node_count;
};

template<typename T> inline void kg_for_each(KmerGraph* graph, T& unary_functor, bool enable_garbage_collect)
{
    if (graph->node_count == 0) return;
    for(uintptr_t i = 0 ; i < graph->buckets ; ++ i) {
        if (enable_garbage_collect && g__NODE_ALLOCATORS->flag_gc_needed_) { g__NODE_ALLOCATORS->GarbageCollect(); }
        for (KmerNode *trav = graph->table[i] ; trav != NULL ; trav = trav->next) {
            if (isNodeDead<KmerNode>(trav)) continue;
            unary_functor(trav);
        }
    }
}

template<typename T> inline void kg_for_each_mutable(KmerGraph* graph, T& unary_functor)
{
    if (graph->node_count == 0) return;
    for(uintptr_t i = 0 ; i < graph->buckets ; ++ i) {
        KmerNode **prev = &graph->table[i];
        for (KmerNode *trav = *prev ; trav != NULL ; trav = *prev) {
            bool remove_now = (trav->connections == 0);
            if (remove_now) {
                *prev = trav->next;
#ifdef VELOUR_TBB
                tls_kmernode_allocator->DeallocateNodeMemory(trav);
#else
                g__KMERNODE_ALLOCATOR->DeallocateNodeMemory(trav);
#endif // VELOUR_TBB
                --(graph->node_count);
                continue;
            }
            if (!isNodeDead<KmerNode>(trav)) {
                unary_functor(trav);
            }
            remove_now = (trav->connections == 0);
            if (remove_now) {
                *prev = trav->next;
#ifdef VELOUR_TBB
                tls_kmernode_allocator->DeallocateNodeMemory(trav);
#else
                g__KMERNODE_ALLOCATOR->DeallocateNodeMemory(trav);
#endif // VELOUR_TBB
                --(graph->node_count);
                continue;
            }
            prev = &(trav->next);
        }
    }
}

#ifdef VELOUR_TBB
namespace {

template<typename T>
struct ForEachInBucket
{
    ForEachInBucket(KmerGraph *kg, T& uf) : kgraph(kg), unary_functor(uf) {}

    void operator()(const tbb::blocked_range<KmerGraph::hash_index_t>& range) const {
        for (KmerGraph::hash_index_t i = range.begin(); i < range.end(); ++i) {
            for (KmerNode *trav = kgraph->table[i]; trav != NULL; trav = trav->next) {
                if (isNodeDead<KmerNode>(trav)) continue; // FIXME?
                unary_functor(trav);
            }
        }
    }

    KmerGraph *kgraph;
    T& unary_functor;
};

} // namespace: anonymous

template<typename T> void kg_for_each_parallel(KmerGraph* graph, T& unary_functor)
{
    if (graph->node_count == 0) return;
    tbb::parallel_for(tbb::blocked_range<KmerGraph::hash_index_t>(0,graph->buckets,1),
            ForEachInBucket<T>(graph, unary_functor)); // XXX OPT: grain size
}
#endif

inline KmerNode* KmerGraph::findNodeOrAllocate(const Kmer kmer)
{
    hash_index_t hash_index = hash(kmer);
    KmerNode **bucket_ptr = &table[hash_index];

#ifdef VELOUR_TBB
    KmerNode *new_node = NULL;
    KmerNode *first_node = NULL;

    for(;;) {
        first_node = ATOMIC_LOAD(*bucket_ptr);
        for (KmerNode *trav = *bucket_ptr ; trav != NULL ; trav = trav->next) {
            if (trav->kmer == kmer) {
                if (new_node != NULL) {
                    tls_kmernode_allocator->DeallocateNodeMemory(new_node);
                }
                ATOMIZE(trav->kmer_occurrences).fetch_and_increment();
                return trav;
            }
        }
        if (new_node == NULL) {
            KmerNode *mem = tls_kmernode_allocator->AllocateNodeMemory();
            new_node = new (mem) KmerNode(kmer);
        }
        new_node->next = first_node;
        KmerNode *old_node = ATOMIZE(*bucket_ptr).compare_and_swap(new_node, first_node);
        if (old_node == first_node) { // was not modified
            ATOMIZE(node_count).fetch_and_increment();
            return new_node;
        }
        // otherwise, loop again
    }
#else
    for (KmerNode *trav = *bucket_ptr ; trav != NULL ; trav = trav->next) {
        if (trav->kmer == kmer) {
            ++ trav->kmer_occurrences;
            return trav;
        }
    }
    KmerNode *mem = g__KMERNODE_ALLOCATOR->AllocateNodeMemory();
    KmerNode *node = new (mem) KmerNode(kmer);
    node->next = *bucket_ptr;
    *bucket_ptr = node;
    ++node_count;
    return node;
#endif
}

inline KmerNode* KmerGraph::findNode(const Kmer kmer) const {
    hash_index_t hash_index = hash(kmer);
/*#ifdef VELOUR_TBB
    KmerNode **bucket = &table[hash_index];
    KmerNode *retval;
    while (1) { // TODO: if NULL, don't do all this atomic jazz
        KmerNode *old_ptr;
        while ((old_ptr = axtomics::fetch_and_nop<KmerNode*>(bucket)) == reinterpret_cast<KmerNode*>(1)) {} // spin waiting for unlock
        if (__sxync_bool_compare_and_swap(bucket, old_ptr, reinterpret_cast<KmerNode*>(1))) { // lock bucket with negative one value
            for (KmerNode *trav = old_ptr ; trav != NULL ; trav = trav->next) {
                if (trav->kmer == kmer) {
                    retval = trav;
                    break;
                }
            }
            retval = NULL;

            // last, restore old pointer to unlock bucket
            bool unlocked = __sxync_bool_compare_and_swap(bucket, reinterpret_cast<KmerNode*>(1), old_ptr);
            assert( unlocked );
            break;
        }
    }
    return retval;
#else*/
    for (KmerNode *trav = table[hash_index] ; trav != NULL ; trav = trav->next) {
        if (trav->kmer == kmer) {
            return trav;
        }
    }
    return NULL;
//#endif // not VELOUR_TBB
}

inline KmerNode* KmerGraph::findNextNode(const KmerNode *node, Nucleotide base, bool moving_right) const
{
    unsigned double_kmer_length = g__FULLKMER_LENGTH << 1;
    Kmer next_kmer;
#ifdef LARGE_KMERS
    Kmer mask;
    mask.createMask(double_kmer_length);
#else
    Kmer mask = (Kmer(1) << double_kmer_length) - 1;
#endif
    if (moving_right) {
        next_kmer = KMER_APPEND(node->kmer, base, double_kmer_length);
    } else {
        next_kmer = KMER_PREPEND(node->kmer, base, double_kmer_length, mask);
    }
    Kmer canonical_kmer = canonicalKmer(next_kmer, g__FULLKMER_LENGTH);
    KmerNode *next = this->findNode(canonical_kmer);
    assert( next != NULL || (g__PSEUDO_NODES_PRESENT && ((moving_right ? node->right_color : node->left_color)[base] != 0)) );
    return next;
}

inline void KmerGraph::insertNode(KmerNode *node)
{
    hash_index_t hash_index = hash(node->kmer);
    KmerNode **bucket_ptr = &table[hash_index];

    node->next = *bucket_ptr;
    *bucket_ptr = node;
    ++node_count;
}

/* inline void KmerGraph::removeNode(const Kmer kmer)
{
    hash_index_t hash_index = hash(kmer);
#ifdef VELOUR_TBB
    bool deleted = false;
    KmerNode **bucket = &table[hash_index];
    while (1) {
        KmerNode *old_ptr;
        while ((old_ptr = axtomics::fetch_and_nop<KmerNode*>(bucket)) == reinterpret_cast<KmerNode*>(1)) {} // spin waiting for unlock
        if (__sxync_bool_compare_and_swap(bucket, old_ptr, reinterpret_cast<KmerNode*>(1))) { // lock bucket with negative one value
            KmerNode **prev = bucket;
            for (KmerNode *trav = old_ptr ; trav != NULL ; trav = *prev) {
                if (trav->kmer == kmer) {
                    *prev = trav->next;
                    tls_kmernode_allocator->DeallocateNodeMemory(trav);
                    --node_count;
                    deleted = true;
                    break;
                }
                prev = &(trav->next);
            }

            // last, restore old pointer to unlock bucket
            bool unlocked = __sxync_bool_compare_and_swap(bucket, reinterpret_cast<KmerNode*>(1), old_ptr);
            assert( unlocked );
            break;
        }
    }
    if (deleted) return; // to not print the error message
#else
    KmerNode *bucket_ptr = table[hash_index];
    KmerNode **prev = &bucket_ptr;
    for (KmerNode *trav = *prev ; trav != NULL ; trav = *prev) {
        if (trav->kmer == kmer) {
            *prev = trav->next;
            g__KMERNODE_ALLOCATOR->DeallocateNodeMemory(trav);
            --node_count;
            return;
        }
        prev = &(trav->next);
    }
#endif
    fprintf(stderr, "INTERNAL ERROR: KmerGraph::removeNode() called for non-existant kmer!\n");
    exit(EXIT_FAILURE);
} */

inline KmerGraph::hash_index_t KmerGraph::hash(const Kmer kmer) const {
    return ((kmer ^ (kmer >> 17) ^ (kmer >> 37) ^ (kmer << 11)) & hashmask);
}

inline void KmerGraph::resetFlags(void)
{
    struct lambda {
        static void f(KmerNode *node) { resetNodeMergingFlags<KmerNode>(node); }
    };
    kg_for_each(this, lambda::f);
}

inline void KmerGraph::PrepareForMajorGarbageCollect(void)
{
    bzero(table, buckets * sizeof(KmerNode *));
    node_count = 0;
}

#ifdef VELOUR_TBB
// hack: for parallel component
static inline KmerNode * findNextNode_helper(KmerGraph *graph, KmerNode *node, int i, int direction)
{
    return graph->findNextNode(node, i, direction);
}
#endif // VELOUR_TBB

#endif // __KMER_GRAPH_H
