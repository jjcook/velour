//
// This file is part of Velour and is distributed under the University of
// Illinois Open Source License.  See LICENSE.txt for details.
//

//
// node_allocators.h
//
//   managed slab allocators for kmer and sequence nodes
//

// XXX: rename DeallocateNodeMemory() since misleading?

// TODO: TBB atomize the update to live memory values? to ensure flushed to memory...

#ifndef VELOUR_NODE_ALLOCATORS_H
#define VELOUR_NODE_ALLOCATORS_H

class NodeAllocatorsBase
{
    public:
        static const size_t     SLAB_SIZE = (16 * 1024 * 1024);
        static const size_t TLS_SLAB_SIZE = ( 1 * 1024 * 1024);
        typedef std::deque<char *> SlabList;
};

class NodeAllocators;
class TLS_KmerNodeAllocator;
class TLS_SeqNodeAllocator;

class KmerNodeAllocator : public NodeAllocatorsBase
{
    friend class NodeAllocators;
    friend class TLS_KmerNodeAllocator;
    public:
        KmerNodeAllocator(NodeAllocators& parent) :
            node_allocators_(parent),
#ifdef VELOUR_TBB
            thread_allocators_(THREADID_MAX+2),
#endif // VELOUR_TBB
            slab_current_(NULL), slab_end_(NULL),
            kmer_graph_(NULL), kmer_node_memory_(0) { kmer_node_dead_memory_ = 0; /* for atomic version */ }

        ~KmerNodeAllocator() { /*for(SlabList::iterator it=used_list_.begin(); it != used_list_.end(); ++it) { free(*it); }*/ }

#ifndef VELOUR_TBB
        KmerNode * AllocateNodeMemory(void);
        void DeallocateNodeMemory(KmerNode *);
#else
        char * AllocateTlsSlab(void);
#endif // VELOUR_TBB

#ifdef VELOUR_TBB
        TLS_KmerNodeAllocator * GetTlsAllocator(threadid_t tid);
#endif // VELOUR_TBB

        void PrepareForGarbageCollect(void); // explicitly mark unallocated nodes as dead, then drop pointers to slab
        void MajorGarbageCollect(void);      // returns GC'd slabs to the parent free list
        //TODO void MinorGarbageCollect(void);
        void BulkFreeListAllSlabs(void);     // immediately move all slabs to the parent free list

        size_t GetKmerLiveMemory(void) { return (kmer_node_memory_ - kmer_node_dead_memory_ - (slab_end_ - slab_current_)); }

        void set_kmergraph(KmerGraph *kg) { kmer_graph_ = kg; }

    //private:
        NodeAllocators& node_allocators_;

#ifdef VELOUR_TBB
        tbb::queuing_mutex  kmer_node_allocator_lock_;
        tbb::concurrent_vector<TLS_KmerNodeAllocator *> thread_allocators_;
#endif // VELOUR_TBB

        char *slab_current_;
        char *slab_end_;
        
        KmerGraph * kmer_graph_;
        size_t kmer_node_memory_;

#ifdef VELOUR_TBB
        tbb::atomic<size_t> kmer_node_dead_memory_;
#else
        size_t kmer_node_dead_memory_;
#endif // VELOUR_TBB

        SlabList used_list_;    // NOTE: contains current slab
};

class SeqNodeAllocator : public NodeAllocatorsBase
{
    friend class NodeAllocators;
    friend class TLS_SeqNodeAllocator;
    public:
        SeqNodeAllocator(NodeAllocators& parent) :
            node_allocators_(parent),
#ifdef VELOUR_TBB
            thread_allocators_(THREADID_MAX+2),
#endif // VELOUR_TBB
            slab_current_(NULL), slab_end_(NULL),
            seq_graph_(NULL), seq_node_memory_(0) { seq_node_dead_memory_ = 0; /* for atomic version */ }

        ~SeqNodeAllocator() { /*for(SlabList::iterator it=used_list_.begin(); it != used_list_.end(); ++it) { free(*it); }*/ }

#ifndef VELOUR_TBB
        SeqNode * AllocateNodeMemory(size_t amount);
        SeqNode * ReallocateNodeMemory(SeqNode * old_node, size_t new_size);
        void DeallocateNodeMemory(SeqNode *);
#else
        char * AllocateTlsSlab(void);
#endif // VELOUR_TBB

#ifdef VELOUR_TBB
        TLS_SeqNodeAllocator * GetTlsAllocator(threadid_t tid);
#endif // VELOUR_TBB

        void PrepareForGarbageCollect(void); // explicitly mark unallocated nodes as dead, then drop pointers to slab
        void MajorGarbageCollect(void);      // returns GC'd slabs to the parent free list
        //TODO void MinorGarbageCollect(void);
        void BulkFreeListAllSlabs(void);     // immediately move all slabs to the parent free list

        size_t GetSeqLiveMemory(void) { return (seq_node_memory_ - seq_node_dead_memory_ - (slab_end_ - slab_current_)); }

        void set_seqgraph(SeqGraph *sg) { seq_graph_ = sg; }

    //private:
        NodeAllocators& node_allocators_;

#ifdef VELOUR_TBB
        tbb::queuing_mutex  seq_node_allocator_lock_;
        tbb::concurrent_vector<TLS_SeqNodeAllocator *> thread_allocators_;
#endif // VELOUR_TBB

        char *slab_current_;
        char *slab_end_;

        SeqGraph * seq_graph_;
        size_t seq_node_memory_;

#ifdef VELOUR_TBB
        tbb::atomic<size_t> seq_node_dead_memory_;
#else
        size_t seq_node_dead_memory_;
#endif // VELOUR_TBB

        SlabList used_list_;    // NOTE: contains current slab
};

class NodeAllocators : public NodeAllocatorsBase
{
    public:
        NodeAllocators(size_t max_memory) : flag_gc_needed_(false),
            kmer_node_allocator_(*this), seq_node_allocator_(*this),
            flag_kmer_gc_needed_(false), flag_seq_gc_needed_(false),
            flag_last_gc_check_below_watermark_(true),
            max_memory_(max_memory), allocated_memory_(0), active_memory_(0),
            freelist_memory_(0), peak_allocated_memory_(0)
        { Init(max_memory); }
            /*mmap_base_(NULL), mmap_size_(0), slab_base_(NULL), slab_next_(NULL), slab_end_(NULL)*/

        ~NodeAllocators();

        char * AllocateSlab(void);
        void AddSlabsToFreeList(SlabList::iterator iter_begin, SlabList::iterator iter_end); // called by child allocators
        // TODO: DeallocateFreeSlabs();

        void PrintStatus(FILE *out=stdout);

        KmerNodeAllocator& get_kmernode_allocator() { return kmer_node_allocator_; }
        SeqNodeAllocator& get_seqnode_allocator() { return seq_node_allocator_; }

        void GarbageCollect(void);

        size_t GetLiveMemory(void)
          { return kmer_node_allocator_.GetKmerLiveMemory() + seq_node_allocator_.GetSeqLiveMemory(); }

        size_t GetMaxSafeAllocation(void) // NOTE: includes free memory in freelists
          {
            if (allocated_memory_ >= max_memory_ || active_memory_ >= max_memory_) {
                return 0;
            } else {
                return ((max_memory_ >> 1) + (max_memory_ >> 2) + (max_memory_ >> 3)) - active_memory_; // 87% - active
            }
          }

        bool flag_gc_needed_;

    private:
        void Init(size_t max_memory);
        void UpdateGarbageCollectFlags(void);

// TODO    private:
    public:
#ifdef VELOUR_TBB
        tbb::queuing_mutex  node_allocators_lock_;
#endif

        KmerNodeAllocator kmer_node_allocator_;
        SeqNodeAllocator  seq_node_allocator_;

        bool flag_kmer_gc_needed_;
        bool flag_seq_gc_needed_;

        bool flag_last_gc_check_below_watermark_;

        size_t max_memory_;

        size_t allocated_memory_;   // = active + freelist
        size_t active_memory_;
        size_t freelist_memory_;

        size_t peak_allocated_memory_;

        /* // for pointer compression: contiguous region allocation
        void * mmap_base_;    // the pointer we later use to munmap the memory
        size_t mmap_size_;    // current size of mmap region

        char * slab_base_;      // first slab aligned for pointer compression
        char * slab_next_;      // next slab to distribute
        char * slab_end_;       // hard limit w.r.t. address space
        */

        SlabList free_list_;    // list of free slabs
};

#ifdef VELOUR_TBB
class TLS_KmerNodeAllocator : public NodeAllocatorsBase
{
    public:
        TLS_KmerNodeAllocator(KmerNodeAllocator& parent) :
            kmer_node_allocator_(parent), slab_current_(NULL), slab_end_(NULL) {}

        ~TLS_KmerNodeAllocator() {}

        KmerNode * AllocateNodeMemory(void);
        void DeallocateNodeMemory(KmerNode *);
        void PrepareForGarbageCollect(void); // explicitly mark unallocated nodes as dead, then drop pointers to slab

    private:
        KmerNodeAllocator& kmer_node_allocator_;

        char *slab_current_;
        char *slab_end_;
};

class TLS_SeqNodeAllocator : public NodeAllocatorsBase
{
    public:
        TLS_SeqNodeAllocator(SeqNodeAllocator& parent) :
            seq_node_allocator_(parent), slab_current_(NULL), slab_end_(NULL) {}

        ~TLS_SeqNodeAllocator() {}

        SeqNode * AllocateNodeMemory(size_t amount);
        SeqNode * ReallocateNodeMemory(SeqNode * old_node, size_t new_size);
        void DeallocateNodeMemory(SeqNode *);
        void PrepareForGarbageCollect(void); // explicitly mark unallocated nodes as dead, then drop pointers to slab

    private:
        SeqNodeAllocator& seq_node_allocator_;

        char *slab_current_;
        char *slab_end_;
};
#endif // VELOUR_TBB

//
//
// inlined member method implementations below
//
//

// explicitly mark unallocated nodes in slab as dead
static inline void KmerNodeSlab_SetUnallocatedNodesDead(char * current, char * end)
{
    for ( ; (current + sizeof(KmerNode)) <= end ; current += sizeof(KmerNode) ) {
        KmerNode *unallocated_node = reinterpret_cast<KmerNode*>(current);
        setNodeDead<KmerNode>(unallocated_node);
    }
}

// explicitly mark unallocated nodes in slab as dead  _AND_ set to zero length as well!
static inline void SeqNodeSlab_SetUnallocatedNodesDead(char * current, char * end)
{
    for ( ; (current + SeqNode::ComputeNodeAllocatedBytes(0)) <= end ;
            current += SeqNode::ComputeNodeAllocatedBytes(0)) {
        SeqNode *unallocated_node = reinterpret_cast<SeqNode*>(current);
        unallocated_node->sequence.zero_length();
        setNodeDead<SeqNode>(unallocated_node);
    }
}


#ifndef VELOUR_TBB
inline KmerNode * KmerNodeAllocator::AllocateNodeMemory(void)
{
    if ((slab_current_ + sizeof(KmerNode)) > slab_end_) {
        // not enough space left, need to grab a new slab
        //WASTED_SPACE += slab_end_ - slab_current_;
        char * new_slab = node_allocators_.AllocateSlab();
        slab_current_ = new_slab;
        slab_end_ = new_slab + SLAB_SIZE;

        used_list_.push_back(new_slab);
        kmer_node_memory_ += SLAB_SIZE;
    }
    assert( (slab_current_ + sizeof(KmerNode)) <= slab_end_ );

    KmerNode *retval = reinterpret_cast<KmerNode*>(slab_current_);
    slab_current_ += sizeof(KmerNode);
    return retval;
}

inline SeqNode * SeqNodeAllocator::AllocateNodeMemory(size_t amount)
{
    if ((slab_current_ + amount) > slab_end_) {
        // not enough space left, need to grab a new slab
        //WASTED_SPACE += slab_end_ - slab_current_;

        // but first, ensure the remaining space is dead
        SeqNodeSlab_SetUnallocatedNodesDead(slab_current_, slab_end_);

        char * new_slab = node_allocators_.AllocateSlab();
        slab_current_ = new_slab;
        slab_end_ = new_slab + SLAB_SIZE;

        used_list_.push_back(new_slab);
        seq_node_memory_ += SLAB_SIZE;
    }
    assert( (slab_current_ + amount) <= slab_end_ );

    SeqNode *retval = reinterpret_cast<SeqNode*>(slab_current_);
    slab_current_ += amount;
    return retval;
}

inline SeqNode * SeqNodeAllocator::ReallocateNodeMemory(SeqNode * old_node, size_t new_size)
{
    size_t old_size = old_node->GetNodeAllocatedBytes();
    assert( new_size > old_size ); // reallocation doesn't consider shrinking
    if ( ((reinterpret_cast<char*>(old_node) + old_size) == slab_current_) &&
         ((slab_current_ + (new_size - old_size)) <= slab_end_) ) { // can reallocate in-place
        slab_current_ += (new_size - old_size);
        return old_node;
    } else { // just allocate new node and copy
        SeqNode *new_node = AllocateNodeMemory(new_size);
        memcpy(new_node, old_node, old_size); // NOTE: does not zero new sequence memory!
        DeallocateNodeMemory(old_node);
        return new_node;
    }
}

inline void KmerNodeAllocator::DeallocateNodeMemory(KmerNode *node)
{
    /*if ((reinterpret_cast<char*>(node) + sizeof(KmerNode)) == slab_current_) { // can undo allocation
        slab_current_ = reinterpret_cast<char*>(node);
    } else*/ {
        setNodeDead<KmerNode>(node);
        kmer_node_dead_memory_ += sizeof(KmerNode);
    }
}

inline void SeqNodeAllocator::DeallocateNodeMemory(SeqNode *node)
{
    size_t node_size = node->GetNodeAllocatedBytes();
    /*if ((reinterpret_cast<char*>(node) + node_size) == slab_current_) { // can undo allocation
        slab_current_ = reinterpret_cast<char*>(node);
    } else*/ {
        setNodeDead<SeqNode>(node);
        seq_node_dead_memory_ += node_size;
    }
}
#endif // not VELOUR_TBB

#ifdef VELOUR_TBB
inline KmerNode * TLS_KmerNodeAllocator::AllocateNodeMemory(void)
{
    if ((slab_current_ + sizeof(KmerNode)) > slab_end_) {
        // not enough space left, need to grab a new slab
        //WASTED_SPACE += slab_end_ - slab_current_;
        char * new_slab = kmer_node_allocator_.AllocateTlsSlab();
        slab_current_ = new_slab;
        slab_end_ = new_slab + TLS_SLAB_SIZE;
    }
    assert( (slab_current_ + sizeof(KmerNode)) <= slab_end_ );

    KmerNode *retval = reinterpret_cast<KmerNode*>(slab_current_);
    slab_current_ += sizeof(KmerNode);
    return retval;
}

inline SeqNode * TLS_SeqNodeAllocator::AllocateNodeMemory(size_t amount)
{
    if ((slab_current_ + amount) > slab_end_) {
        // not enough space left, need to grab a new slab
        //WASTED_SPACE += slab_end_ - slab_current_;

        // but first, ensure the remaining space is dead
        SeqNodeSlab_SetUnallocatedNodesDead(slab_current_, slab_end_);

        char * new_slab = seq_node_allocator_.AllocateTlsSlab();
        slab_current_ = new_slab;
        slab_end_ = new_slab + TLS_SLAB_SIZE;
    }
    assert( (slab_current_ + amount) <= slab_end_ );

    SeqNode *retval = reinterpret_cast<SeqNode*>(slab_current_);
    slab_current_ += amount;
    return retval;
}

inline SeqNode * TLS_SeqNodeAllocator::ReallocateNodeMemory(SeqNode * old_node, size_t new_size)
{
    size_t old_size = old_node->GetNodeAllocatedBytes();
    assert( new_size > old_size ); // reallocation doesn't consider shrinking
    if ( ((reinterpret_cast<char*>(old_node) + old_size) == slab_current_) &&
         ((slab_current_ + (new_size - old_size)) <= slab_end_) ) { // can reallocate in-place
        slab_current_ += (new_size - old_size);
        return old_node;
    } else { // just allocate new node and copy
        SeqNode *new_node = AllocateNodeMemory(new_size);
        memcpy(new_node, old_node, old_size); // NOTE: does not zero new sequence memory!
        DeallocateNodeMemory(old_node);
        return new_node;
    }
}

inline void TLS_KmerNodeAllocator::DeallocateNodeMemory(KmerNode *node)
{
    /*if ((reinterpret_cast<char*>(node) + sizeof(KmerNode)) == slab_current_) { // can undo allocation
        slab_current_ = reinterpret_cast<char*>(node);
    } else*/ {
        setNodeDead<KmerNode>(node);
        kmer_node_allocator_.kmer_node_dead_memory_ += sizeof(KmerNode); // NOTE: safe atomic variable
    }
}

inline void TLS_SeqNodeAllocator::DeallocateNodeMemory(SeqNode *node)
{
    size_t node_size = node->GetNodeAllocatedBytes();
    /*if ((reinterpret_cast<char*>(node) + node_size) == slab_current_) { // can undo allocation
        slab_current_ = reinterpret_cast<char*>(node);
    } else*/ {
        setNodeDead<SeqNode>(node);
        seq_node_allocator_.seq_node_dead_memory_ += node_size; // NOTE: safe atomic variable
    }
}
#endif // VELOUR_TBB

#endif // VELOUR_NODE_ALLOCATORS_H
