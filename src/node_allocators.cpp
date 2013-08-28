//
// This file is part of Velour and is distributed under the University of
// Illinois Open Source License.  See LICENSE.txt for details.
//

//
// node_allocators.cpp
//

#include "types.h"
#include "node_allocators.h"

void NodeAllocators::Init(size_t max_memory)
{
    assert( is_power_of_2(SLAB_SIZE) );
#ifdef VELOUR_TBB
    assert( is_power_of_2(TLS_SLAB_SIZE) );
#endif
    
    /*
    // pointer compression:
    //   over-allocate virtual memory to guarantee the upper bits are constant
    size_t max_slab_memory = max_memory - (max_memory % SLAB_SIZE); // round to multiple of SLAB_SIZE
    size_t alloc_amt = 2 * max_slab_memory;

    mmap_base_ = mmap(0, alloc_amt, PROT_READ | PROT_WRITE, MAP_PRIVATE | MAP_NORESERVE | MAP_ANON, -1, 0);
    if (mmap_base_ == MAP_FAILED) {
        perror("ERROR: Velour failed to mmap allocate node virtual memory");
        exit(EXIT_FAILURE);
    }
    mmap_size_ = alloc_amt;

    size_t mask = ~(round_up_power_of_2(max_slab_memory) - 1); // FIXME
    char * masked_address = ( static_cast<char *>(malloc_base_) + max_slab_memory ) & mask;

    slab_next_ = slab_base_ = (masked_address > malloc_base_ ? masked_address : malloc_base_);
    slab_end_ = slab_base_ + max_slab_memory;

    assert( (slab_base_ & mask) == (slab_end_ & mask) );
    */
}
        
NodeAllocators::~NodeAllocators()
{
    /* pointer compression:
     if(mmap_base_ != NULL) munmap(mmap_base_, mmap_size_);*/

    //for (SlabList::iterator it = free_list_.begin(); it != free_list_.end(); ++it) { free( *it ); }

    printf("NodeAllocators: peak allocated memory %zu MB\n", peak_allocated_memory_ / (1024 * 1024));
}

char * NodeAllocators::AllocateSlab(void)
{
#ifdef VELOUR_TBB
    // lock the allocator -- automatically released upon return
    tbb::queuing_mutex::scoped_lock lock(node_allocators_lock_);
#endif

    // pull from the free list
    if (!free_list_.empty()) {
        char * slab = free_list_.front(); free_list_.pop_front();
        active_memory_ += SLAB_SIZE;
        freelist_memory_ -= SLAB_SIZE;
        assert( (active_memory_ + freelist_memory_) == allocated_memory_ );
        return slab;
    }
    
    // otherwise, we need to allocate a new slab

    /*// pointer compression: check not violating maximum allocation limit
    if (slab_next_ >= slab_end_) {
        fprintf(stderr, "ERROR: insufficient memory for graph nodes.\n");
        exit(EXIT_FAILURE);
    }*/

    // XXX: notify memory manager we are allocating memory: SLAB_SIZE

    /*// pointer compression: allocate a new slab
    char * slab = slab_next_;
    slab_next_ += SLAB_SIZE; */

    char * slab = static_cast<char *>( malloc(SLAB_SIZE) );
    if (slab == NULL) {
        perror("ERROR: node slab allocation failed");
        PrintStatus(stderr);
        exit(EXIT_FAILURE);
    }

    bool prior_overcommit_memory = (allocated_memory_ > max_memory_);
    bool prior_super_overcommit_memory = (allocated_memory_ > (1.2 * max_memory_));  // over 20%

    allocated_memory_ += SLAB_SIZE;
    active_memory_ += SLAB_SIZE;

    if (allocated_memory_ > peak_allocated_memory_) {
        peak_allocated_memory_ = allocated_memory_;
    }

    if ((allocated_memory_ > max_memory_) && !prior_overcommit_memory) {
        printf("WARNING: NodeAllocator memory overcommit: crossed %zu MB threshold.\n",
                max_memory_ / (1024 * 1024));
        PrintStatus();
        fflush(stdout);
    }

    if ((allocated_memory_ > (1.2 * max_memory_)) && !prior_super_overcommit_memory) {
        printf("WARNING: NodeAllocator memory overcommit over 20%% : crossed %zu MB threshold.\n",
                (size_t)((1.2 * max_memory_) / (1024 * 1024)));
        PrintStatus();
        fflush(stdout);
    }

    // SANITY: to avoid weird segfaults due to virtual memory overcommit,
    //   we should walk the pages to verify each will not segfault

    // update flags that may initiate garbage collection
    UpdateGarbageCollectFlags();

    return slab;
}

void NodeAllocators::AddSlabsToFreeList(SlabList::iterator iter_begin, SlabList::iterator iter_end)
{
#ifdef VELOUR_TBB
    // lock the allocator -- automatically released upon return
    tbb::queuing_mutex::scoped_lock lock(node_allocators_lock_);
#endif

    for (SlabList::iterator it = iter_begin ; it != iter_end ; ++ it) {
        char * free_slab = *( it );
        free_list_.push_back(free_slab);
        active_memory_ -= SLAB_SIZE;
        freelist_memory_ += SLAB_SIZE;
        assert( (active_memory_ + freelist_memory_) == allocated_memory_ );
    
        // XXX: notify memory manager we are free-list'ing memory: SLAB_SIZE
    }
}

void NodeAllocators::PrintStatus(FILE *out)
{
    fprintf(out, "  NodeAllocator state summary:\n");
    fprintf(out, "    %6zu MB set memory limit\n", max_memory_ / (1024*1024));
    fprintf(out, "    %6zu MB peak allocated memory\n", peak_allocated_memory_ / (1024*1024));

    fprintf(out, "\n");

    fprintf(out, "    %6zu MB allocated memory\n", allocated_memory_ / (1024*1024));
    fprintf(out, "    %6zu MB active memory\n", active_memory_ / (1024*1024));
    fprintf(out, "    %6zu MB freelist memory\n", freelist_memory_ / (1024*1024));
    fprintf(out, "    %6zu MB live memory\n",
        (active_memory_ - kmer_node_allocator_.kmer_node_dead_memory_ - seq_node_allocator_.seq_node_dead_memory_) /
        (1024*1024));
    
    fprintf(out, "\n");

    fprintf(out, "    %6zu MB live kmer memory = %6zu MB active - %6zu MB dead\n",
        (kmer_node_allocator_.kmer_node_memory_ - kmer_node_allocator_.kmer_node_dead_memory_) / (1024*1024),
        (kmer_node_allocator_.kmer_node_memory_) / (1024*1024),
        (kmer_node_allocator_.kmer_node_dead_memory_) / (1024*1024) );

    fprintf(out, "    %6zu MB live  seq memory = %6zu MB active - %6zu MB dead\n",
        (seq_node_allocator_.seq_node_memory_ - seq_node_allocator_.seq_node_dead_memory_) / (1024*1024),
        (seq_node_allocator_.seq_node_memory_) / (1024*1024),
        (seq_node_allocator_.seq_node_dead_memory_) / (1024*1024) );
    fprintf(out, "  END NodeAllocator state summary.\n");
}

#ifdef VELOUR_TBB
char * KmerNodeAllocator::AllocateTlsSlab(void)
{
    // lock the allocator -- automatically released upon return
    tbb::queuing_mutex::scoped_lock lock(kmer_node_allocator_lock_);

    if ((slab_current_ + TLS_SLAB_SIZE) > slab_end_) {
        // not enough space left, need to grab a new slab
        //WASTED_SPACE += slab_end_ - slab_current_;
        char * new_slab = node_allocators_.AllocateSlab();
        slab_current_ = new_slab;
        slab_end_ = new_slab + SLAB_SIZE;

        used_list_.push_back(new_slab);
        kmer_node_memory_ += SLAB_SIZE;
    }

    char * slab = slab_current_;
    slab_current_ += TLS_SLAB_SIZE;
    return slab;
}

char * SeqNodeAllocator::AllocateTlsSlab(void)
{
    // lock the allocator -- automatically released upon return
    tbb::queuing_mutex::scoped_lock lock(seq_node_allocator_lock_);

    if ((slab_current_ + TLS_SLAB_SIZE) > slab_end_) {
        // not enough space left, need to grab a new slab
        //WASTED_SPACE += slab_end_ - slab_current_;
        char * new_slab = node_allocators_.AllocateSlab();
        slab_current_ = new_slab;
        slab_end_ = new_slab + SLAB_SIZE;

        used_list_.push_back(new_slab);
        seq_node_memory_ += SLAB_SIZE;
    }

    char * slab = slab_current_;
    slab_current_ += TLS_SLAB_SIZE;
    return slab;
}
#endif // VELOUR_TBB

#ifdef VELOUR_TBB
TLS_KmerNodeAllocator * KmerNodeAllocator::GetTlsAllocator(threadid_t tid)
{
    if (thread_allocators_[tid] == NULL) {
        thread_allocators_[tid] = new TLS_KmerNodeAllocator(*this);
    }

    return thread_allocators_[tid];
}

TLS_SeqNodeAllocator * SeqNodeAllocator::GetTlsAllocator(threadid_t tid)
{
    if (thread_allocators_[tid] == NULL) {
        thread_allocators_[tid] = new TLS_SeqNodeAllocator(*this);
    }

    return thread_allocators_[tid];
}
#endif // VELOUR_TBB

//
// garbage collection
//

void NodeAllocators::UpdateGarbageCollectFlags(void)
{
    bool old_flag_gc_needed = flag_gc_needed_;
        
    assert( (active_memory_ + freelist_memory_) == allocated_memory_ );

    // XXX: counts unallocated slab memory as live
    size_t live_memory =
        (active_memory_ - kmer_node_allocator_.kmer_node_dead_memory_ - seq_node_allocator_.seq_node_dead_memory_);

    // condition: active memory > about 87% of max memory limit
    bool active_memory_above_watermark =
        (active_memory_ >= ((max_memory_ >> 1) + (max_memory_ >> 2) + (max_memory_ >> 3)));

    // condition: active minus dead <= about 97% of max memory limit
    bool live_memory_below_limit = (live_memory <= ((max_memory_ >> 1) + (max_memory_ >> 2) + (max_memory_ >> 3) +
             (max_memory_ >> 4) + (max_memory_ >> 5)));

    if (active_memory_above_watermark && live_memory_below_limit)
    {
        // kmer allocator condition:
        //   dead kmer node memory > about 3% of kmer memory
        flag_kmer_gc_needed_ =
            (kmer_node_allocator_.kmer_node_dead_memory_ >= (kmer_node_allocator_.kmer_node_memory_ >> 5));

        // seq allocator condition:
        //   dead seq node memory > about 3% of seq memory
        flag_seq_gc_needed_ =
            (seq_node_allocator_.seq_node_dead_memory_ >= (seq_node_allocator_.seq_node_memory_ >> 5));

        // update single flag tested by external functions
        flag_gc_needed_ |= flag_kmer_gc_needed_ || flag_seq_gc_needed_;
    } /* else {
        assert( !flag_gc_needed_ && !flag_kmer_gc_needed_ && !flag_seq_gc_needed_ ); TODO ???
    }*/

    // if just decided to gc, print status
    if (!old_flag_gc_needed && flag_gc_needed_) {
        printf("  NodeAllocator: garbage collection need identified.\n");
        PrintStatus();
    }
}

void KmerNodeAllocator::PrepareForGarbageCollect(void)
{
#ifdef VELOUR_TBB
    for (tbb::concurrent_vector<TLS_KmerNodeAllocator*>::iterator it = thread_allocators_.begin();
            it != thread_allocators_.end(); ++it) {
        if (*it != NULL) {
            (*it)->PrepareForGarbageCollect();
        }
    }
#endif // VELOUR_TBB

    KmerNodeSlab_SetUnallocatedNodesDead(slab_current_, slab_end_);

    // drop reference to slab
    slab_current_ = NULL;
    slab_end_ = NULL;
}

void SeqNodeAllocator::PrepareForGarbageCollect(void)
{
#ifdef VELOUR_TBB
    for (tbb::concurrent_vector<TLS_SeqNodeAllocator*>::iterator it = thread_allocators_.begin();
            it != thread_allocators_.end(); ++it) {
        if (*it != NULL) {
            (*it)->PrepareForGarbageCollect();
        }
    }
#endif // VELOUR_TBB

    SeqNodeSlab_SetUnallocatedNodesDead(slab_current_, slab_end_);

    // drop reference to slab
    slab_current_ = NULL;
    slab_end_ = NULL;
}

#ifdef VELOUR_TBB
void TLS_KmerNodeAllocator::PrepareForGarbageCollect(void)
{
    KmerNodeSlab_SetUnallocatedNodesDead(slab_current_, slab_end_);

    // drop reference to TLS slab
    slab_current_ = NULL;
    slab_end_ = NULL;
}

void TLS_SeqNodeAllocator::PrepareForGarbageCollect(void)
{
    SeqNodeSlab_SetUnallocatedNodesDead(slab_current_, slab_end_);

    // drop reference to TLS slab
    slab_current_ = NULL;
    slab_end_ = NULL;
}
#endif // VELOUR_TBB

void NodeAllocators::GarbageCollect(void)
{
    assert( flag_gc_needed_ && (flag_kmer_gc_needed_ || flag_seq_gc_needed_) );

    if (flag_kmer_gc_needed_) {
        if (kmer_node_allocator_.kmer_graph_ != NULL) {
            kmer_node_allocator_.kmer_graph_->PrepareForMajorGarbageCollect();
            kmer_node_allocator_.PrepareForGarbageCollect();
            kmer_node_allocator_.MajorGarbageCollect();
        } else {
            kmer_node_allocator_.BulkFreeListAllSlabs();
        }
        flag_kmer_gc_needed_ = false;
    }

    if (flag_seq_gc_needed_) {
        if (seq_node_allocator_.seq_graph_ != NULL) {
            seq_node_allocator_.seq_graph_->PrepareForMajorGarbageCollect();
            seq_node_allocator_.PrepareForGarbageCollect();
            seq_node_allocator_.MajorGarbageCollect();
        } else {
            seq_node_allocator_.BulkFreeListAllSlabs();
        }
        flag_seq_gc_needed_ = false;
    }

    flag_gc_needed_ = false;
}

// major garbage collection: compact all nodes and rebuild hashtable
void KmerNodeAllocator::MajorGarbageCollect(void)
{
    assert( kmer_graph_ != NULL );
    assert( kmer_graph_->node_count == 0 );

    // ensures that PrepareForGarbageCollect was called
    assert( slab_current_ == NULL );
    assert( slab_end_ == NULL );

    if (used_list_.empty()) { return; }

    printf("KmerNodeAllocator: major garbage collection begins: %zu MB / %zu MB is dead\n",
            (kmer_node_dead_memory_ / (1024*1024)), (kmer_node_memory_ / (1024*1024)) ); fflush(stdout);

    size_t old_kmer_node_memory = kmer_node_memory_;

    SlabList::iterator fill_iter = used_list_.begin();
    char * fill_slab = *( fill_iter );
    char * fill_current = fill_slab;
    char * fill_end = fill_slab + SLAB_SIZE;

    for (SlabList::iterator load_iter = used_list_.begin(); load_iter != used_list_.end() ; ++ load_iter) {
        char * load_slab = *( load_iter );
        char * load_current = load_slab; 
        char * load_end = load_slab + SLAB_SIZE;

        for ( ; (load_current + sizeof(KmerNode)) <= load_end ; load_current += sizeof(KmerNode) ) {
            KmerNode *load_node = reinterpret_cast<KmerNode*>(load_current);
            if (!isNodeDead<KmerNode>(load_node)) {
                KmerNode *fill_node = reinterpret_cast<KmerNode*>(fill_current);
                if (fill_node != load_node) { // don't self copy
                    *fill_node = *load_node; // copy node contents
                }
                fill_current += sizeof(KmerNode);

                kmer_graph_->insertNode( fill_node );
            }

            if ((fill_current + sizeof(KmerNode)) > fill_end) { // bump the fill iterator if full
                ++ fill_iter;
                assert( fill_iter != used_list_.end() );

                fill_slab = *( fill_iter );
                fill_current = fill_slab;
                fill_end = fill_slab + SLAB_SIZE;
            }
        }
    }

    slab_current_ = fill_current;
    slab_end_ = fill_end;

    SlabList::iterator free_begin = fill_iter + 1;
    SlabList::iterator free_end = used_list_.end();

    if (free_begin != free_end) {
        node_allocators_.AddSlabsToFreeList(free_begin, free_end);
        kmer_node_memory_ -= (free_end - free_begin) * SLAB_SIZE;
        used_list_.erase(free_begin, free_end);
    }

    kmer_node_dead_memory_ = 0;

    printf("KmerNodeAllocator: major garbage collection finished: reclaimed %zu MB\n",
            ((old_kmer_node_memory - kmer_node_memory_) / (1024*1024)) ); fflush(stdout);
}

// major garbage collection: compact all nodes and rebuild hashtables
void SeqNodeAllocator::MajorGarbageCollect(void)
{
    assert( seq_graph_ != NULL );
    assert( seq_graph_->node_count == 0 );

    // ensures that PrepareForGarbageCollect was called
    assert( slab_current_ == NULL );
    assert( slab_end_ == NULL );

    if (used_list_.empty()) { return; }

    printf("SeqNodeAllocator: major garbage collection begins: %zu MB / %zu MB is dead\n",
            (seq_node_dead_memory_ / (1024*1024)), (seq_node_memory_ / (1024*1024)) ); fflush(stdout);

    size_t old_seq_node_memory = seq_node_memory_;

    seq_node_dead_memory_ = 0;

    SlabList::iterator fill_iter = used_list_.begin();
    char * fill_slab = *( fill_iter );
    char * fill_current = fill_slab;
    char * fill_end = fill_slab + SLAB_SIZE;

    for (SlabList::iterator load_iter = used_list_.begin(); load_iter != used_list_.end() ; ++ load_iter) {
        char * slab_start = *( load_iter );
        char * slab_current = slab_start;
        char * slab_end = slab_start + SLAB_SIZE;

#ifdef VELOUR_TBB
        // handle TLS slab's wasted space when walking slabs
        for ( ; slab_current < slab_end ; slab_current += TLS_SLAB_SIZE ) {
            assert( slab_current + TLS_SLAB_SIZE <= slab_end );
#endif

        char * load_slab = slab_current;
        char * load_current = load_slab;
#ifdef VELOUR_TBB
        char * load_end = load_slab + TLS_SLAB_SIZE;
#else
        char * load_end = load_slab + SLAB_SIZE;
#endif

        for ( ; (load_current + SeqNode::ComputeNodeAllocatedBytes(0) <= load_end) ;
                load_current += reinterpret_cast<SeqNode*>(load_current)->GetNodeAllocatedBytes() )
        {
            assert( load_current + reinterpret_cast<SeqNode*>(load_current)->GetNodeAllocatedBytes() <= load_end );
            SeqNode *load_node = reinterpret_cast<SeqNode*>(load_current);
            if (!isNodeDead<SeqNode>(load_node)) {
                size_t node_size = load_node->GetNodeAllocatedBytes();
                if (fill_current + node_size > fill_end) { // pad wasted space, bump the fill iterator
                    SeqNodeSlab_SetUnallocatedNodesDead(fill_current, fill_end);

                    ++ fill_iter;
                    assert( fill_iter != used_list_.end() );

                    fill_slab = *( fill_iter );
                    fill_current = fill_slab;
                    fill_end = fill_slab + SLAB_SIZE;
                }
                SeqNode *fill_node = reinterpret_cast<SeqNode*>(fill_current);
                if (fill_node != load_node) { // don't self copy
                    memmove(fill_node, load_node, node_size); // copy node contents
                }
                fill_current += node_size;

                seq_graph_->insertNode( fill_node );
            }
        }

#ifdef VELOUR_TBB
        }
#endif

    }

#ifdef VELOUR_TBB
    // align for future TLS allocations
    if (fill_current + TLS_SLAB_SIZE < fill_end) { // one or more TLS slabs fit in current slab
        char * tls_end = fill_current + (TLS_SLAB_SIZE - ((fill_current - fill_slab) % TLS_SLAB_SIZE));
        SeqNodeSlab_SetUnallocatedNodesDead(fill_current, tls_end);

        slab_current_ = tls_end;
        slab_end_ = fill_end;
    } else {
        SeqNodeSlab_SetUnallocatedNodesDead(fill_current, fill_end);

        slab_current_ = NULL;
        slab_end_ = NULL;
    }
#else
    slab_current_ = fill_current;
    slab_end_ = fill_end;
#endif

    SlabList::iterator free_begin = fill_iter + 1;
    SlabList::iterator free_end = used_list_.end();

    if (free_begin != free_end) {
        node_allocators_.AddSlabsToFreeList(free_begin, free_end);
        seq_node_memory_ -= (free_end - free_begin) * SLAB_SIZE;
        used_list_.erase(free_begin, free_end);
    }

    printf("SeqNodeAllocator: major garbage collection finished: reclaimed %zu MB\n",
            ((old_seq_node_memory - seq_node_memory_) / (1024*1024)) ); fflush(stdout);
}

/*// minor garbage collection: compact all nodes in-place and fix up hashtable
void KmerNodeAllocator::MinorGarbageCollect(void)
{
    // TODO
}*/

/*// minor garbage collection: compact all nodes in-place and fix up hashtable
void SeqNodeAllocator::MinorGarbageCollect(void)
{
    // TODO
}*/


void KmerNodeAllocator::BulkFreeListAllSlabs(void)
{
    node_allocators_.AddSlabsToFreeList(used_list_.begin(), used_list_.end());
    kmer_node_memory_ = 0;
    kmer_node_dead_memory_ = 0;
    used_list_.clear();
}

void SeqNodeAllocator::BulkFreeListAllSlabs(void)
{
    node_allocators_.AddSlabsToFreeList(used_list_.begin(), used_list_.end());
    seq_node_memory_ = 0;
    seq_node_dead_memory_ = 0;
    used_list_.clear();
}

