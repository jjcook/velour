//
// This file is part of Velour and is distributed under the University of
// Illinois Open Source License.  See LICENSE.txt for details.
//

//
// flowing.cpp
//

#include "types.h"

static uintptr_t p__peakKmerNodes = 0;
static uintptr_t p__peakSeqNodes = 0;

static size_t p__peakKmerLiveMemory = 0;
static size_t p__peakSeqLiveMemory = 0;
size_t g__peakLiveMemory = 0;

//static size_t p__peakFragmentNodes = 0;
//static size_t p__peakFragmentMemory = 0;

static uintptr_t p__preRedistributionSequenceNodes = 0;
//static size_t p__preRedistributionSequenceNodeMemory = 0;

#ifdef VELOUR_MPI
static size_t *p__MPI_SAFE_LENGTH_ARRAY;
MPI_Win g__MPI_WINDOW_SAFE_LENGTHS;
#endif // VELOUR_MPI

template<typename GraphType, typename NodeType, typename NodeList>
static uintptr_t SerialMemoizeComponent(GraphType *graph, NodeType *root, NodeList *nodelist)
{
    uintptr_t nodes_memoized = 0;

    if (!isNodeDead<NodeType>(root) && root->claim_tid == 0) {
        std::deque<NodeType *> worklist;

        root->claim_tid = 1;
        worklist.push_front(root);
        nodelist->push_back(root);
        while(!worklist.empty()) {
            NodeType *node = worklist.back();
            worklist.pop_back();

            // check right side
            counter_t *counters = node->right_count;
            Color *colors = node->right_color;
            for (int i = 0; i < 4; ++ i) {
                if (counters[i] != 0) {
                    if (colors[i] == 0) {
                        NodeType *next = graph->findNextNode(node, i, GO_RIGHT);
                        assert( !isNodeDead<NodeType>(next) );
                        if (next->claim_tid == 0) {
                            next->claim_tid = 1;
                            worklist.push_front(next);
                            nodelist->push_back(next);
                            ++ nodes_memoized;
                        }
                    }
                }
            }

            // check left side
            counters = node->left_count;
            colors = node->left_color;
            for (int i = 0; i < 4; ++ i) {
                if (counters[i] != 0) {
                    if (colors[i] == 0) {
                        NodeType *next = graph->findNextNode(node, i, GO_LEFT);
                        assert( !isNodeDead<NodeType>(next) );
                        if (next->claim_tid == 0) {
                            next->claim_tid = 1;
                            worklist.push_front(next);
                            nodelist->push_back(next);
                            ++ nodes_memoized;
                        }
                    }
                }
            }
        }
    }

    return nodes_memoized;
}

#ifdef VELOUR_TBB
static inline bool coopClaimNode(SeqNode *node)
{
    bool claimed = (ATOMIZE(node->claim_tid).compare_and_swap(1, 0) == 0);
    return claimed;
}

template<typename GraphType, typename NodeType, typename NodeList>
static uintptr_t ParallelCoopMemoizeComponent(GraphType *graph, NodeType *root, NodeList *nodelist)
{
    uintptr_t nodes_memoized = 0;

    if (!atomic_isNodeDead<NodeType>(root) && coopClaimNode(root)) {
        std::deque<NodeType *> worklist;

        worklist.push_front(root);
        nodelist->push_back(root);
        while(!worklist.empty()) {
            NodeType *node = worklist.back();
            worklist.pop_back();

            // check right side
            counter_t *counters = node->right_count;
            Color *colors = node->right_color;
            for (int i = 0; i < 4; ++ i) {
                if (counters[i] != 0) {
                    if (colors[i] == 0) {
                        NodeType *next = graph->findNextNode(node, i, GO_RIGHT);
                        assert( !atomic_isNodeDead<NodeType>(next) );
                        if (coopClaimNode(next)) {
                            worklist.push_front(next);
                            nodelist->push_back(next);
                            ++ nodes_memoized;
                        }
                    }
                }
            }

            // check left side
            counters = node->left_count;
            colors = node->left_color;
            for (int i = 0; i < 4; ++ i) {
                if (counters[i] != 0) {
                    if (colors[i] == 0) {
                        NodeType *next = graph->findNextNode(node, i, GO_LEFT);
                        assert( !atomic_isNodeDead<NodeType>(next) );
                        if (coopClaimNode(next)) {
                            worklist.push_front(next);
                            nodelist->push_back(next);
                            ++ nodes_memoized;
                        }
                    }
                }
            }
        }
    }

    return nodes_memoized;
}
#endif // VELOUR_TBB

template<typename NodeList>
uintptr_t sg_serial_chain_memoize_related_components(SeqNode *seed_chain, SeqGraph *source_graph, NodeList *nodelist)
{
    uintptr_t nodes_memoized = 0;

    SeqNode *iterator = seed_chain;
    while (iterator != NULL) {
        SeqNode *root = iterator;
        for (int i = 0 ; i < 4 ; ++ i) {
            // check on the left side
            if (root->left_color[i] == g__PARTITION_INDEX) {
                assert( root->left_count[i] != 0 );
                SeqNode *next = source_graph->findNextNode(root, i, GO_LEFT);
                if (next != NULL) {
                    nodes_memoized += SerialMemoizeComponent<SeqGraph, SeqNode, NodeList>(source_graph, next, nodelist);
                }
            }
            // check on the right side
            if (root->right_color[i] == g__PARTITION_INDEX) {
                assert( root->right_count[i] != 0 );
                SeqNode *next = source_graph->findNextNode(root, i, GO_RIGHT);
                if (next != NULL) {
                    nodes_memoized += SerialMemoizeComponent<SeqGraph, SeqNode, NodeList>(source_graph, next, nodelist);
                }
            }
        }
        iterator = iterator->head_next;
    }

    return nodes_memoized;
}

#ifdef VELOUR_TBB
namespace {

struct lambda_parallel_memoize_components {
    typedef tbb::concurrent_vector<SeqNode *>::const_range_type::iterator iterator;
    lambda_parallel_memoize_components(SeqGraph *source_graph, flow_nodelist_t *nodelist, tbb::atomic<uintptr_t> &nodes_memoized) :
        source_graph(source_graph), nodelist(nodelist), nodes_memoized(nodes_memoized) {}
    void operator()(const tbb::concurrent_vector<SeqNode *>::const_range_type& range) const {
        for (iterator it = range.begin() ; it != range.end() ; ++it) {
            SeqNode *root = *it;
            for (int i = 0 ; i < 4 ; ++ i) {
                // check on the left side
                if (root->left_color[i] == g__PARTITION_INDEX) {
                    assert( root->left_count[i] != 0 );
                    SeqNode *next = source_graph->findNextNode(root, i, GO_LEFT);
                    if (next != NULL) {
                        nodes_memoized += ParallelCoopMemoizeComponent<SeqGraph, SeqNode, flow_nodelist_t>(source_graph, next, nodelist);
                    }
                }
                // check on the right side
                if (root->right_color[i] == g__PARTITION_INDEX) {
                    assert( root->right_count[i] != 0 );
                    SeqNode *next = source_graph->findNextNode(root, i, GO_RIGHT);
                    if (next != NULL) {
                        nodes_memoized += ParallelCoopMemoizeComponent<SeqGraph, SeqNode, flow_nodelist_t>(source_graph, next, nodelist);
                    }
                }
            }
        }
    }
    private:
        SeqGraph *source_graph;
        flow_nodelist_t *nodelist;
        tbb::atomic<uintptr_t> &nodes_memoized;
};

} // namespace: anonymous

uintptr_t sg_parallel_chain_memoize_related_components(SeqNode *seed_chain, SeqGraph *source_graph, flow_nodelist_t *nodelist)
{
    tbb::atomic<uintptr_t> nodes_memoized;
    nodes_memoized = 0;

    flow_nodelist_t input_nodelist;

    // XXX HACK FOR PARALLELISM: build vector that can be processed by parallel_for
    SeqNode *iterator = seed_chain;
    while (iterator != NULL) {
        input_nodelist.push_back(iterator);
        iterator = iterator->head_next;
    }

    // do it in parallel
    //tbb::parallel_for(tbb::blocked_range<size_t>(0,input_nodelist.size(),1), // XXX: configurable grain size
    tbb::parallel_for(input_nodelist.range(), // XXX: configurable grain size
            lambda_parallel_memoize_components(source_graph, nodelist, nodes_memoized), tbb::auto_partitioner());

    return nodes_memoized;
}
#endif // VELOUR_TBB

static void flow_private_graph(SeqGraph *private_sgraph, SeqGraph *resident_sgraph, SplitBuckets *buckets, uintptr_t round)
{
    // TODO: histogram of initial private_sgraph size, etc etc

    /*if (g__DOTGRAPH) {
      char dot_filename[PATH_MAX+1];
      sprintf(dot_filename, "%s/SplitBucket-%u-%u-resident.dot", g__WORK_BASE_DIRECTORY, g__PARTITION_INDEX, round);
      emit_graphviz(resident_sgraph, dot_filename);

      char dot_filename2[PATH_MAX+1];
      sprintf(dot_filename2, "%s/SplitBucket-%u-%u-work.dot", g__WORK_BASE_DIRECTORY, g__PARTITION_INDEX, round);
      emit_graphviz(private_sgraph, dot_filename2);
      }*/

    /*if( g__FULL_STATISTICS ) {
      sg_stat_components(private_sgraph, stdout);
      }*/

#ifdef VERIFY_THREADSAFE
    private_sgraph->verify(true);
#endif

    sg_remove_tips(private_sgraph, true);

#ifdef VERIFY_THREADSAFE
    private_sgraph->verify(true);
#endif

    sg_concatenate(private_sgraph, true);

#ifdef VERIFY_THREADSAFE
    private_sgraph->verify(true);
#endif

    ATOMIC_ADD(p__preRedistributionSequenceNodes,private_sgraph->node_count);

    /*
    struct lambda1 {
        static void compute_sequence_sizes(SeqNode *node) {
            if (!isNodeDead<SeqNode>(node)) {
                todo_allocated_size += node->sequence.GetAllocatedBytes();
                todo_actual_size += node->sequence.GetLengthInBytes();
            }
        }
    };
    sg_for_each(private_sgraph, lambda1::compute_sequence_sizes);
    p__preRedistributionSequenceNodeMemory += ...;
    */


    /*if (g__DOTGRAPH) {
      char dot_filename[PATH_MAX+1];
      sprintf(dot_filename, "%s/SplitBucket-%u-%u-simplify.dot", g__WORK_BASE_DIRECTORY, currentPartitionIndex, round);
      emit_graphviz(private_sgraph, dot_filename);
      }*/

    // emit sub-components that are no longer relevant to
    if (g__SLICING) { slice2_graph(private_sgraph, g__PARTITION_INDEX); }

    /*if (g__DOTGRAPH) {
      char dot_filename[PATH_MAX+1];
      sprintf(dot_filename, "%s/SplitBucket-%u-%u-worksliced.dot", g__WORK_BASE_DIRECTORY, g__PARTITION_INDEX, round);
      emit_graphviz(private_sgraph, dot_filename);
      }*/

    // emit components that are no longer relevant to the working graph
    buckets->split(private_sgraph);

    /*if (g__DOTGRAPH) {
      char dot_filename[PATH_MAX+1];
      sprintf(dot_filename, "%s/SplitBucket-%u-%u-resident-postemit.dot", g__WORK_BASE_DIRECTORY, g__PARTITION_INDEX, round);
      emit_graphviz(resident_sgraph, dot_filename);

      char dot_filename2[PATH_MAX+1];
      sprintf(dot_filename2, "%s/SplitBucket-%u-%u-work-postemit.dot", g__WORK_BASE_DIRECTORY, g__PARTITION_INDEX, round);
      emit_graphviz(private_sgraph, dot_filename2);
      }*/

    // move remaining nodes from working graph into resident graph
    private_sgraph->bulkMoveAllNodes(resident_sgraph);
    assert( private_sgraph->node_count == 0 );
}

namespace {
    // build sequence components from kmer components, split & slice, remainder dumped into resident graph
    // TODO: update peak live memory etc etc in here
    struct lambda_flowKmerComponent {
        lambda_flowKmerComponent(KmerGraph *kgraph, SplitBuckets *buckets, SeqGraph *resident_sgraph) :
            kgraph_(kgraph), buckets_(buckets), resident_sgraph_(resident_sgraph) {}

        void operator()(KmerNode *seed) {
            if (isNodeDead<KmerNode>(seed)) return; // HACK

#ifdef VELOUR_TBB
            ParallelComponent<KmerGraph, KmerNode> component(kgraph_, seed);
#  ifdef VERIFY_THREADSAFE
            // VERIFY we own the whole component
            for(std::deque<KmerNode*>::iterator it=component.nodes_.begin(); it != component.nodes_.end(); ++it) {
                assert( (*it)->claim_tid == tls_thread_index );
            }
#  endif // VERIFY_THREADSAFE
#else
            SerialComponent<KmerGraph, KmerNode> component(kgraph_, seed);
#endif
            if (!component.nodes_.empty()) {

#ifdef VELOUR_TBB
                SeqGraph * private_sgraph = tls_private_sgraph_set->GetPrivateGraph(component.nodes_.size());
#else
                SeqGraph * private_sgraph = pgset.GetPrivateGraph(component.nodes_.size());
#endif // VELOUR_TBB

                assert( private_sgraph->node_count == 0 );
                sg_build_component(kgraph_, component.nodes_, private_sgraph);
                flow_private_graph(private_sgraph, resident_sgraph_, buckets_, 0);
                assert( private_sgraph->node_count == 0 );
            }
        }
        private:
            KmerGraph *kgraph_;
            SplitBuckets *buckets_;
            SeqGraph *resident_sgraph_;
#ifndef VELOUR_TBB
            PrivateSeqGraphSet pgset;
#endif
    };
} // namespace: anonymous

void newflow_inbox(file_object_vector &inbox_file_objects, SeqGraph *resident_sgraph, SplitBuckets *buckets,
        uintptr_t &total_bucket_nodes)
{

#ifdef VELOUR_TBB
    tbb::tick_count::interval_t time_inboxload, time_relatedcomponents, time_insertnodes;
    tbb::tick_count::interval_t time_removetips, time_concatenate;
    tbb::tick_count::interval_t time_slicing, time_splitting, time_cleanup;
#endif

    if (!inbox_file_objects.empty()) { // no inbox for the first partition

        // incrementally load components from the inbox files
        for (file_object_vector::iterator itr=inbox_file_objects.begin() ; itr != inbox_file_objects.end() ; ++itr) {
            if (itr->filetype != BUCKET) {
                fprintf(stderr, "ERROR: Cannot use %s file %s as flow input.  Exiting...\n",
                    FILE_TYPES[itr->filetype], itr->filename);
                exit(EXIT_FAILURE);
            }

#ifdef VELOUR_MPI
            unsigned inbox_partition = itr->fileindex;
#endif // VELOUR_MPI

            printf("    flowing inbox: %s\n", itr->filename); fflush(stdout);

            uintptr_t round = 1;
            size_t file_offset = 0;

#ifdef VELOUR_MPI
            int producer_finished = 0;
            size_t safe_length = 0;
            do {
            //if (safe_length < p__MPI_SAFE_LENGTH_ARRAY[inbox_partition])
            //    printf("DBG: Partition %u read new safe_length %lli\n", g__PARTITION_INDEX, safe_length); fflush(stdout);
            size_t new_safe_length = * static_cast<volatile size_t *>(&p__MPI_SAFE_LENGTH_ARRAY[inbox_partition]);
            if (new_safe_length == safe_length) {
                MPI_Status mpi_status;
                MPI_Iprobe(PART_TO_RANK(inbox_partition), 42, MPI_COMM_WORLD, &producer_finished, &mpi_status);
                // TODO: do MPI_Recv() on the message()
                continue;
            } else {
                safe_length = new_safe_length;
            }
#endif // VELOUR_MPI

            int filedes = open(itr->filename, O_RDONLY);
            if (filedes == -1) {
                fprintf(stderr, "ERROR: failed to open file: %s\n", itr->filename);
                perror("REASON: ");
                exit(EXIT_FAILURE);
            }

            struct stat file_stat;
            if (fstat(filedes, &file_stat) != 0) {
                fprintf(stderr,"ERROR: failed to stat file: %s\n", itr->filename);
                exit(EXIT_FAILURE);
            }

            size_t file_length = file_stat.st_size;
            //itr->length = file_stat.st_size;

#ifdef VELOUR_MPI
            file_length = min(file_length, safe_length);
#endif

            if (file_length == 0) { // otherwise, mmap will fail
                if (close(filedes) == -1) {
                    fprintf(stderr, "ERROR: failed to close file: %s\n", itr->filename);
                    perror("REASON: ");
                    exit(EXIT_FAILURE);
                }
#ifdef VELOUR_MPI
                MPI_Status mpi_status;
                MPI_Iprobe(PART_TO_RANK(inbox_partition), 42, MPI_COMM_WORLD, &producer_finished, &mpi_status);
                if (producer_finished && safe_length == * static_cast<volatile size_t *>(&p__MPI_SAFE_LENGTH_ARRAY[inbox_partition])) {
                    // TODO: do MPI_Recv() on the message()
                    //printf("DBG: Partition %u producer finished break\n", g__PARTITION_INDEX); fflush(stdout);
                    break; // break out of do loop to start next file
                } else {
                    //printf("DBG: Partition %u not finished continue\n", g__PARTITION_INDEX); fflush(stdout);
                    continue; // continue do loop
                }
#else
                continue; // continue to start next file
#endif // VELOUR_MPI
            }

            char *file_mmap = static_cast<char *>( mmap( 0, file_length, PROT_READ, MAP_PRIVATE, filedes, 0) );
            if (file_mmap == reinterpret_cast<char *>(-1)) {
                fprintf(stderr, "ERROR: failed to mmap file: %s\n", itr->filename);
                perror("REASON: ");
                exit(EXIT_FAILURE);
            }

            flow_nodelist_t flowlist;
#ifdef VELOUR_TBB
            flowlist.reserve(1000000); // XXX: constant, should relate to value below
#endif

            while (file_offset < file_length) {
                uintptr_t round_inbox_nodes = 0;
                uintptr_t round_import_nodes = 0;

#ifdef VELOUR_TBB
                tbb::tick_count time0 = tbb::tick_count::now();
#endif
                SeqNode *inbox_chain = NULL;
                uintptr_t inbox_chain_length = 0;
                uintptr_t min_inbox_chain_length = max(100000, resident_sgraph->node_count >> 5);  // about 3%

                //
                // CHOOSE WHICH MODE: PERFORMANCE or MINIMUM PEAK FOOTPRINT MODE FOR PUBLICATIONS
                //
                uintptr_t max_inbox_chain_length;
                if (g__MINIMIZE_FOOTPRINT) {
                    max_inbox_chain_length = min_inbox_chain_length;
                } else {
                    max_inbox_chain_length = max(min_inbox_chain_length, (g__NODE_ALLOCATORS->GetMaxSafeAllocation() / 100)); // XXX: max average seqnode size constant
                }

                while (file_offset < file_length && inbox_chain_length < max_inbox_chain_length) {
                    size_t component_size;
                    uintptr_t component_length = sg_load_mmap_stream_component((file_mmap+file_offset), &component_size, &inbox_chain);
                    file_offset += component_size;
                    inbox_chain_length += component_length;
                }

                total_bucket_nodes += inbox_chain_length;
                round_inbox_nodes += inbox_chain_length;

                p__peakSeqNodes = max(p__peakSeqNodes, resident_sgraph->node_count + inbox_chain_length);
                p__peakSeqLiveMemory = max(p__peakSeqLiveMemory, g__SEQNODE_ALLOCATOR->GetSeqLiveMemory());

#ifdef VELOUR_TBB
                tbb::tick_count time1 = tbb::tick_count::now();
#endif

#ifdef VELOUR_TBB
                uintptr_t nodes_memoized = sg_parallel_chain_memoize_related_components(inbox_chain, resident_sgraph, &flowlist);
                //uintptr_t nodes_memoized = sg_serial_chain_memoize_related_components(inbox_chain, resident_sgraph, &flowlist);
#else
                uintptr_t nodes_memoized = sg_serial_chain_memoize_related_components(inbox_chain, resident_sgraph, &flowlist);
#endif

                round_import_nodes += nodes_memoized;

#ifdef VELOUR_TBB
                tbb::tick_count time2 = tbb::tick_count::now();
#endif

                while (inbox_chain != NULL) {
                    SeqNode *node = inbox_chain;
                    inbox_chain = inbox_chain->head_next;
                    node->head_next = NULL; // XXX
                    flowlist.push_back(node);
                    resident_sgraph->insertNodeAndUpdateColors(node);
                }

#ifdef VELOUR_TBB
                tbb::tick_count time3 = tbb::tick_count::now();

                time_inboxload += time1 - time0;
                time_relatedcomponents += time2 - time1;
                time_insertnodes += time3 - time2;
#endif

                //printf("(flow) inbox component nodes: %"PRIuPTR"\n", round_inbox_nodes); // XXX: remove me
                //printf("(flow) inbox component related nodes: %"PRIuPTR"\n", round_import_nodes); fflush(stdout); // XXX: remove me

#ifdef VERIFY
                resident_sgraph->verify(true);
#endif

#ifdef VELOUR_TBB
                tbb::tick_count time4 = tbb::tick_count::now();
#endif

#ifdef VELOUR_TBB
                //TODO sg_parallel_nodelist_remove_tips(resident_sgraph, true, &flowlist);
                sg_nodelist_remove_tips(resident_sgraph, true, &flowlist);
#else
                sg_nodelist_remove_tips(resident_sgraph, true, &flowlist);
#endif

#ifdef VELOUR_TBB
                tbb::tick_count time5 = tbb::tick_count::now();
#endif

#ifdef VERIFY
                resident_sgraph->verify(true);
#endif

                //sg_concatenate(resident_sgraph, true);
                sg_nodelist_concatenate(resident_sgraph, true, &flowlist);
#ifdef VELOUR_TBB
                tbb::tick_count time6 = tbb::tick_count::now();
#endif

#ifdef VERIFY
                resident_sgraph->verify(true);
#endif

                // FIXME: execution time is counted with concatenation
                //sg_nodelist_pop_bubbles(resident_sgraph, true, &flowlist);
#ifdef VELOUR_TBB
                //tbb::tick_count time6 = tbb::tick_count::now();
#endif

#ifdef VERIFY
                //resident_sgraph->verify(true);
#endif


                // emit sub-components that are no longer relevant
                //slice2_graph(resident_sgraph, g__PARTITION_INDEX);
                if (g__SLICING) slice2_nodelist(resident_sgraph, g__PARTITION_INDEX, &flowlist);

#ifdef VELOUR_TBB
                tbb::tick_count time7 = tbb::tick_count::now();
#endif

                // emit components that are no longer relevant to the working graph
                //buckets->split(resident_sgraph);
                buckets->split_nodelist(resident_sgraph, &flowlist);

#ifdef VELOUR_TBB
                tbb::tick_count time8 = tbb::tick_count::now();
#endif

                // lastly, reset claim tid on nodes for future iterations
                for (flow_nodelist_t::iterator it = flowlist.begin(); it != flowlist.end() ; ++it) {
                    SeqNode *node = *it;
                    //assert( node->claim_tid == 1 ); // don't assert this as not true for grown nodes and for inbox nodes
                    node->claim_tid = 0;
                }

                flowlist.clear();

#ifdef VELOUR_TBB
                tbb::tick_count time9 = tbb::tick_count::now();

                time_removetips += time5 - time4;
                time_concatenate += time6 - time5;
                time_slicing += time7 - time6;
                time_splitting += time8 - time7;
                time_cleanup += time9 - time8;
#endif

                // initiate garbage collection if needed
                if (g__NODE_ALLOCATORS->flag_gc_needed_) { g__NODE_ALLOCATORS->GarbageCollect(); }

                ++ round;
            }

            if (munmap(file_mmap, file_length) == -1) {
                fprintf(stderr, "ERROR: failed to munmap file: %s\n", itr->filename);
                perror("REASON: ");
                exit(EXIT_FAILURE);
            }

            if (close(filedes) == -1) {
                fprintf(stderr, "ERROR: failed to close file: %s\n", itr->filename);
                perror("REASON: ");
                exit(EXIT_FAILURE);
            }

#ifdef VELOUR_MPI
checkdone:
            MPI_Status mpi_status;
            MPI_Iprobe(PART_TO_RANK(inbox_partition), 42, MPI_COMM_WORLD, &producer_finished, &mpi_status);
            // TODO: do MPI_Recv() on the message()
            } while (!producer_finished || safe_length != * static_cast<volatile size_t *>(&p__MPI_SAFE_LENGTH_ARRAY[inbox_partition]));
#endif // VELOUR_MPI
        }

        /*if( g__FULL_STATISTICS ) {
          sg_stat_components(resident_sgraph, stdout);
          }*/

        // last, split the resident graph to outboxes // FIXME: should this be empty!?
        printf("FLOW: %zu nodes (FIXME maybe dead) in resident graph at end of flowing.... slicing and splitting resident graph.\n",
                resident_sgraph->node_count); fflush(stdout); // XXX: HACK
        if (g__SLICING) { slice2_graph(resident_sgraph, g__PARTITION_INDEX); } // XXX: HACK HACK HACK
        buckets->split(resident_sgraph);
        //assert( resident_sgraph->node_count == 0 );
    }

#ifdef VELOUR_TBB
    tbb::tick_count::interval_t time_total = time_inboxload + time_relatedcomponents + time_insertnodes +
        time_removetips + time_concatenate + time_slicing + time_splitting + time_cleanup;

    printf("  flow time: %lfs  inbox: %02.1lf%%  related: %02.1lf%%  insert: %02.1lf%%"
           "  remove: %02.1lf%%  concat: %02.1lf%%  slice: %02.1lf%%  split: %02.1lf%%"
           "  clean: %02.1lf%%\n",
            time_total.seconds(), 
            100.0 * (time_inboxload.seconds() / time_total.seconds()),
            100.0 * (time_relatedcomponents.seconds() / time_total.seconds()),
            100.0 * (time_insertnodes.seconds() / time_total.seconds()),
            100.0 * (time_removetips.seconds() / time_total.seconds()),
            100.0 * (time_concatenate.seconds() / time_total.seconds()),
            100.0 * (time_slicing.seconds() / time_total.seconds()),
            100.0 * (time_splitting.seconds() / time_total.seconds()),
            100.0 * (time_cleanup.seconds() / time_total.seconds()) );
    fflush(stdout);
#endif
}

void runFlow(KmerGraph *resident_kgraph, SeqGraph *resident_sgraph)
{
    char loom_filename[PATH_MAX+1];
    file_object_t next_file;

    sprintf(loom_filename, "%s/Subsequences-%u.loom", g__WORK_LOOM_DIRECTORY, g__PARTITION_INDEX);
    next_file.filetype = LOOM;
    next_file.filename = loom_filename;

    file_object_vector loom_file_object;
    loom_file_object.push_back(next_file);

#ifdef VELOUR_TBB
    tbb::tick_count::interval_t loom_time;
    tbb::tick_count::interval_t kmerseq_time;
    tbb::tick_count::interval_t inbox_time;

    tbb::tick_count time0, time1;
    time0 = tbb::tick_count::now();
#endif // VELOUR_TBB

    SplitBuckets *outbox_buckets = new SplitBuckets(true);

#ifdef VELOUR_MPI
    MPI_Barrier(MPI_COMM_WORLD); // NOTE: barrier to ensure the SplitBuckets are all created (above)

    MPI_Alloc_mem(g__PARTITION_COUNT * sizeof(size_t), MPI_INFO_NULL, &p__MPI_SAFE_LENGTH_ARRAY);
    MPI_Win_create(p__MPI_SAFE_LENGTH_ARRAY, g__PARTITION_COUNT*sizeof(size_t), sizeof(size_t), MPI_INFO_NULL, MPI_COMM_WORLD, &g__MPI_WINDOW_SAFE_LENGTHS);

    for (unsigned u=0; u < g__PARTITION_COUNT; ++u) {
        p__MPI_SAFE_LENGTH_ARRAY[u] = 0;
    }

    MPI_Barrier(MPI_COMM_WORLD);
#endif // VELOUR_MPI

    // next, load the subsequences into kmer graph
    load_loom_files(loom_file_object);

#ifdef VELOUR_TBB
    time1 = tbb::tick_count::now();
    loom_time = time1 - time0;
    printf("  loom time: %lfs\n", loom_time.seconds()); fflush(stdout);
#endif // VELOUR_TBB

    if (g__DOTGRAPH) {
      char dot_filename[PATH_MAX+1];
      sprintf(dot_filename, "%s/Flow-initial-kmergraph-%u.dot", g__WORK_BASE_DIRECTORY, g__PARTITION_INDEX);
      emit_graphviz(resident_kgraph, dot_filename);
    }

    printf("%" PRIuPTR " kmer nodes built from Loom file(s)\n", resident_kgraph->node_count); fflush(stdout);

    p__peakKmerNodes = resident_kgraph->node_count;
    p__peakKmerLiveMemory = g__KMERNODE_ALLOCATOR->GetKmerLiveMemory();

    if( g__FULL_STATISTICS ) {
        kg_stat_components(resident_kgraph, stdout);
    }

    /*remove_tips(resident_kgraph);
    if (g__DOTGRAPH) {
      char dot_filename[PATH_MAX+1];
      sprintf(dot_filename, "%s/Flow-clipped-kmergraph-%u.dot", g__WORK_BASE_DIRECTORY, g__PARTITION_INDEX);
      emit_graphviz(resident_kgraph, dot_filename);
    }*/

    //sg_dump_pregraph_from_kmergraph(resident_kgraph, "ResidentGraph-after-tipclip.pregraph");

    // XXX: allocate sequence graph here instead?


#ifdef VELOUR_TBB
    tbb::tick_count time2, time3;
    time2 = tbb::tick_count::now();
#endif // VELOUR_TBB

    lambda_flowKmerComponent f(resident_kgraph, outbox_buckets, resident_sgraph);
#ifdef VELOUR_TBB
    kg_for_each_parallel(resident_kgraph, f); // TODO: GC version of this iterator
#else
    kg_for_each(resident_kgraph, f, true); // TODO: enumerate the true constant for GC
#endif

#ifdef VELOUR_TBB
    time3 = tbb::tick_count::now();
    kmerseq_time = time3 - time2;
    printf("  kmerseq time: %lfs\n", kmerseq_time.seconds()); fflush(stdout);
#endif // VELOUR_TBB

    if (g__DOTGRAPH) {
      char dot_filename[PATH_MAX+1];
      sprintf(dot_filename, "%s/Flow-initial-seqgraph-%u.dot", g__WORK_BASE_DIRECTORY, g__PARTITION_INDEX);
      emit_graphviz(resident_sgraph, dot_filename);
    }

    /*char str[300];
    sprintf(str, "ResidentGraph-initial-seqgraph-%u.pregraph", g__PARTITION_INDEX);
    sg_dump_pregraph(resident_sgraph, str);*/

    // TODO: this block of code should be a function
    g__KMERNODE_ALLOCATOR->BulkFreeListAllSlabs();
    g__NODE_ALLOCATORS->flag_kmer_gc_needed_ = false;
    delete resident_kgraph;
    g__KMERNODE_ALLOCATOR->set_kmergraph(NULL);
    // END BLOCK

    printf("FLOW: kmerGraph -> seqGraph phase finished.\n");  fflush(stdout);

    // to compute initial reduction
    printf("%"PRIuPTR" total sequence nodes pre-redistribution (kmer flow).\n", p__preRedistributionSequenceNodes);
    //printf("%zu MB total sequence node memory pre-redistribution (kmer flow).\n", p__preRedistributionSequenceNodeMemory);

    printf("%"PRIuPTR" initial sequence nodes in resident graph (flow).\n", resident_sgraph->node_count);

    p__peakSeqNodes = resident_sgraph->node_count; // initializes

    uintptr_t total_bucket_nodes = 0;

    p__peakSeqLiveMemory = g__SEQNODE_ALLOCATOR->GetSeqLiveMemory(); // initializes -- but XXX: doesn't count flowKmerComponent stuff

#ifdef VELOUR_TBB
    tbb::tick_count time4, time5;
    time4 = tbb::tick_count::now();
#endif // VELOUR_TBB

#ifdef VERIFY
    resident_sgraph->verify(false);
#endif

    file_object_vector inbox_file_objects;

    next_file.filetype = BUCKET;
    for (unsigned i=1; i < g__PARTITION_INDEX; ++i) {
        char *filename = (char *)malloc((PATH_MAX+1) * sizeof(char));
        sprintf(filename, "%s/%u/InboxBucket-from-%u.bucket", g__WORK_INBOX_ROOT_DIRECTORY, g__PARTITION_INDEX, i);
        next_file.filename = filename;
        next_file.fileindex = i; // MPI: record the inbox source partition index
        next_file.length = 0; // MPI: memoize progress as we rotate through the inboxes
        inbox_file_objects.push_back(next_file);
    }

#ifdef VELOUR_TBB
    newflow_inbox(inbox_file_objects, resident_sgraph, outbox_buckets, total_bucket_nodes); // XXX
    //parallel_process_inbox(inbox_file_objects, resident_sgraph, outbox_buckets, peak_nodes, total_bucket_nodes); 
#else
    newflow_inbox(inbox_file_objects, resident_sgraph, outbox_buckets, total_bucket_nodes); // XXX
    //serial_process_inbox(inbox_file_objects, resident_sgraph, outbox_buckets, peak_nodes, total_bucket_nodes);
#endif

#ifdef VELOUR_TBB
    time5 = tbb::tick_count::now();
    inbox_time = time5 - time4;
    printf("  inbox time: %lfs\n", inbox_time.seconds());
    printf("  total time: %lfs  loom: %02.1lf%%  kmerseq: %02.1lf%%  inbox: %02.1lf%%\n",
            loom_time.seconds() + kmerseq_time.seconds() + inbox_time.seconds(),
            100.0 * (loom_time.seconds()/(loom_time.seconds() + kmerseq_time.seconds() + inbox_time.seconds())),
            100.0 * (kmerseq_time.seconds()/(loom_time.seconds() + kmerseq_time.seconds() + inbox_time.seconds())),
            100.0 * (inbox_time.seconds()/(loom_time.seconds() + kmerseq_time.seconds() + inbox_time.seconds())) );
    fflush(stdout);
#endif // VELOUR_TBB

    /*if (g__DOTGRAPH) {
      char dot_filename[PATH_MAX+1];
      sprintf(dot_filename, "%s/Flow-postinboxes-seqgraph-%u.dot", g__WORK_BASE_DIRECTORY, g__PARTITION_INDEX);
      emit_graphviz(resident_sgraph, dot_filename);
    }*/

    printf("%"PRIuPTR" peak kmer nodes (flow).\n", p__peakKmerNodes);
    printf("%"PRIuPTR" peak sequence nodes (flow).\n", p__peakSeqNodes);
    printf("%zu MB peak kmer node live memory (flow).\n", p__peakKmerLiveMemory / (1024 * 1024));
    printf("%zu MB peak sequence node live memory (flow).\n", p__peakSeqLiveMemory / (1024 * 1024));
    printf("%zu MB peak live memory (flow).\n", max(g__peakLiveMemory, max(p__peakKmerLiveMemory, p__peakSeqLiveMemory)) / (1024 * 1024));
    printf("%"PRIuPTR" total sequence nodes loaded from stream bucket(s).\n", total_bucket_nodes);

    printf("%"PRIuPTR" nodes sliced out to final.\n", g__SLICE2_FINAL_NODE_COUNT);
    printf("%"PRIuPTR" nodes sliced out to others.\n", g__SLICE2_NODE_COUNT);

    outbox_buckets->printStatistics();
    delete outbox_buckets;

    //delete inbox_file_objects; and free the filenames

#ifdef VELOUR_MPI
    // signal completion to other partitions
    unsigned consumer_count = 0;
    unsigned mpi_message = 0; // NOTE: same message used in all Isend below
    MPI_Request mpi_requests[g__PARTITION_COUNT+1];
    MPI_Status  mpi_status[g__PARTITION_COUNT+1];

    for (unsigned u=g__PARTITION_INDEX+1; u <= g__PARTITION_COUNT; ++u) {
        int ret = MPI_Isend(&mpi_message, 1, MPI_UNSIGNED, PART_TO_RANK(u), 42, MPI_COMM_WORLD, &mpi_requests[consumer_count]);
        assert( ret == MPI_SUCCESS );
        ++consumer_count;
    }
    int ret = MPI_Waitall(consumer_count, mpi_requests, mpi_status);
    assert( ret == MPI_SUCCESS );

    // TODO: check the mpi_status?

    MPI_Barrier(MPI_COMM_WORLD); // needed?

    MPI_Win_free(&g__MPI_WINDOW_SAFE_LENGTHS);
    MPI_Free_mem(p__MPI_SAFE_LENGTH_ARRAY);
#endif // VELOUR_MPI

}
