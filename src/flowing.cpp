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

static size_t p__peakFragmentNodes = 0;
static size_t p__peakFragmentMemory = 0;

static uintptr_t p__preRedistributionSequenceNodes = 0;
//static size_t p__preRedistributionSequenceNodeMemory = 0;

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

namespace {
struct lambda_serial_import_components {
    lambda_serial_import_components(SeqGraph *working_graph, SeqGraph *resident_graph) :
        working_graph(working_graph), resident_graph(resident_graph) {}
    void operator()(SeqNode *node) {
        for (int i = 0 ; i < 4 ; ++ i) {
            // check on the left side
            if (node->left_color[i] == g__PARTITION_INDEX) {
                assert( node->left_count[i] != 0 );
                SeqNode *next = resident_graph->findNextNode(node, i, GO_LEFT);
                assert( next != NULL );
                SerialComponent<SeqGraph, SeqNode> component(resident_graph, next);
                if (!component.nodes_.empty()) { // empty if I already claimed it with another edge
                    for (std::deque<SeqNode*>::iterator it=component.nodes_.begin(); it != component.nodes_.end(); ++it) {
                        SeqNode *node = *it;
                        resident_graph->removeNode(node);
                        working_graph->insertNodeAndUpdateColors(node);
                    }
                }
            }
            // check on the right side
            if (node->right_color[i] == g__PARTITION_INDEX) {
                assert( node->right_count[i] != 0 );
                SeqNode *next = resident_graph->findNextNode(node, i, GO_RIGHT);
                assert( next != NULL );
                SerialComponent<SeqGraph, SeqNode> component(resident_graph, next);
                if (!component.nodes_.empty()) { // empty if I already claimed it with another edge
                    for (std::deque<SeqNode*>::iterator it=component.nodes_.begin(); it != component.nodes_.end(); ++it) {
                        SeqNode *node = *it;
                        resident_graph->removeNode(node);
                        working_graph->insertNodeAndUpdateColors(node);
                    }
                }
            }
        }
    }
  private:
    SeqGraph *working_graph;
    SeqGraph *resident_graph;
};
} // namespace: anonymous

void sg_serial_import_related_components(SeqGraph *working_graph, SeqGraph *resident_graph)
{
    lambda_serial_import_components f(working_graph, resident_graph);
    sg_for_each(working_graph, f);
}

uintptr_t sg_serial_chain_import_related_components(SeqNode *seed_chain, SeqGraph *source_graph, SeqNode **dest_chain)
{
    uintptr_t nodes_imported = 0;

    SeqNode *iterator = seed_chain;
    while (iterator != NULL) {
        SeqNode *root = iterator;
        for (int i = 0 ; i < 4 ; ++ i) {
            // check on the left side
            if (root->left_color[i] == g__PARTITION_INDEX) {
                assert( root->left_count[i] != 0 );
                SeqNode *next = source_graph->findNextNode(root, i, GO_LEFT);
                if (next != NULL) {
                    SerialComponent<SeqGraph, SeqNode> component(source_graph, next);
                    if (!component.nodes_.empty()) { // empty if I already claimed it with another edge
                        for (std::deque<SeqNode*>::iterator it=component.nodes_.begin(); it != component.nodes_.end(); ++it) {
                            SeqNode *node = *it;
                            source_graph->removeNode(node);

                            // chain it
                            node->head_next = *dest_chain;
                            *dest_chain = node;
                            ++ nodes_imported;
                        }
                    }
                } else {
                    // XXX: VERIFY that NULL target is in dest_chain -- as we expected to find it
                }
            }
            // check on the right side
            if (root->right_color[i] == g__PARTITION_INDEX) {
                assert( root->right_count[i] != 0 );
                SeqNode *next = source_graph->findNextNode(root, i, GO_RIGHT);
                if (next != NULL) {
                    SerialComponent<SeqGraph, SeqNode> component(source_graph, next);
                    if (!component.nodes_.empty()) { // empty if I already claimed it with another edge
                        for (std::deque<SeqNode*>::iterator it=component.nodes_.begin(); it != component.nodes_.end(); ++it) {
                            SeqNode *node = *it;
                            source_graph->removeNode(node);

                            // chain it
                            node->head_next = *dest_chain;
                            *dest_chain = node;
                            ++ nodes_imported;
                        }
                    }
                } else {
                    // XXX: VERIFY that NULL target is in dest_chain -- as we expected to find it
                }
            }
        }
        iterator = iterator->head_next;
    }

    return nodes_imported;
}

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

#ifdef VELOUR_TBB
namespace {
struct boundary_edge {
    SeqNode    *node;
    Nucleotide  nucleotide;
    unsigned    direction;
};

typedef std::deque<boundary_edge> edgelist_t;

typedef std::deque<SeqNode *> nodelist_t;

void gdb_list_nodelist_seq(nodelist_t &nodelist)
{
    for (nodelist_t::iterator nit = nodelist.begin(); nit != nodelist.end(); ++nit) {
        SeqNode *node = *nit;
        printf("node: %p  claim_tid: %u\n", node, node->claim_tid);
    }
    fflush(stdout);
}

struct lambda_enumerateAdjacentBoundaryEdges {
    lambda_enumerateAdjacentBoundaryEdges(edgelist_t &el) : edgelist(el) {}
    void operator()(SeqNode *node) {
        for (int i = 0 ; i < 4 ; ++ i) {
            if (node->left_color[i] == g__PARTITION_INDEX) {
                assert( node->left_count[i] != 0 );
                boundary_edge b;
                b.node = node;
                b.nucleotide = i;
                b.direction = GO_LEFT;
                edgelist.push_back(b);
            }
            if (node->right_color[i] == g__PARTITION_INDEX) {
                assert( node->right_count[i] != 0 );
                boundary_edge b;
                b.node = node;
                b.nucleotide = i;
                b.direction = GO_RIGHT;
                edgelist.push_back(b);
            }
        }
    }
    private:
      edgelist_t &edgelist;
};
} // namespace: anonymous

static tbb::atomic<uintptr_t> component_counter; // hack to detect when NULL nodes might be available again

static void unclaim_nodes(nodelist_t &nodelist)
{
    for (nodelist_t::iterator nit = nodelist.begin(); nit != nodelist.end(); ++nit) {
        SeqNode *node = *nit;
        //atomic_resetNodeClaimed(node);
        ATOMIZE(node->claim_tid).compare_and_swap(ParallelComponent<SeqGraph,SeqNode>::UNCLAIMED, tls_thread_index);
    }
    nodelist.clear();
}

static bool atomic_claim_components(SeqGraph *resident_graph, SeqGraph *inbox_graph, edgelist_t &edgelist, nodelist_t &nodelist)
{
    uintptr_t old_component_counter = component_counter;

    assert( nodelist.empty() );
    for (edgelist_t::iterator eit = edgelist.begin() ; eit != edgelist.end() ; ++eit) {

        boundary_edge edge = *eit;

        bool claimed_component = false;
        while(!claimed_component) {
            SeqNode *next = resident_graph->atomic_findNextNode(edge.node, edge.nucleotide, edge.direction);

            if (next != NULL) {
                if (ATOMIC_LOAD(next->claim_tid) == tls_thread_index) { // already claimed by me
                    claimed_component = true;
                    continue;
                }
                ParallelComponent<SeqGraph, SeqNode> component(resident_graph, next);
                if (!component.nodes_.empty()) { // successfully claimed component
                    for (std::deque<SeqNode*>::iterator it=component.nodes_.begin(); it != component.nodes_.end(); ++it) {
                        SeqNode *node = *it;
                        //atomic_setNodeClaimed(node);
                        nodelist.push_back(node);
                    }
                    claimed_component = true;
                    continue;
                } else { // didn't claim component
                    if (ATOMIC_LOAD(*tls_thread_abort)) { // check if I should abort
                        unclaim_nodes(nodelist);
                        ATOMIC_STORE(*tls_thread_abort) = false; // reset abort flag
                        return false;
                    }
                    threadid_t other_thread = ATOMIC_LOAD(next->claim_tid);
                    if (other_thread > 0 && other_thread < tls_thread_index) { // abort other thread -- okay if false abort if claim changes
                        if (ATOMIC_LOAD(*thread_aborts[other_thread]) == false) { // check if false to avoid cache-line ping-pong
                            ATOMIC_STORE(*thread_aborts[other_thread]) = true;
                        }
                    }
                }
            } else if (inbox_graph->findNextNode(edge.node, edge.nucleotide, edge.direction) != NULL) { // already have it
                claimed_component = true;
                continue;
            } else { // next == NULL -- unclaim all, wait for someone to finish, then restart
                unclaim_nodes(nodelist);
                while (component_counter == old_component_counter) { tbb::this_tbb_thread::yield(); }
                if (ATOMIC_LOAD(*tls_thread_abort)) {
                    ATOMIC_STORE(*tls_thread_abort) = false; // reset abort flag
                }
                return false;
            }
            if (ATOMIC_LOAD(*tls_thread_abort)) { // check if I should abort
                    unclaim_nodes(nodelist);
                    ATOMIC_STORE(*tls_thread_abort) = false; // reset abort flag
                    return false;
            }
        }
    }
    return true;
}


uintptr_t sg_atomic_chain_import_related_components(SeqGraph *inbox_graph, SeqGraph *resident_graph, SeqNode **dest_chain)
{
    uintptr_t node_counter = 0;

    edgelist_t edgelist;
    lambda_enumerateAdjacentBoundaryEdges e(edgelist);
    sg_for_each(inbox_graph, e);

    nodelist_t nodelist;  // nodes of all components that have been claimed

    ATOMIC_STORE(*tls_thread_abort) = false; // reset abort flag to avoid unrelated aborts

    while(!atomic_claim_components(resident_graph, inbox_graph, edgelist, nodelist)) ; // retry until all components claimed

    // remove nodes from the resident graph and insert in private graph
    for (nodelist_t::iterator nit = nodelist.begin(); nit != nodelist.end(); ++nit) {
        SeqNode *node = *nit;
        assert( node->claim_tid == tls_thread_index );
        resident_graph->atomic_removeNode(node);

        node->head_next = *dest_chain;
        *dest_chain = node;

        ++ node_counter;
    }
    return node_counter;
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

    p__preRedistributionSequenceNodes += private_sgraph->node_count;

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

#ifdef VELOUR_TBB
typedef struct inbox_token {
    char * component_start;
    size_t num_components;
} inbox_token_t;

class InboxLoadFilter : public tbb::filter
{
    public:
        static const size_t n_buffer = 64;  // XXX: controls the concurrency level

    private:
        file_object_t **file_object_;

        int filedes_;
        off_t file_offset_;
        char * file_mmap_;
        size_t file_length_;

        uintptr_t touch_total_;  // fool the optimizer

        inbox_token_t buffer_[n_buffer];
        size_t next_buffer_;

        bool open_file()
        {
            /*if (filedes_ != -1) { // close current file
                //if (file_length_ > 0) { FIXME FIXME FIXME can't munmap here
                //    if (munmap(file_mmap_, file_length_) == -1) {
                //        fprintf(stderr, "ERROR: failed to munmap file: %s\n", files_to_load_[file_index_].filename);
                //        perror("REASON: ");
                //        exit(EXIT_FAILURE);
                //    }
                //}
                if (close(filedes_) == -1) {
                    fprintf(stderr, "ERROR: failed to close file: %s\n", files_object_.filename);
                    perror("REASON: ");
                    exit(EXIT_FAILURE);
                }
                filedes_ = -1;
            }*/

            file_offset_ = 0;   // reset offset

            if ((*file_object_)->filetype != BUCKET) {
                fprintf(stderr, "ERROR: Cannot use %s file %s as flow input.  Exiting...\n",
                    FILE_TYPES[(*file_object_)->filetype], (*file_object_)->filename);
                exit(EXIT_FAILURE);
            }

            filedes_ = open((*file_object_)->filename, O_RDONLY);
            if (filedes_ == -1) {
                fprintf(stderr, "ERROR: failed to open file: %s\n", (*file_object_)->filename);
                perror("REASON: ");
                exit(EXIT_FAILURE);
            }

            struct stat file_stat;
            if (fstat(filedes_, &file_stat) != 0) {
                fprintf(stderr,"ERROR: failed to stat file: %s\n", (*file_object_)->filename);
                exit(EXIT_FAILURE);
            }

            file_length_ = file_stat.st_size;
            (*file_object_)->length = file_stat.st_size;
            
            if (file_length_ == 0) { // otherwise, mmap will fail
                return false;
            }

            file_mmap_ = static_cast<char *>( mmap( 0, file_length_, PROT_READ, MAP_PRIVATE, filedes_, 0) );
            if (file_mmap_ == reinterpret_cast<char *>(-1)) {
                fprintf(stderr, "ERROR: failed to mmap file: %s\n", (*file_object_)->filename);
                perror("REASON: ");
                exit(EXIT_FAILURE);
            }

            // TODO: walk the pages???

            return true;
        }

    public:
        InboxLoadFilter(file_object_t **inbox_file) : tbb::filter(serial_in_order),
            file_object_(inbox_file), filedes_(-1), file_offset_(0),
            file_mmap_(NULL), file_length_(0), touch_total_(0), next_buffer_(0)
        { }

        void* operator()(void* __dummy)
        {
            //printf("inbox-load %u\n", tls_thread_index); fflush(stdout);
            if (file_offset_ >= file_length_) { // finished file or just starting file
                if ((*file_object_) != NULL) {
                    if (!open_file()) {
                        return NULL;
                    }
                    (*file_object_) = NULL;
                } else {
                    // TODO TODO TODO : close and munmap file sometime
                    return NULL;
                }
            }

            inbox_token_t *buf = &buffer_[next_buffer_];
            next_buffer_ = (next_buffer_ + 1) % n_buffer;

            buf->component_start = &file_mmap_[file_offset_];
            buf->num_components = 0;

            // grab components and walk pages
            off_t length_to_grab = min( 4*4*16384, (file_length_ - file_offset_));
            assert( file_offset_ + length_to_grab <= file_length_ );
                
            //printf("buf %p length_to_grab = %zd\n", buf, length_to_grab);

            while (length_to_grab > 0) {
                size_t serialized_bytes;
                memcpy(&serialized_bytes, (file_mmap_ + file_offset_), sizeof(size_t));
                assert( serialized_bytes <= file_length_ - file_offset_ );

                //printf("buf %p serialized_bytes = %zu\n", buf, serialized_bytes);

                //printf("touching %u\n", tls_thread_index); fflush(stdout);
                // touch the pages for the component
                for (size_t touch = 0 ; touch < serialized_bytes ; touch += 4096) { // XXX: constant
                    touch_total_ += file_mmap_[file_offset_ + touch];
                }

                ++ (buf->num_components);
                file_offset_ += serialized_bytes;
                length_to_grab -= serialized_bytes;
            }
            assert( buf->num_components > 0 );

            //printf("buf %p component_start = %p\n", buf, buf->component_start);
            //printf("buf %p  num_components = %zu\n", buf, buf->num_components);
            //printf("buf %p next_component  = %p\n", buf, file_mmap_ + file_offset_);
            //fflush(stdout);
            return buf;
        }
};

class InboxProcessFilter : public tbb::filter
{
    public:
        InboxProcessFilter(SeqGraph *resident_sgraph, SplitBuckets *buckets,
            tbb::atomic<uintptr_t> &peak_nodes,tbb::atomic<uintptr_t> &total_bucket_nodes) :
            tbb::filter(parallel), resident_sgraph(resident_sgraph), buckets(buckets),
            peak_nodes(peak_nodes), total_bucket_nodes(total_bucket_nodes) {}

        void * operator()(void *item)
        {
            inbox_token_t& token = * static_cast<inbox_token_t *>(item);

            char * next_component_start = token.component_start;

            size_t components_left = token.num_components;

            while (components_left > 0) {
                SeqNode *inbox_chain = NULL;
                SeqNode *import_chain = NULL;

                uintptr_t inbox_components = 0;
                uintptr_t inbox_nodes = 0;
                while (components_left > 0 && inbox_nodes < (100000 / 2)) { // XXX: configurable constant
                    size_t component_size;
                    uintptr_t inbox_chain_length = sg_load_mmap_stream_component(next_component_start, &component_size, &inbox_chain);
                    next_component_start += component_size;
                    ++ inbox_components;
                    inbox_nodes += inbox_chain_length;

                    // NOTE: can't atomic import components here

                    -- components_left;
                }

                // TODO: peak nodes
                //uintptr_t new_peak = max(old_peak, resident_sgraph->node_count + private_sgraph->node_count);

                SeqGraph *inbox_sgraph = tls_private_sgraph_set->GetPrivateGraph(inbox_nodes);
                assert( inbox_sgraph->node_count == 0 );

                while (inbox_chain != NULL) {
                    SeqNode *node = inbox_chain;
                    inbox_chain = inbox_chain->head_next;
                    node->head_next = NULL; // XXX: needed?
                    inbox_sgraph->insertNodeAndUpdateColors(node); // XXX: updating colors too since multiple components!
                }
                assert( inbox_chain == NULL );

                uintptr_t import_chain_length = sg_atomic_chain_import_related_components(inbox_sgraph, resident_sgraph, &import_chain);

                uintptr_t total_private_nodes = inbox_nodes + import_chain_length;
                SeqGraph *private_sgraph = tls_private_sgraph_set->GetPrivateGraph(total_private_nodes);

                if (private_sgraph != inbox_sgraph) {
                    assert( private_sgraph->node_count == 0 );
                    inbox_sgraph->bulkMoveAllNodes(private_sgraph);
                }

                // then, add imported nodes and fixup colors
                while (import_chain != NULL) {
                    SeqNode *node = import_chain;
                    import_chain = import_chain->head_next;
                    node->head_next = NULL; // XXX: needed?
                    private_sgraph->insertNodeAndUpdateColors(node);
                }

                uintptr_t round = (token.num_components - components_left) + 1;
                flow_private_graph(private_sgraph, resident_sgraph, buckets, round);
                assert( private_sgraph->node_count == 0 );

                component_counter.fetch_and_add(inbox_components);
            }

            return NULL;
        }
    private:
        SeqGraph *resident_sgraph;
        SplitBuckets *buckets;
        tbb::atomic<uintptr_t> &peak_nodes;
        tbb::atomic<uintptr_t> &total_bucket_nodes;
};

void parallel_process_inbox(file_object_vector &inbox_file_objects, SeqGraph *resident_sgraph, SplitBuckets *buckets,
        uintptr_t &peak_nodes, uintptr_t &total_bucket_nodes)
{
    tbb::tick_count time0, time1;
    tbb::pipeline pipeline;

    file_object_t *current_inbox = NULL;

    printf("  parallel inbox: %zu tokens\n", InboxLoadFilter::n_buffer); fflush(stdout);

    time0 = tbb::tick_count::now();

    InboxLoadFilter lf(&current_inbox);
    pipeline.add_filter( lf );

    tbb::atomic<uintptr_t> atomic_peak_nodes;
    atomic_peak_nodes = peak_nodes;
    tbb::atomic<uintptr_t> atomic_total_bucket_nodes;
    atomic_total_bucket_nodes = total_bucket_nodes;

    InboxProcessFilter pf(resident_sgraph, buckets, atomic_peak_nodes, atomic_total_bucket_nodes);
    pipeline.add_filter( pf );

    // run the pipeline for each inbox file, since multiple inboxes can't in general be run in parallel
    for (file_object_vector::iterator it=inbox_file_objects.begin(); it != inbox_file_objects.end(); ++it) {
        current_inbox = &(*it);
        printf("    flowing inbox: %s\n", it->filename); fflush(stdout);
        pipeline.run( InboxLoadFilter::n_buffer );
    }

    time1 = tbb::tick_count::now();
    
    peak_nodes = atomic_peak_nodes;
    total_bucket_nodes = atomic_total_bucket_nodes;

    tbb::tick_count::interval_t loom_time = time1 - time0;

    printf("  inbox time: %lfs\n", loom_time.seconds()); fflush(stdout);
}
#endif // VELOUR_TBB

} // namespace: anonymous
    
void serial_process_inbox(file_object_vector &inbox_file_objects, SeqGraph *resident_sgraph, SplitBuckets *buckets,
        uintptr_t &peak_nodes, uintptr_t &total_bucket_nodes)
{
    if (!inbox_file_objects.empty()) { // no inbox for the first partition

        PrivateSeqGraphSet psg_set;

        // incrementally load components from the inbox files
        for (file_object_vector::iterator itr=inbox_file_objects.begin() ; itr != inbox_file_objects.end() ; ++itr) {
            if (itr->filetype != BUCKET) {
                fprintf(stderr, "ERROR: Cannot use %s file %s as flow input.  Exiting...\n",
                    FILE_TYPES[itr->filetype], itr->filename);
                exit(EXIT_FAILURE);
            }

            printf("    flowing inbox: %s\n", itr->filename); fflush(stdout);
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

            if (file_length == 0) continue;  // otherwise, mmap will fail

            char *file_mmap = static_cast<char *>( mmap( 0, file_length, PROT_READ, MAP_PRIVATE, filedes, 0) );
            if (file_mmap == reinterpret_cast<char *>(-1)) {
                fprintf(stderr, "ERROR: failed to mmap file: %s\n", itr->filename);
                perror("REASON: ");
                exit(EXIT_FAILURE);
            }

            size_t file_offset = 0;

            uintptr_t round = 1;
            while (file_offset < file_length) {
                SeqNode *bulk_inbox_chain = NULL;
                SeqNode *import_chain = NULL;

                uintptr_t round_inbox_nodes = 0;
                uintptr_t round_import_nodes = 0;
                uintptr_t total_private_nodes = 0;
                while (file_offset < file_length && total_private_nodes < 100000) { // XXX: configurable parameter

                    SeqNode *inbox_chain = NULL;
                    size_t component_size;
                    uintptr_t inbox_chain_length = sg_load_mmap_stream_component((file_mmap+file_offset), &component_size, &inbox_chain);
                    file_offset += component_size;
                    total_bucket_nodes += inbox_chain_length;
                    round_inbox_nodes += inbox_chain_length;

                    peak_nodes = max(peak_nodes, resident_sgraph->node_count + inbox_chain_length);

                    uintptr_t import_chain_length = sg_serial_chain_import_related_components(inbox_chain, resident_sgraph, &import_chain);

                    total_private_nodes += inbox_chain_length + import_chain_length;
                    round_import_nodes += import_chain_length;

                    while (inbox_chain != NULL) {
                        SeqNode *node = inbox_chain;
                        inbox_chain = inbox_chain->head_next;
                        node->head_next = bulk_inbox_chain;
                        bulk_inbox_chain = node;
                    }
                }
                //printf("(flow) inbox component nodes: %"PRIuPTR"\n", round_inbox_nodes); // XXX: remove me
                //printf("(flow) inbox component related nodes: %"PRIuPTR"\n", round_import_nodes); fflush(stdout); // XXX: remove me

                SeqGraph *private_sgraph = psg_set.GetPrivateGraph(total_private_nodes);
                assert( private_sgraph->node_count == 0 );

                // first, add inbox nodes
                while (bulk_inbox_chain != NULL) {
                    SeqNode *node = bulk_inbox_chain;
                    bulk_inbox_chain = bulk_inbox_chain->head_next;
                    node->head_next = NULL; // XXX: needed?
                    private_sgraph->insertNodeAndUpdateColors(node); // XXX: updating colors too since multiple components!
                }

                // then, add imported nodes and fixup colors
                while (import_chain != NULL) {
                    SeqNode *node = import_chain;
                    import_chain = import_chain->head_next;
                    node->head_next = NULL; // XXX: needed?
                    private_sgraph->insertNodeAndUpdateColors(node);
                }

                flow_private_graph(private_sgraph, resident_sgraph, buckets, round);
                assert( private_sgraph->node_count == 0 );

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
        }

        /*if( g__FULL_STATISTICS ) {
          sg_stat_components(resident_sgraph, stdout);
          }*/

        // last, split the resident graph to outboxes // FIXME: should this be empty!?
        assert( resident_sgraph->node_count == 0 );
        //buckets->split(resident_sgraph);
    }
}

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

            printf("    flowing inbox: %s\n", itr->filename); fflush(stdout);
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

            if (file_length == 0) continue;  // otherwise, mmap will fail

            char *file_mmap = static_cast<char *>( mmap( 0, file_length, PROT_READ, MAP_PRIVATE, filedes, 0) );
            if (file_mmap == reinterpret_cast<char *>(-1)) {
                fprintf(stderr, "ERROR: failed to mmap file: %s\n", itr->filename);
                perror("REASON: ");
                exit(EXIT_FAILURE);
            }

            size_t file_offset = 0;
                
            flow_nodelist_t flowlist;
#ifdef VELOUR_TBB
            flowlist.reserve(1000000); // XXX: constant, should relate to value below
#endif

            uintptr_t round = 1;
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

void initializePartitionIndexFromInputFilename(file_object_vector *file_objects); // forward decl

void runFlow(file_object_vector *file_objects, KmerGraph *resident_kgraph, SeqGraph *resident_sgraph)
{
    file_object_vector loom_file_object;
    loom_file_object.push_back( file_objects->front() );
    assert( loom_file_object.front().filetype == LOOM );

    file_object_vector inbox_file_objects = *file_objects; // duplicate
    assert( inbox_file_objects.front().filetype == LOOM );
    inbox_file_objects.erase( inbox_file_objects.begin() ); // erase loom file entry
    assert( inbox_file_objects.empty() || inbox_file_objects.front().filetype == BUCKET );

    initializePartitionIndexFromInputFilename(&loom_file_object);

#ifdef VELOUR_TBB
    tbb::tick_count::interval_t loom_time;
    tbb::tick_count::interval_t kmerseq_time;
    tbb::tick_count::interval_t inbox_time;

    tbb::tick_count time0, time1;
    time0 = tbb::tick_count::now();
#endif // VELOUR_TBB

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

    SplitBuckets *outbox_buckets = new SplitBuckets(true);

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

    // TODO: when running with NDEBUG, delete this partition's inbox file
}


