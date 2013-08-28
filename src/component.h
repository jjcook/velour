//
// This file is part of Velour and is distributed under the University of
// Illinois Open Source License.  See LICENSE.txt for details.
//

//
// component.h
//

#ifndef VELOUR_COMPONENT_H
#define VELOUR_COMPONENT_H

template<typename GraphType, typename NodeType>
class SerialComponent {
    public:
        typedef std::deque<NodeType *> nodelist_t;

        SerialComponent(GraphType *graph, NodeType *root) :
            min_frontier_index_(g__PARTITION_COUNT+1), frontier_edge_count_(0)
            { formComponent(graph, root); }

        unsigned get_min_frontier_index() { return min_frontier_index_; }
        uintptr_t get_node_count() { return nodes_.size(); }
        uintptr_t get_frontier_edge_count() { return frontier_edge_count_; }

    public:
        nodelist_t nodes_;

    private:
        void formComponent(GraphType *graph, NodeType *root);
        unsigned min_frontier_index_;
        uintptr_t frontier_edge_count_;
};

#ifdef VELOUR_TBB
template<typename GraphType, typename NodeType>
class ParallelComponent {
    public:
        static const unsigned UNCLAIMED = 0;
        static const unsigned RELEASING;
        typedef std::deque<NodeType *> nodelist_t;

        ParallelComponent(GraphType *graph, NodeType *root) { formComponent(graph, root); }

    public:
        nodelist_t nodes_;

    private:
        void formComponent(GraphType *graph, NodeType *root);
        bool claimNode(NodeType *node, Color old_value);
        void abortComponent(void);
};

template<typename GraphType, typename NodeType>
const unsigned ParallelComponent<GraphType, NodeType>::RELEASING = THREADID_MAX+1;

#endif // VELOUR_TBB


template<typename GraphType, typename NodeType>
void SerialComponent<GraphType, NodeType>::formComponent(GraphType *graph, NodeType *root)
{
    if (!isNodeDead<NodeType>(root) && !isNodeMerged<NodeType>(root)) {
        nodelist_t worklist;

        setNodeMerged<NodeType>(root);
        worklist.push_front(root);
        nodes_.push_front(root);
        while(!worklist.empty()) {
            NodeType *node = worklist.back();
            worklist.pop_back();

            // check right side
            counter_t *counters = node->right_count;
            Color *colors = node->right_color;
            for (int i = 0; i < 4; ++ i) {
                if (counters[i] != 0) {
                    if (colors[i] != 0) {
                        if (node->slice_color != 0) {
                            min_frontier_index_ = min(node->slice_color, min_frontier_index_);
                        } else {
                            min_frontier_index_ = min(colors[i], min_frontier_index_);
                        }
                        ++ frontier_edge_count_;
                    } else {
                        NodeType *next = graph->findNextNode(node, i, GO_RIGHT);
                        assert( !isNodeDead<NodeType>(next) );
                        if (!isNodeMerged<NodeType>(next)) {
                            setNodeMerged<NodeType>(next);
                            worklist.push_front(next);
                            nodes_.push_front(next);
                        }
                    }
                }
            }

            // check left side
            counters = node->left_count;
            colors = node->left_color;
            for (int i = 0; i < 4; ++ i) {
                if (counters[i] != 0) {
                    if (colors[i] != 0) {
                        if (node->slice_color != 0) {
                            min_frontier_index_ = min(node->slice_color, min_frontier_index_);
                        } else {
                            min_frontier_index_ = min(colors[i], min_frontier_index_);
                        }
                        ++ frontier_edge_count_;
                    } else {
                        NodeType *next = graph->findNextNode(node, i, GO_LEFT);
                        assert( !isNodeDead<NodeType>(next) );
                        if (!isNodeMerged<NodeType>(next)) {
                            setNodeMerged<NodeType>(next);
                            worklist.push_front(next);
                            nodes_.push_front(next);
                        }
                    }
                }
            }
        }
        assert( min_frontier_index_ > 0 && min_frontier_index_ <= (g__PARTITION_COUNT+1) );
    }
}

#ifdef VELOUR_TBB
template<typename GraphType, typename NodeType>
inline bool ParallelComponent<GraphType, NodeType>::claimNode(NodeType *node, Color old_value)
{
    bool claimed = (ATOMIZE(node->claim_tid).compare_and_swap(tls_thread_index, old_value) == old_value);
    return claimed;
}

template<typename GraphType, typename NodeType>
inline void ParallelComponent<GraphType, NodeType>::abortComponent(void)
{
    // reset claim_tid in case we re-enter the same component again
    for (typename nodelist_t::iterator it = nodes_.begin() ; it != nodes_.end() ; ++it) {
        NodeType *node = *it;
        ATOMIZE(node->claim_tid).compare_and_swap(UNCLAIMED, tls_thread_index); // reset if still claimed by me
        assert( ATOMIC_LOAD(node->claim_tid) != tls_thread_index ); // TODO: upgrade this to a do-while loop since so important???
    }
    nodes_.clear();
    ATOMIC_STORE(*tls_thread_abort) = false; // reset abort flag
}

// overloaded (instead of templatized) find node function
// -- kmer graph find is non-atomic
// -- seq graph find is atomic
static inline KmerNode * findNextNode_helper(KmerGraph *graph, KmerNode *node, int i, int direction);
static inline SeqNode * findNextNode_helper(SeqGraph *graph, SeqNode *node, int i, int direction);

template<typename GraphType, typename NodeType>
void ParallelComponent<GraphType, NodeType>::formComponent(GraphType *graph, NodeType *root)
{
    assert( nodes_.empty() );
    if (!atomic_isNodeDead<NodeType>(root) && (ATOMIC_LOAD(root->claim_tid) == UNCLAIMED)) {
           
        if (!claimNode(root, UNCLAIMED)) { return; } // root stolen
        nodes_.push_front(root);

        if (atomic_isNodeDead<NodeType>(root)) { abortComponent(); return; }

        nodelist_t worklist;

        worklist.push_front(root);
        while(!worklist.empty()) {
            NodeType *node = worklist.back();
            worklist.pop_back();

            // check right side
            counter_t *counters = node->right_count;
            Color *colors = node->right_color;
            for (int i = 0; i < 4; ++ i) {
                if (counters[i] != 0) {
                    if (colors[i] == 0) {
                        NodeType *next = findNextNode_helper(graph, node, i, GO_RIGHT);
                        if (next == NULL || atomic_isNodeDead<NodeType>(next) || ATOMIC_LOAD(*tls_thread_abort)) { abortComponent(); return; }
                        if (ATOMIC_LOAD(next->claim_tid) != tls_thread_index) {
                            Color old_value = UNCLAIMED;
                            while(1) {
                                if (claimNode(next, old_value)) {
#ifdef VERIFY_THREADSAFE
                                    bool duplicate = ( find(nodes_.begin(), nodes_.end(), next) != nodes_.end() );
                                    assert( !duplicate );
#endif // VERIFY_THREADSAFE
                                    nodes_.push_front(next);
                                    if (ATOMIC_LOAD(*tls_thread_abort) || atomic_isNodeDead<NodeType>(next)) { abortComponent(); return; }
                                    worklist.push_front(next);
                                    break;
                                } else {
                                    old_value = ATOMIC_LOAD(next->claim_tid);
                                    assert( old_value != tls_thread_index );
                                    if (old_value == RELEASING) {
                                        while (ATOMIC_LOAD(next->claim_tid) == RELEASING) {} // spin until released
                                        old_value = UNCLAIMED; // reset and try to claim the node
                                    } else if (old_value > tls_thread_index || ATOMIC_LOAD(*tls_thread_abort) ) { // abort: release claimed nodes
                                        abortComponent();
                                        return;
                                    } else if (old_value != UNCLAIMED) { // implicit: old_value < tls_thread_index
                                        // abort other thread
                                        if (ATOMIC_LOAD(*thread_aborts[old_value]) == false) { // check if false to avoid cache-line ping-pong
                                            ATOMIC_STORE(*thread_aborts[old_value]) = true;
                                        }
                                    }
                                    // loop waiting for other thread to release claimed nodes
                                }
                            }
                        }
#ifdef VERIFY_THREADSAFE
                        else {
                            bool found = ( find(nodes_.begin(), nodes_.end(), next) != nodes_.end() );
                            assert( found );
                        }
#endif // VERIFY_THREADSAFE
                    }
                }
            }

            // check left side
            counters = node->left_count;
            colors = node->left_color;
            for (int i = 0; i < 4; ++ i) {
                if (counters[i] != 0) {
                    if (colors[i] == 0) {
                        NodeType *next = findNextNode_helper(graph, node, i, GO_LEFT);
                        if (next == NULL || atomic_isNodeDead<NodeType>(next) || ATOMIC_LOAD(*tls_thread_abort)) { abortComponent(); return; }
                        if (ATOMIC_LOAD(next->claim_tid) != tls_thread_index) {
                            Color old_value = UNCLAIMED;
                            while(1) {
                                if (claimNode(next, old_value)) {
#ifdef VERIFY_THREADSAFE
                                    bool duplicate = ( find(nodes_.begin(), nodes_.end(), next) != nodes_.end() );
                                    assert( !duplicate );
#endif // VERIFY_THREADSAFE
                                    nodes_.push_front(next);
                                    if (ATOMIC_LOAD(*tls_thread_abort) || atomic_isNodeDead<NodeType>(next)) { abortComponent(); return; }
                                    worklist.push_front(next);
                                    break;
                                } else {
                                    old_value = ATOMIC_LOAD(next->claim_tid);
                                    assert( old_value != tls_thread_index );
                                    if (old_value == RELEASING) {
                                        while (ATOMIC_LOAD(next->claim_tid) == RELEASING) {} // spin until released
                                        old_value = UNCLAIMED; // reset and try to claim the node
                                    } else if (old_value > tls_thread_index || ATOMIC_LOAD(*tls_thread_abort) ) { // abort: release claimed nodes
                                        abortComponent();
                                        return;
                                    } else if (old_value != UNCLAIMED) { // implicit: old_value < tls_thread_index
                                        // abort other thread
                                        if (ATOMIC_LOAD(*thread_aborts[old_value]) == false) { // check if false to avoid cache-line ping-pong
                                            ATOMIC_STORE(*thread_aborts[old_value]) = true;
                                        }
                                    }
                                    // loop waiting for other thread to release claimed nodes
                                }
                            }
                        }
#ifdef VERIFY_THREADSAFE
                        else {
                            bool found = ( find(nodes_.begin(), nodes_.end(), next) != nodes_.end() );
                            assert( found );
                        }
#endif // VERIFY_THREADSAFE
                    }
                }
            }
        }
    }
}
#endif // VELOUR_TBB

#endif // VELOUR_COMPONENT_H
