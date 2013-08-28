//
// This file is part of Velour and is distributed under the University of
// Illinois Open Source License.  See LICENSE.txt for details.
//

//
// slicing2.cpp
//
//   identifies those nodes that are relevant to the frontier nodes ** of the preceeding partitions **,
//     thus indirectly identifying nodes that can be sent directly to other inboxes or the final bucket
//

#include "types.h"

uintptr_t g__SLICE2_FINAL_NODE_COUNT = 0;
uintptr_t g__SLICE2_NODE_COUNT = 0;

namespace {

static inline void slice_node(SeqGraph *graph, SeqNode *node, bool moving_right, unsigned incoming, Color frontier_color)
{
    if (node->slice_color <= frontier_color && (moving_right ? getNodeLeftSliceFlag(node,incoming) : getNodeRightSliceFlag(node,incoming))) {
        return;
    }

    //  TODO? reset needs resetNodeLeftSliceFlag etc
    //if (node->slice_color > frontier_color && (moving_right ? getNodeLeftSliceFlag(node,incoming) : getNodeRightSliceFlag(node,incoming))) { ... }

    // OPT? if both slice_color and frontier_color preceed currentPartitionIndex, then redundant?

    // claim the node for this frontier if we are a predecessor
    if (node->slice_color == 0 || node->slice_color > frontier_color) {
        if (node->slice_color > frontier_color) { resetNodeSliceFlags(node); } // FIXME: only reset some flags, somehow?
        node->slice_color = frontier_color;
    }

    unsigned left_count, right_count;
    node->getNeighborCounts(&left_count, &right_count);
  
    // first, slice the side we came in from
    unsigned back_count = (moving_right ? left_count : right_count);
    if (back_count > 1) {
        if (moving_right) {
            setNodeLeftSliceFlag(node, incoming);
        } else {
            setNodeRightSliceFlag(node, incoming);
        }
        counter_t *back_counters = (moving_right ? node->left_count : node->right_count);
        Color *back_colors = (moving_right ? node->left_color : node->right_color);
        counter_t back_value = cnorm(back_counters[incoming]);
        //Color back_color = back_colors[incoming];
        if (back_value < get_counter_max(back_counters) || back_value == 1) { // if minority (could be clipped), slice other edges
            for (unsigned j=0; j < 4; ++j) {
                bool isFlagged = (moving_right ? getNodeLeftSliceFlag(node,j) : getNodeRightSliceFlag(node,j));
                if (!isFlagged) {
                    if (back_colors[j] == 0) {
                        if (moving_right) { // don't mark a frontier edge
                            setNodeLeftSliceFlag(node,j);
                        } else {
                            setNodeRightSliceFlag(node,j);
                        }
                    }
                    if (back_counters[j] != 0 && back_colors[j] == 0) {
                        bool sense_changed;
                        SeqNode *next = graph->findNextNode(node, j, !moving_right, &sense_changed);
                        bool next_moving_right = !moving_right ^ sense_changed;
                        Nucleotide head_rightmost_base = node->sequence.GetHeadKmerRightmostBase(g__FULLKMER_LENGTH);
                        Nucleotide tail_leftmost_base = node->sequence.GetTailKmerLeftmostBase(g__FULLKMER_LENGTH);
                        int next_back_valid = !moving_right ? tail_leftmost_base : head_rightmost_base;
                        if (sense_changed) {
                            next_back_valid = COMPLEMENT(next_back_valid);
                        }
                        slice_node(graph,next,next_moving_right,next_back_valid,frontier_color);
                    }
                }
            }
        }
    }
    // useful for majority and single edges
    if (moving_right) {
        setNodeLeftSliceFlag(node, incoming);
    } else {
        setNodeRightSliceFlag(node, incoming);
    }

    // slice the flow-through direction if its a single edge
    unsigned next_count = (moving_right ? right_count : left_count);
    bool isNextFlagged = (moving_right ? getNodeRightSliceFlags(node) : getNodeLeftSliceFlags(node));
    if (!isNextFlagged && next_count == 1) {
        counter_t *counters = (moving_right ? node->right_count : node->left_count);
        int valid_direction = valid_single_successor(counters);
        Color *colors = (moving_right ? node->right_color : node->left_color);
        if (colors[valid_direction] == 0) {
            if (moving_right) { // don't mark a frontier edge
                setNodeRightSliceFlag(node, valid_direction);
            } else {
                setNodeLeftSliceFlag(node,valid_direction);
            }
            bool sense_changed;
            SeqNode *next = graph->findNextNode(node, valid_direction, moving_right, &sense_changed);
            bool next_moving_right = moving_right ^ sense_changed;
            Nucleotide head_rightmost_base = node->sequence.GetHeadKmerRightmostBase(g__FULLKMER_LENGTH);
            Nucleotide tail_leftmost_base = node->sequence.GetTailKmerLeftmostBase(g__FULLKMER_LENGTH);
            int next_back_valid = moving_right ? tail_leftmost_base : head_rightmost_base;
            if (sense_changed) {
                next_back_valid = COMPLEMENT(next_back_valid);
            }
            slice_node(graph,next,next_moving_right,next_back_valid,frontier_color);
        }
    }
 }

struct lambda_slice_marking {
    lambda_slice_marking(SeqGraph *graph, unsigned idx) : graph(graph), currentPartitionIndex(idx) {}
    void operator()(SeqNode *node)
    {
        if (node->connections == 0) return;

        /*if (getNodeLeftSliceFlags(node) != 0 && getNodeRightSliceFlags(node) != 0) { // marked on both sides
            return;
        }*/

        // check left side
        for (unsigned i=0; i < 4; ++i) {
            if (node->left_count[i] != 0 && node->left_color[i] != 0 && !getNodeLeftSliceFlag(node,i) ) {
                if (node->left_color[i] != (g__PARTITION_COUNT+1)) { // XXX: avoid final edges
                    slice_node(graph, node, true, i, node->left_color[i]);
                }
            }
        }

        // check right side
        for (unsigned i=0; i < 4; ++i) {
            if (node->right_count[i] != 0 && node->right_color[i] != 0 && !getNodeRightSliceFlag(node,i)) {
                if (node->right_color[i] != (g__PARTITION_COUNT+1)) { // XXX: avoid final edges
                    slice_node(graph, node, false, i, node->right_color[i]);
                }
            }
        }
    }
  private:
    SeqGraph *graph;
    unsigned currentPartitionIndex;
};

struct lambda_slice_transform_final {
    lambda_slice_transform_final(SeqGraph *graph, unsigned idx) : graph(graph), currentPartitionIndex(idx) {}
    void operator()(SeqNode *node)
    {
        if (node->slice_flags == 0) { // final node
            assert( node->slice_color == 0 );
            // check left side
            bool moving_right = false;
            for (unsigned i=0; i < 4; ++i) {
                if (node->left_count[i] != 0 && node->left_color[i] == 0) {
                    bool sense_changed;
                    SeqNode *next = graph->findNextNode(node, i, moving_right, &sense_changed);
                    if (next->slice_flags != 0) { // frontier of final
                        bool next_moving_right = moving_right ^ sense_changed;
                        Nucleotide head_rightmost_base = node->sequence.GetHeadKmerRightmostBase(g__FULLKMER_LENGTH);
                        Nucleotide tail_leftmost_base = node->sequence.GetTailKmerLeftmostBase(g__FULLKMER_LENGTH);
                        int next_back_valid = moving_right ? tail_leftmost_base : head_rightmost_base;
                        if (sense_changed) {
                            next_back_valid = COMPLEMENT(next_back_valid);
                        }
                        Color *next_back_colors = (next_moving_right ? next->left_color : next->right_color);
                        assert( next_back_colors[next_back_valid] == 0 );
                        // set to final bucket
                        node->left_color[i] = (g__PARTITION_COUNT+1);
                        counter_t *next_back_count = (next_moving_right ? next->left_count : next->right_count);
                        assert( next_back_count[next_back_valid] != 0 );
                        next_back_colors[next_back_valid] = (g__PARTITION_COUNT+1);
                    } else {
                        assert( next->slice_color == 0 );
                    }
                }
            }

            // check right side
            moving_right = true;
            for (unsigned i=0; i < 4; ++i) {
                if (node->right_count[i] != 0 && node->right_color[i] == 0) {
                    bool sense_changed;
                    SeqNode *next = graph->findNextNode(node, i, moving_right, &sense_changed);
                    if (next->slice_flags != 0) { // frontier of final
                        bool next_moving_right = moving_right ^ sense_changed;
                        Nucleotide head_rightmost_base = node->sequence.GetHeadKmerRightmostBase(g__FULLKMER_LENGTH);
                        Nucleotide tail_leftmost_base = node->sequence.GetTailKmerLeftmostBase(g__FULLKMER_LENGTH);
                        int next_back_valid = moving_right ? tail_leftmost_base : head_rightmost_base;
                        if (sense_changed) {
                            next_back_valid = COMPLEMENT(next_back_valid);
                        }
                        Color *next_back_colors = (next_moving_right ? next->left_color : next->right_color);
                        assert( next_back_colors[next_back_valid] == 0 );
                        // set to final bucket
                        node->right_color[i] = (g__PARTITION_COUNT+1);
                        counter_t *next_back_count = (next_moving_right ? next->left_count : next->right_count);
                        assert( next_back_count[next_back_valid] != 0 );
                        next_back_colors[next_back_valid] = (g__PARTITION_COUNT+1);
                    } else {
                        assert( next->slice_color == 0 );
                    }
                }
            }

            ++ g__SLICE2_FINAL_NODE_COUNT; // FIXME counting correctly?
        }
    }
  private:
    SeqGraph *graph;
    unsigned currentPartitionIndex;
};

struct lambda_slice_transform_others {
    lambda_slice_transform_others(SeqGraph *graph, unsigned idx) : graph(graph), currentPartitionIndex(idx) {}
    void operator()(SeqNode *node)
    {
        if (node->slice_color != 0) { // improved slicing: cut arc if slice colors disagree
            assert( node->slice_color != 0 );
            assert( node->slice_color != currentPartitionIndex);
            bool sliced = false;
            // check left side
            bool moving_right = false;
            for (unsigned i=0; i < 4; ++i) {
                if (node->left_count[i] != 0 && node->left_color[i] == 0) {
                    bool sense_changed;
                    SeqNode *next = graph->findNextNode(node, i, moving_right, &sense_changed);
                    assert( next != NULL );
                    assert( next->slice_color != 0 );
                    assert( next->slice_color != currentPartitionIndex);
                    if (next->slice_color != node->slice_color && (next->slice_color > currentPartitionIndex || node->slice_color > currentPartitionIndex)) { // frontier of slice -- checking cpidx s.t. color set only if node dne after split-emit
                        sliced = true;
                        bool next_moving_right = moving_right ^ sense_changed;
                        Nucleotide head_rightmost_base = node->sequence.GetHeadKmerRightmostBase(g__FULLKMER_LENGTH);
                        Nucleotide tail_leftmost_base = node->sequence.GetTailKmerLeftmostBase(g__FULLKMER_LENGTH);
                        int next_back_valid = moving_right ? tail_leftmost_base : head_rightmost_base;
                        if (sense_changed) {
                            next_back_valid = COMPLEMENT(next_back_valid);
                        }
                        Color *next_back_colors = (next_moving_right ? next->left_color : next->right_color);
                        assert( next_back_colors[next_back_valid] == 0 );
                        // slice edge by updating colors
                        node->left_color[i] = next->slice_color;
                        counter_t *next_back_count = (next_moving_right ? next->left_count : next->right_count);
                        assert( next_back_count[next_back_valid] != 0 );
                        next_back_colors[next_back_valid] = node->slice_color;
                    }
                }
            }

            // check right side
            moving_right = true;
            for (unsigned i=0; i < 4; ++i) {
                if (node->right_count[i] != 0 && node->right_color[i] == 0) {
                    bool sense_changed;
                    SeqNode *next = graph->findNextNode(node, i, moving_right, &sense_changed);
                    assert( next != NULL );
                    assert( next->slice_color != 0 );
                    assert( next->slice_color != currentPartitionIndex);
                    if (next->slice_color != node->slice_color && (next->slice_color > currentPartitionIndex || node->slice_color > currentPartitionIndex)) { // frontier of slice
                        sliced = true;
                        bool next_moving_right = moving_right ^ sense_changed;
                        Nucleotide head_rightmost_base = node->sequence.GetHeadKmerRightmostBase(g__FULLKMER_LENGTH);
                        Nucleotide tail_leftmost_base = node->sequence.GetTailKmerLeftmostBase(g__FULLKMER_LENGTH);
                        int next_back_valid = moving_right ? tail_leftmost_base : head_rightmost_base;
                        if (sense_changed) {
                            next_back_valid = COMPLEMENT(next_back_valid);
                        }
                        Color *next_back_colors = (next_moving_right ? next->left_color : next->right_color);
                        assert( next_back_colors[next_back_valid] == 0 );
                        // slice edge by updating colors
                        node->right_color[i] = next->slice_color;
                        counter_t *next_back_count = (next_moving_right ? next->left_count : next->right_count);
                        assert( next_back_count[next_back_valid] != 0 );
                        next_back_colors[next_back_valid] = node->slice_color;
                    }
                }
            }

            if (sliced) {
                ++ g__SLICE2_NODE_COUNT; // FIXME counting correctly?
            }
        }
    }
  private:
    SeqGraph *graph;
    unsigned currentPartitionIndex;
};

} // namespace: anonymous

/*void slice_component(SeqGraph *graph, Component& component)
{
    // reset the merged flags on the component nodes
    for (std::deque<SeqNode*>::iterator it=component.node_list.begin(); it != component.node_list.end(); ++it) {
        resetNodeMergingFlags<SeqNode>(*it);
    }

    // slice what we can -- start ONLY from frontier edges
    lambda_slice_marking slice_marking(graph);
    for (std::deque<SeqNode*>::iterator it=component.node_list.begin(); it != component.node_list.end(); ++it) {
        SeqNode *node = *it;
        slice_marking(node);
    }

    // set edge color to FINAL for all nodes without any merginging flags (and other end of edge)
    lambda_slice_transform slice_transform(graph);
    for (std::deque<SeqNode*>::iterator it=component.node_list.begin(); it != component.node_list.end(); ++it) {
        SeqNode *node = *it;
        slice_transform(node);
    }

    // restore the merged flags on the component nodes
    for (std::deque<SeqNode*>::iterator it=component.node_list.begin(); it != component.node_list.end(); ++it) {
        SeqNode *node = *it;
        resetNodeMergingFlags<SeqNode>(node);
        setNodeMerged<SeqNode>(node);
    }
}*/

void slice2_graph(SeqGraph *graph, unsigned currentPartitionIndex)
{
    graph->resetSliceFlagsAndSliceColor(); // XXX OPT: is this necessary, or do I keep these fields clean elsewhere?

    // slice what we can -- start ONLY from frontier edges
    lambda_slice_marking slice_marking(graph, currentPartitionIndex);
    sg_for_each(graph, slice_marking);

    // set edge color to FINAL for all nodes without any merginging flags
    lambda_slice_transform_final t_final(graph, currentPartitionIndex);
    sg_for_each(graph, t_final);

    // set edge colors for sliced out non-final nodes
    lambda_slice_transform_others t_others(graph, currentPartitionIndex);
    sg_for_each(graph, t_others);
}

void slice2_nodelist(SeqGraph *graph, unsigned currentPartitionIndex, flow_nodelist_t *nodelist)
{
    // XXX OPT: is this necessary, or do I keep these fields clean elsewhere?
    for (flow_nodelist_t::iterator it = nodelist->begin(); it != nodelist->end(); ++it) {
        resetNodeSliceFlags(*it);
        (*it)->slice_color = 0;
    }

    {
        // slice what we can -- start ONLY from frontier edges
        lambda_slice_marking slice_marking(graph, currentPartitionIndex);
        for (flow_nodelist_t::iterator it = nodelist->begin(); it != nodelist->end(); ++it) {
            SeqNode *node = *it;
            if (!isNodeDead(node)) { slice_marking(node); }
        }
    }

    {
        // set edge color to FINAL for all nodes without any merginging flags
        lambda_slice_transform_final t_final(graph, currentPartitionIndex);
        for (flow_nodelist_t::iterator it = nodelist->begin(); it != nodelist->end(); ++it) {
            SeqNode *node = *it;
            if (!isNodeDead(node)) { t_final(node); }
        }
    }

    {
        // set edge colors for sliced out non-final nodes
        lambda_slice_transform_others t_others(graph, currentPartitionIndex);
        for (flow_nodelist_t::iterator it = nodelist->begin(); it != nodelist->end(); ++it) {
            SeqNode *node = *it;
            if (!isNodeDead(node)) { t_others(node); }
        }
    }
}

#ifdef VELOUR_TBB
namespace {
struct lambda_slice_reset {
    typedef flow_nodelist_t::const_range_type::iterator iterator;
    lambda_slice_reset() {}
    void operator()( const flow_nodelist_t::const_range_type& range ) const {
        for (iterator it = range.begin() ; it != range.end(); ++it) {
            resetNodeSliceFlags(*it);
            (*it)->slice_color = 0;
        }
    }
};
} // namespace: anonymous

void slice2_parallel_nodelist(SeqGraph *graph, unsigned currentPartitionIndex, flow_nodelist_t *nodelist)
{
    // XXX OPT: is this necessary, or do I keep these fields clean elsewhere?
    tbb::parallel_for(nodelist->range(), lambda_slice_reset(), tbb::auto_partitioner());

    {
        // slice what we can -- start ONLY from frontier edges
        lambda_slice_marking slice_marking(graph, currentPartitionIndex);
        for (flow_nodelist_t::iterator it = nodelist->begin(); it != nodelist->end(); ++it) {
            SeqNode *node = *it;
            if (!isNodeDead(node)) { slice_marking(node); }
        }
    }

    {
        // set edge color to FINAL for all nodes without any merginging flags
        lambda_slice_transform_final t_final(graph, currentPartitionIndex);
        for (flow_nodelist_t::iterator it = nodelist->begin(); it != nodelist->end(); ++it) {
            SeqNode *node = *it;
            if (!isNodeDead(node)) { t_final(node); }
        }
    }

    {
        // set edge colors for sliced out non-final nodes
        lambda_slice_transform_others t_others(graph, currentPartitionIndex);
        for (flow_nodelist_t::iterator it = nodelist->begin(); it != nodelist->end(); ++it) {
            SeqNode *node = *it;
            if (!isNodeDead(node)) { t_others(node); }
        }
    }
}
#endif // VELOUR_TBB

