//
// This file is part of Velour and is distributed under the University of
// Illinois Open Source License.  See LICENSE.txt for details.
//

//
// slicing.cpp
//
//   identifies those nodes that are relevant to the frontier nodes,
//     thus indirectly identifying nodes that can be sent directly to the final bucket
//

#include "types.h"

uintptr_t g__SLICE_NODE_COUNT = 0;

namespace {

static inline void slice_node(SeqGraph *graph, SeqNode *node, bool moving_right, unsigned incoming)
{
    if ((moving_right ? getNodeLeftFlag<SeqNode>(node,incoming) : getNodeRightFlag<SeqNode>(node,incoming))) {
        return;
    }

    unsigned left_count, right_count;
    node->getNeighborCounts(&left_count, &right_count);
  
    // first, slice the side we came in from
    unsigned back_count = (moving_right ? left_count : right_count);
    if (back_count > 1) {
        if (moving_right) {
            setNodeLeftFlag<SeqNode>(node, incoming);
        } else {
            setNodeRightFlag<SeqNode>(node, incoming);
        }
        counter_t *back_counters = (moving_right ? node->left_count : node->right_count);
        Color *back_colors = (moving_right ? node->left_color : node->right_color);
        counter_t back_value = cnorm(back_counters[incoming]);
        Color back_color = back_colors[incoming];
        if (back_value < get_counter_max(back_counters) || back_value == 1) { // if minority (could be clipped), slice other edges
            for (unsigned j=0; j < 4; ++j) {
                bool isFlagged = (moving_right ? getNodeLeftFlag<SeqNode>(node,j) : getNodeRightFlag<SeqNode>(node,j));
                if (!isFlagged) {
                    if (moving_right) {
                        setNodeLeftFlag<SeqNode>(node,j);
                    } else {
                        setNodeRightFlag<SeqNode>(node,j);
                    }
                    if (back_counters[j] != 0 && back_colors[j] == 0) {
                        bool sense_changed;
                        SeqNode *next = graph->findNextNode(node, j, !moving_right, &sense_changed);
                        bool next_moving_right = !moving_right ^ sense_changed;
                        Nucleotide head_rightmost_base = seq_get_head_rightmost_base(&node->sequence, g__FULLKMER_LENGTH);
                        Nucleotide tail_leftmost_base = seq_get_tail_leftmost_base(&node->sequence, g__FULLKMER_LENGTH);
                        int next_back_valid = !moving_right ? tail_leftmost_base : head_rightmost_base;
                        if (sense_changed) {
                            next_back_valid = COMPLEMENT(next_back_valid);
                        }
                        slice_node(graph,next,next_moving_right,next_back_valid);
                    }
                }
            }
        }
    }
    // useful for majority and single edges
    if (moving_right) {
        setNodeLeftFlag<SeqNode>(node, incoming);
    } else {
        setNodeRightFlag<SeqNode>(node, incoming);
    }

    // slice the flow-through direction if its a single edge
    unsigned next_count = (moving_right ? right_count : left_count);
    bool isNextFlagged = (moving_right ? getNodeRightFlags<SeqNode>(node) : getNodeLeftFlags<SeqNode>(node));
    if (!isNextFlagged && next_count == 1) {
        counter_t *counters = (moving_right ? node->right_count : node->left_count);
        int valid_direction = valid_single_successor(counters);
        if (moving_right) {
            setNodeRightFlag<SeqNode>(node, valid_direction);
        } else {
            setNodeLeftFlag<SeqNode>(node,valid_direction);
        }
        Color *colors = (moving_right ? node->right_color : node->left_color);
        if (colors[valid_direction] == 0) {
            bool sense_changed;
            SeqNode *next = graph->findNextNode(node, valid_direction, moving_right, &sense_changed);
            bool next_moving_right = moving_right ^ sense_changed;
            Nucleotide head_rightmost_base = seq_get_head_rightmost_base(&node->sequence, g__FULLKMER_LENGTH);
            Nucleotide tail_leftmost_base = seq_get_tail_leftmost_base(&node->sequence, g__FULLKMER_LENGTH);
            int next_back_valid = moving_right ? tail_leftmost_base : head_rightmost_base;
            if (sense_changed) {
                next_back_valid = COMPLEMENT(next_back_valid);
            }
            slice_node(graph,next,next_moving_right,next_back_valid);
        }
    }
 }

struct lambda_slice_marking {
    lambda_slice_marking(SeqGraph *graph) : graph(graph) {}
    void operator()(SeqNode *node)
    {
        if (node->connections == 0) return;

        /*if (getNodeLeftFlags<SeqNode>(node) != 0 && getNodeRightFlags<SeqNode>(node) != 0) { // marked on both sides
            return;
        }*/

        // check left side
        for (unsigned i=0; i < 4; ++i) {
            if (node->left_count[i] != 0 && node->left_color[i] != 0 && !getNodeLeftFlag<SeqNode>(node,i) ) {
                slice_node(graph, node, true, i);
            }
        }

        // check right side
        for (unsigned i=0; i < 4; ++i) {
            if (node->right_count[i] != 0 && node->right_color[i] != 0 && !getNodeRightFlag<SeqNode>(node,i)) {
                slice_node(graph, node, false, i);
            }
        }
    }
  private:
    SeqGraph *graph;
};

struct lambda_slice_transform {
    lambda_slice_transform(SeqGraph *graph) : graph(graph) {}
    void operator()(SeqNode *node)
    {
        if (node->flags == 0) { // final node
            // check left side
            bool moving_right = false;
            for (unsigned i=0; i < 4; ++i) {
                if (node->left_count[i] != 0 && node->left_color[i] == 0) {
                    bool sense_changed;
                    SeqNode *next = graph->findNextNode(node, i, moving_right, &sense_changed);
                    if (next->flags != 0) { // frontier of final
                        bool next_moving_right = moving_right ^ sense_changed;
                        Nucleotide head_rightmost_base = seq_get_head_rightmost_base(&node->sequence, g__FULLKMER_LENGTH);
                        Nucleotide tail_leftmost_base = seq_get_tail_leftmost_base(&node->sequence, g__FULLKMER_LENGTH);
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
                    }
                }
            }

            // check right side
            moving_right = true;
            for (unsigned i=0; i < 4; ++i) {
                if (node->right_count[i] != 0 && node->right_color[i] == 0) {
                    bool sense_changed;
                    SeqNode *next = graph->findNextNode(node, i, moving_right, &sense_changed);
                    if (next->flags != 0) { // frontier of final
                        bool next_moving_right = moving_right ^ sense_changed;
                        Nucleotide head_rightmost_base = seq_get_head_rightmost_base(&node->sequence, g__FULLKMER_LENGTH);
                        Nucleotide tail_leftmost_base = seq_get_tail_leftmost_base(&node->sequence, g__FULLKMER_LENGTH);
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
                    }
                }
            }

            ++ g__SLICE_NODE_COUNT;
        }
    }
  private:
    SeqGraph *graph;
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

void slice_graph(SeqGraph *graph)
{
    graph->resetFlags();

    // slice what we can -- start ONLY from frontier edges
    lambda_slice_marking slice_marking(graph);
    sg_for_each(graph, slice_marking);

    // set edge color to FINAL (and other side of edge) for all nodes without any merginging flags
    lambda_slice_transform slice_transform(graph);
    sg_for_each(graph, slice_transform);
}

