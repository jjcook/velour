//
// This file is part of Velour and is distributed under the University of
// Illinois Open Source License.  See LICENSE.txt for details.
//

// 
// sg_tipClipping.cpp
//
//   sequence graph tip clipping
//

#include "types.h"

// TODO: configurable parameter?
#define TIP_REMOVAL_LENGTH (2 * g__FULLKMER_LENGTH)


static counter_t g_minimum_edge_weight = 2; // TODO: configurable parameter?


//********************************************************************************
//******************************  Tip Clipping  **********************************
//********************************************************************************


//static unsigned EXTENSIONS = 0;
static uintptr_t TIP_DETACHED = 0;
//static unsigned TIP_DEALLOCATED = 0;
static uintptr_t NODES_DEALLOCATED = 0;
static uintptr_t DISCONNECTED_CHUNKS = 0;
static uintptr_t BOUNDARY_CUTS = 0;

static int REMOVE_NOW_VALUE = 2;
static int REMOVE_NEXT_VALUE = 2; // 4;

static void 
flip_to_remove_value() {
  int temp = REMOVE_NOW_VALUE;
  REMOVE_NOW_VALUE = REMOVE_NEXT_VALUE;
  REMOVE_NEXT_VALUE = temp;
}

// TODO: small nodes vs large nodes

// we know that this "node" is a node that has either no neighbors on
// one side -or- connects to a chain of nodes that leads to a single node.
static bool
check_for_tip_removal(SeqGraph *graph, SeqNode *node, bool moving_right, bool singlecopy_tip, unsigned length)
{
    // XXX: performance optimization not correctness, right?
#ifdef VELOUR_TBB
    if (atomic_isNodeDead<SeqNode>(node)) { return false; }
#else
    if (isNodeDead<SeqNode>(node)) { return false; }
#endif

	length += node->sequence.get_length() - g__FULLKMER_LENGTH + 1;
    singlecopy_tip &= (isNodeWeight_None(node) || isNodeWeight_AllSingleCopy(node));

	if (length >= TIP_REMOVAL_LENGTH && !singlecopy_tip) { // FIXME? singlecopy_tip might be incorrect as didn't consider forward arc???
		return false;  // can't remove a tip that is too long
	}

	// find out if this tip can be extended.  This has two conditions:
	// CONDITION 1: there is a single neighbor in the extension direction
#ifdef VELOUR_TBB
    union {
        counter_t count[4];
        four_counter_t side;
    } local_counters;
    local_counters.side = moving_right ? ATOMIC_LOAD(node->right_side) : ATOMIC_LOAD(node->left_side);
	counter_t   *counters = local_counters.count;
#else
	counter_t   *counters = moving_right ? node->right_count : node->left_count;
#endif

	int valid_dir = valid_single_successor(counters);

	if (valid_dir == MULTIPLE_SUCCESSORS) {  
		return false; // can't remove a tip with two outgrowths (the base of a Y)
	}

	if (valid_dir == NO_SUCCESSORS) {
#ifdef VELOUR_TBB
        atomic_setNodeDead<SeqNode>(node);
        ATOMIC_STORE(node->connections) = 0;
        ATOMIC_STORE(node->left_side_colors) = 0;
        ATOMIC_STORE(node->right_side_colors) = 0;
        ATOMIZE(DISCONNECTED_CHUNKS).fetch_and_increment();
#else
        setNodeDead<SeqNode>(node);
		node->connections = 0;
        node->left_side_colors = 0;
        node->right_side_colors = 0;
        ++ DISCONNECTED_CHUNKS;
#endif
		// printf("a disconnected chunk of length %d\n", length);
		return true;
	}

	// a single successor
	// CONDITION 2: the node in the extension direction has a single
	// neighbor in the direction facing this node
	assert((valid_dir >= 0) && (valid_dir < 4));
	bool sense_changed;
#ifdef VELOUR_TBB
    union {
        color_t color[4];
        four_color_t side;
    } local_colors;
    local_colors.side = moving_right ? ATOMIC_LOAD(node->right_side_colors) : ATOMIC_LOAD(node->left_side_colors);
    color_t *colors = local_colors.color;
#else
    color_t *colors = moving_right ? node->right_color : node->left_color;
#endif
    SeqNode *next = NULL;
    if (colors[valid_dir] == 0) {
        next = graph->findNextNode(node, valid_dir, moving_right, &sense_changed);
    } else {
        return false;
	}
    assert( next != NULL );

	bool next_moving_right = moving_right ^ sense_changed;
#ifdef VELOUR_TBB
    union {
        counter_t count[4];
        four_counter_t side;
    } local_back_counters;
    local_back_counters.side = next_moving_right ? ATOMIC_LOAD(next->left_side) : ATOMIC_LOAD(next->right_side);
	counter_t *next_back_counters = local_back_counters.count;
#else
	counter_t *next_back_counters = next_moving_right ? next->left_count : next->right_count;
#endif

	int next_back_valid = valid_single_successor(next_back_counters);

	if (next_back_valid != MULTIPLE_SUCCESSORS) {
        bool forceclip = (abs(next_back_counters[next_back_valid]) == CLIP_SINGLECOPY_COUNTER_VALUE);
        singlecopy_tip &= cnorm(counters[valid_dir]) == 1;
		bool remove = check_for_tip_removal(graph, next, next_moving_right, singlecopy_tip, length);
		if (remove || singlecopy_tip || forceclip) {
#ifdef VELOUR_TBB
            atomic_setNodeDead<SeqNode>(node);
            ATOMIC_STORE(node->connections) = 0;
            ATOMIC_STORE(node->left_side_colors) = 0;
            ATOMIC_STORE(node->right_side_colors) = 0;
#else
            setNodeDead<SeqNode>(node);
			node->connections = 0;
            node->left_side_colors = 0;
            node->right_side_colors = 0;
#endif
        }
        
        if (!remove && (singlecopy_tip || forceclip)) {
#ifdef VELOUR_TBB
            ATOMIC_STORE( (next_moving_right ? next->left_count : next->right_count)[next_back_valid] ) = 0;
            ATOMIC_STORE( (next_moving_right?next->left_color:next->right_color)[next_back_valid] ) = 0;
#else
            next_back_counters[next_back_valid] = 0;
            (next_moving_right?next->left_color:next->right_color)[next_back_valid] = 0;
#endif
        }
        return remove || singlecopy_tip || forceclip;
	}

	// can't extend.  Need to determine if we should trim.
	counter_t max_value = get_counter_max(next_back_counters);
	counter_t value = abs(counters[valid_dir]);
    if (value == CLIP_SINGLECOPY_COUNTER_VALUE) value = 1;
	if ((value >= g_minimum_edge_weight) && (value == max_value)) {
		// printf("dominant connection: length %d, weight %d\n", length, abs(counters[valid_dir]));
		return false;   // can't remove the dominant connection
	}

  // check for multiple single-copy predecessors in next node; mark those arcs as potential clips
  if (value == 1 && value == max_value) {
    for (int i = 0 ; i < 4 ; ++ i) {
        if (next_back_counters[i] == 1) {
#ifdef VELOUR_TBB
            ATOMIC_STORE( (next_moving_right ? next->left_count : next->right_count)[i] ) = CLIP_SINGLECOPY_COUNTER_VALUE;
#else
            next_back_counters[i] = CLIP_SINGLECOPY_COUNTER_VALUE;
#endif
        } else if (next_back_counters[i] == -1) {
#ifdef VELOUR_TBB
            ATOMIC_STORE( (next_moving_right ? next->left_count : next->right_count)[i] ) = - CLIP_SINGLECOPY_COUNTER_VALUE;
#else
            next_back_counters[i] = - CLIP_SINGLECOPY_COUNTER_VALUE;
#endif
        }
    }
  }

	// if we are clipping the tip, then mark the current "node" for
	// removal and disconnect the "next" node from the tip.  

/*
#ifdef VERIFY
  verify_node(node, HASHTABLE, g__FULLKMER_LENGTH);     // FINE, BEFORE WE LET IT GET INCONSISTENT.
  verify_node(next, HASHTABLE, g__FULLKMER_LENGTH);
#endif
*/

	// Note: we remove the link from "next", but not from "node".  The
	// reason why is that we don't want code to start on the other side
	// thinking that it is a tip that might need removal.  This means 
	// that from "node's" perspective, the graph looks inconsistent, but
	// that is okay, since we're going to remove "node" shortly.
#ifdef VELOUR_TBB
    atomic_setNodeDead<SeqNode>(node);
    ATOMIC_STORE(node->connections) = 0;
    ATOMIC_STORE(node->left_side_colors) = 0;
    ATOMIC_STORE(node->right_side_colors) = 0;
#else
    setNodeDead<SeqNode>(node);
	node->connections = 0;
    node->left_side_colors = 0;
    node->right_side_colors = 0;
#endif
	// printf("removal candidate: length %d, weight %d\n", length, abs(counters[valid_dir]));

	Nucleotide head_rightmost_base = node->sequence.GetHeadKmerRightmostBase(g__FULLKMER_LENGTH);
	Nucleotide tail_leftmost_base = node->sequence.GetTailKmerLeftmostBase(g__FULLKMER_LENGTH);
	next_back_valid = moving_right ? tail_leftmost_base : head_rightmost_base;

	if (sense_changed) { next_back_valid ^= 0x3; } // FIXME: i.e. complement??
           
#ifdef VELOUR_TBB 
    ATOMIC_STORE( (next_moving_right?next->left_color:next->right_color)[next_back_valid] ) = 0;
#else
    (next_moving_right?next->left_color:next->right_color)[next_back_valid] = 0;
#endif
  
#ifdef VELOUR_TBB
    ATOMIC_STORE( (next_moving_right ? next->left_count : next->right_count)[next_back_valid] ) = 0;
#else
	next_back_counters[next_back_valid] = 0;
#endif

#ifdef VELOUR_TBB
    ATOMIZE(TIP_DETACHED).fetch_and_increment();
#else
	++ TIP_DETACHED;
#endif

/*
#ifdef VERIFY
	// verify_node(node, g__FULLKMER_LENGTH);  // THIS VALIDATION WOULD FAIL, SINCE WE DON'T MODIFY NODE
	verify_node(next, HASHTABLE, g__FULLKMER_LENGTH);
#endif
*/

  	return true;
}

namespace {

struct serial_remove_tips_functor {
    serial_remove_tips_functor(SeqGraph *graph, bool *modified) : graph(graph), modified(modified) {}
    void operator()(SeqNode *node) {
#ifdef SMALL_NODES
        unsigned left_count = 0, right_count = 0;
        node->getNeighborCounts(&left_count, &right_count);
#else
        // count the number of left and right neighbors
        unsigned connections = trav->connections;
        unsigned left_count = CONNECTION_COUNT_MAP[LEFT(connections)];
        unsigned right_count = CONNECTION_COUNT_MAP[RIGHT(connections)];
#endif
        unsigned total_count = left_count + right_count;

        if (total_count <= 1) {  // is a tip; check if it should be removed
            *modified |= check_for_tip_removal(graph, node, (left_count == 0), (isNodeWeight_None(node) || isNodeWeight_AllSingleCopy(node)), 0);
        }
    }
  private:
    SeqGraph *graph;
    bool *modified;
};

#ifdef VELOUR_TBB
struct parallel_remove_tips_functor {
    typedef flow_nodelist_t::const_range_type::iterator iterator;
    parallel_remove_tips_functor(SeqGraph *graph, bool *modified) : graph(graph), modified(modified) {}
    void operator()( const flow_nodelist_t::const_range_type& range ) const {
        for (iterator it = range.begin() ; it != range.end() ; ++it) {
            SeqNode *node = *it;

#ifdef SMALL_NODES
            unsigned left_count = 0, right_count = 0;
            node->getNeighborCounts(&left_count, &right_count);
#else
            // count the number of left and right neighbors
            unsigned connections = trav->connections;
            unsigned left_count = CONNECTION_COUNT_MAP[LEFT(connections)];
            unsigned right_count = CONNECTION_COUNT_MAP[RIGHT(connections)];
#endif
            unsigned total_count = left_count + right_count;

            if (total_count <= 1) {  // is a tip; check if it should be removed
                bool was_modified = check_for_tip_removal(graph, node, (left_count == 0), (isNodeWeight_None(node) || isNodeWeight_AllSingleCopy(node)), 0);
                if (was_modified) { ATOMIC_STORE(*modified) = true; }
            }
        }
    }
  private:
    SeqGraph *graph;
    bool *modified;
};
#endif // VELOUR_TBB

} // namespace: anonymous

// serial remove_tips
void sg_remove_tips(SeqGraph *graph, bool silent)
{
	bool modified = true;
	unsigned pass = 1;

    if (g__NO_TIP_CLIPPING) return;

	while (modified) {
        uintptr_t start_node_count = graph->node_count;
	   	modified = false;

        if (!silent) {
            printf("==================== SG: TIP REMOVAL PASS %d ====================\n", pass);
        }
     
        serial_remove_tips_functor f(graph, &modified);
        sg_for_each_mutable(graph, f);

        uintptr_t nodes_deallocated_this_time = (start_node_count - graph->node_count);

        if (!silent) {
            printf("sg_tip_removal: %"PRIuPTR" tips detached, %"PRIuPTR" nodes deallocated, %"PRIuPTR" total nodes deallocated, %"PRIuPTR" floating chunks, %"PRIuPTR" boundary arcs cut\n", 
			  TIP_DETACHED, nodes_deallocated_this_time, NODES_DEALLOCATED, DISCONNECTED_CHUNKS, BOUNDARY_CUTS);
        }

        ++pass;

#ifndef SMALL_PRENODES
	    flip_to_remove_value();
#endif
    }
}

// serial remove_tips
void sg_nodelist_remove_tips(SeqGraph *graph, bool silent, flow_nodelist_t *nodelist)
{
	bool modified = true;
	unsigned pass = 1;

    if (g__NO_TIP_CLIPPING) return;

	while (modified) {
        uintptr_t start_node_count = graph->node_count;
        modified = false;

        if (!silent) {
            printf("==================== SG: TIP REMOVAL PASS %d ====================\n", pass);
        }

        serial_remove_tips_functor f(graph, &modified);
        std::for_each(nodelist->begin(), nodelist->end(), f);

        uintptr_t nodes_deallocated_this_time = (start_node_count - graph->node_count);

        if (!silent) {
            printf("sg_tip_removal: %"PRIuPTR" tips detached, %"PRIuPTR" nodes deallocated, %"PRIuPTR" total nodes deallocated, %"PRIuPTR" floating chunks, %"PRIuPTR" boundary arcs cut\n", 
			  TIP_DETACHED, nodes_deallocated_this_time, NODES_DEALLOCATED, DISCONNECTED_CHUNKS, BOUNDARY_CUTS);
        }

        ++pass;

#ifndef SMALL_PRENODES
	    flip_to_remove_value();
#endif
    }
}

#ifdef VELOUR_TBB
// parallel remove_tips
void sg_parallel_nodelist_remove_tips(SeqGraph *graph, bool silent, flow_nodelist_t *nodelist)
{
	bool modified = true;
	unsigned pass = 1;

    if (g__NO_TIP_CLIPPING) return;

	while (modified) {
        uintptr_t start_node_count = graph->node_count;
        modified = false;

        if (!silent) {
            printf("==================== SG: TIP REMOVAL PASS %d ====================\n", pass);
        }

        parallel_remove_tips_functor pf(graph, &modified);
        tbb::parallel_for(nodelist->range(), pf, tbb::auto_partitioner());

        uintptr_t nodes_deallocated_this_time = (start_node_count - graph->node_count);

        if (!silent) {
            printf("sg_tip_removal: %"PRIuPTR" tips detached, %"PRIuPTR" nodes deallocated, %"PRIuPTR" total nodes deallocated, %"PRIuPTR" floating chunks, %"PRIuPTR" boundary arcs cut\n", 
			  TIP_DETACHED, nodes_deallocated_this_time, NODES_DEALLOCATED, DISCONNECTED_CHUNKS, BOUNDARY_CUTS);
        }

        ++pass;

#ifndef SMALL_PRENODES
	    flip_to_remove_value();
#endif
    }
}
#endif // VELOUR_TBB

