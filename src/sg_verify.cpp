//
// This file is part of Velour and is distributed under the University of
// Illinois Open Source License.  See LICENSE.txt for details.
//

//
// sg_verify.cpp
//
//   comprehensive sequence graph consistentcy verification
//

#include "types.h"
#include "seqGraph.h"

static inline void verify_node_helper(SeqGraph *graph, SeqNode *node);

struct VerifyNode_functor {
    VerifyNode_functor(SeqGraph *graph) : graph(graph) {}
    void operator()(SeqNode *node) { verify_node_helper(graph, node); }
  private:
    SeqGraph *graph;
};

void SeqGraph::verify(bool silent)
{
    if (!silent) {
        printf("======== SEQUENCE GRAPH VERIFY ========\n");
        fflush(stdout);
    }

    VerifyNode_functor f(this);
    sg_for_each(this, f);
}

void SeqGraph::verify_node(SeqNode *node)
{
    VerifyNode_functor f(this);
    f(node);
}


static inline void verify_node_helper(SeqGraph *graph, SeqNode *node)
{
    assert( !isNodeDead<SeqNode>(node) );
    // TODO check other fields, etc -- e.g. flags consistent?, slice_color consistent?

    //
    // check connectivity and colors
    //
	Nucleotide head_rightmost_base = node->sequence.GetHeadKmerRightmostBase(g__FULLKMER_LENGTH);
	Nucleotide tail_leftmost_base = node->sequence.GetTailKmerLeftmostBase(g__FULLKMER_LENGTH);
	for (int i = 0 ; i < 4 ; ++ i) {
	   	// check on the left side
		int count = node->left_count[i];
		if (count != 0) {
            bool sense_changed;
			SeqNode *next = graph->findNextNode(node, i, GO_LEFT, &sense_changed);
			assert( next != NULL || (g__PSEUDO_NODES_PRESENT && node->left_color[i] != 0) );
			if (next != NULL) {
                assert( node->left_color[i] == 0 );
                int next_back = sense_changed ? COMPLEMENT(head_rightmost_base) : head_rightmost_base;
                counter_t *next_count = sense_changed ? next->left_count : next->right_count;
                Color *next_color = sense_changed ? next->left_color : next->right_color;
                assert( cnorm(next_count[next_back]) == cnorm(node->left_count[i]) );
                assert( next_color[next_back] == 0 );

                assert( !isNodeDead<SeqNode>(next) );
			}
		} else {
            assert( node->left_color[i] == 0 );
        }

		// check on the right side
		count = node->right_count[i];
		if (count != 0) {
            bool sense_changed;
			SeqNode *next = graph->findNextNode(node, i, GO_RIGHT, &sense_changed);
			assert( next != NULL || (g__PSEUDO_NODES_PRESENT && node->right_color[i] != 0) );
			if (next != NULL) {
                assert( node->right_color[i] == 0 );
                int next_back = sense_changed ? COMPLEMENT(tail_leftmost_base) : tail_leftmost_base;
                counter_t *next_count = sense_changed ? next->right_count : next->left_count;
                Color *next_color = sense_changed ? next->right_color : next->left_color;
                assert( cnorm(next_count[next_back]) == cnorm(node->right_count[i]) );
                assert( next_color[next_back] == 0 );

                assert( !isNodeDead<SeqNode>(next) );
			}
		} else {
            assert( node->right_color[i] == 0 );
        }
	}
}

