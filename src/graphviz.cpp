//
// This file is part of Velour and is distributed under the University of
// Illinois Open Source License.  See LICENSE.txt for details.
//

//
// graphviz.cpp
//
//   graphviz output of kmer & sequence graph
//

#include "types.h"

static void emitNode(KmerNode *node, FILE *file)
{
    bool has_pseudo_arc = false;

    // node identifier
    fprint_kmer(node->kmer, g__FULLKMER_LENGTH, file);
    fprintf(file, " [label=\"");

    // left connections
    fprintf(file, "{");
    for (unsigned i=0; i < 4; ++i) {
        //counter_t count = cnorm(node->left_count[i]);
        counter_t count = node->left_count[i];
        fprintf(file, "<L%d>", i);
        if (count != 0) fprintf(file, "%d", count);
        if (i != 3) fprintf(file,"|");
    }
    fprintf(file, "} |");

    // left colors
    fprintf(file, "{");
    for (unsigned i=0; i < 4; ++i) {
        Color color = node->left_color[i];
        if (color != 0) { fprintf(file, "%d", color); has_pseudo_arc = true; }
        if (i != 3) fprintf(file,"|");
    }
    fprintf(file, "} |");

    // middle fields: kmer, rc_kmer, claim_tid, flags
    fprintf(file, "{");
    fprint_kmer(node->kmer, g__FULLKMER_LENGTH, file);
    fprintf(file, "|");
    fprint_kmer(reverseComplement(node->kmer, g__FULLKMER_LENGTH), g__FULLKMER_LENGTH, file);
    fprintf(file, "|");
    fprintf(file, "{tid=%u|flags=0x%01x}", node->claim_tid, node->flags);
    fprintf(file, "} |");

    // right colors
    fprintf(file, "{");
    for (unsigned i=0; i < 4; ++i) {
        Color color = node->right_color[i];
        if (color != 0) { fprintf(file, "%d", color); has_pseudo_arc = true; }
        if (i != 3) fprintf(file,"|");
    }
    fprintf(file, "} |");

    // right connections
    fprintf(file, "{");
    for (unsigned i=0; i < 4; ++i) {
        //counter_t count = cnorm(node->right_count[i]);
        counter_t count = node->right_count[i];
        fprintf(file, "<R%d>", i);
        if (count != 0) fprintf(file, "%d", count);
        if (i != 3) fprintf(file,"|");
    }
    fprintf(file, "}");

    fprintf(file, "\"];\n");

    if (has_pseudo_arc) {
        // node identifier
        fprint_kmer(node->kmer, g__FULLKMER_LENGTH, file);
        fprintf(file, " [style=filled,color=tomato];\n");
    }
}

static void emitNode(SeqNode *node, FILE *file)
{
    bool has_pseudo_arc = false;

    Kmer head_kmer = node->sequence.GetHeadKmer(g__FULLKMER_LENGTH);
    Kmer tail_kmer = node->sequence.GetTailKmer(g__FULLKMER_LENGTH);

    // node identifier
    fprint_kmer(head_kmer, g__FULLKMER_LENGTH, file);
    fprintf(file, " [label=\"");

    // left connections
    fprintf(file, "{");
    for (unsigned i=0; i < 4; ++i) {
        //counter_t count = cnorm(node->left_count[i]);
        counter_t count = node->left_count[i];
        fprintf(file, "<L%d>", i);
        if (count != 0) fprintf(file, "%d", count);
        if (i != 3) fprintf(file,"|");
    }
    fprintf(file, "} |");

    // left colors
    fprintf(file, "{");
    for (unsigned i=0; i < 4; ++i) {
        Color color = node->left_color[i];
        if (color != 0) { fprintf(file, "%d", color); has_pseudo_arc = true; }
        if (i != 3) fprintf(file,"|");
    }
    fprintf(file, "} |");

    // middle fields: head_kmer, tail_kmer, length-min-max, flags
    fprintf(file, "{");
    fprint_kmer(head_kmer, g__FULLKMER_LENGTH, file);
    fprintf(file, "|");
    fprint_kmer(tail_kmer, g__FULLKMER_LENGTH, file);
    fprintf(file, "|");
    fprintf(file, "{length=%u|flags=0x%01x}", node->sequence.get_length(), node->flags);
    fprintf(file, "|");
    fprintf(file, "{slice=%u|sflags=0x%01x|tid=%u}", node->slice_color, node->slice_flags, node->claim_tid);
    fprintf(file, "} |");

    // right colors
    fprintf(file, "{");
    for (unsigned i=0; i < 4; ++i) {
        Color color = node->right_color[i];
        if (color != 0) { fprintf(file, "%d", color); has_pseudo_arc = true; }
        if (i != 3) fprintf(file,"|");
    }
    fprintf(file, "} |");

    // right connections
    fprintf(file, "{");
    for (unsigned i=0; i < 4; ++i) {
        //counter_t count = cnorm(node->right_count[i]);
        counter_t count = node->right_count[i];
        fprintf(file, "<R%d>", i);
        if (count != 0) fprintf(file, "%d", count);
        if (i != 3) fprintf(file,"|");
    }
    fprintf(file, "}");

    fprintf(file, "\"];\n");

    // full sequence
    fprint_kmer(head_kmer, g__FULLKMER_LENGTH, file);
    fprintf(file, " [comment=\"");
    node->sequence.PrintToFile(file);
    fprintf(file, "\"];\n");

    if (has_pseudo_arc) {
        // node identifier
        fprint_kmer(head_kmer, g__FULLKMER_LENGTH, file);
        fprintf(file, " [style=filled,color=tomato];\n");
    }
}

static void emitNodeConnections(KmerGraph *graph, KmerNode *node, FILE *file)
{
	for (int i = 0 ; i < 4 ; ++ i) {
	   	// check on the left side
		if (node->left_count[i] != 0 && node->left_color[i] == 0) {
			KmerNode *next = graph->findNextNode(node, i, false);

            if (node <= next) {
                fprint_kmer(node->kmer, g__FULLKMER_LENGTH, file);
                fprintf(file, ":L%d -- ", i);
                fprint_kmer(next->kmer, g__FULLKMER_LENGTH, file);
                if (node->left_count[i] < 0) {
                    fprintf(file, ":L%u ", static_cast<unsigned>(COMPLEMENT(KMER_GET_TAIL_BASE(node->kmer, g__FULLKMER_LENGTH))));
                } else {
                    fprintf(file, ":R%u ", static_cast<unsigned>(KMER_GET_TAIL_BASE(node->kmer, g__FULLKMER_LENGTH)));
                }
                fprintf(file, "[weight=%d] ;\n", cnorm(node->left_count[i]));
            }
		}

		// check on the right side
		if (node->right_count[i] != 0 && node->right_color[i] == 0) {
			KmerNode *next = graph->findNextNode(node, i, true);

            if (node <= next) {
                fprint_kmer(node->kmer, g__FULLKMER_LENGTH, file);
                fprintf(file, ":R%d -- ", i);
                fprint_kmer(next->kmer, g__FULLKMER_LENGTH, file);
                if (node->right_count[i] < 0) {
                    fprintf(file, ":R%u ", static_cast<unsigned>(COMPLEMENT(KMER_GET_HEAD_BASE(node->kmer, g__FULLKMER_LENGTH))));
                } else {
                    fprintf(file, ":L%u ", static_cast<unsigned>(KMER_GET_HEAD_BASE(node->kmer, g__FULLKMER_LENGTH)));
                }
                fprintf(file, "[weight=%d] ;\n", cnorm(node->right_count[i]));
            }
        }
	}
}

static void emitNodeConnections(SeqGraph *graph, SeqNode *node, FILE *file)
{
    Kmer head_kmer = node->sequence.GetHeadKmer(g__FULLKMER_LENGTH);

	Nucleotide head_rightmost_base = node->sequence.GetHeadKmerRightmostBase(g__FULLKMER_LENGTH);
	Nucleotide tail_leftmost_base = node->sequence.GetTailKmerLeftmostBase(g__FULLKMER_LENGTH);
	for (int i = 0 ; i < 4 ; ++ i) {
	   	// check on the left side
		if (node->left_count[i] != 0 && node->left_color[i] == 0) {
            bool sense_changed;
			SeqNode *next = graph->findNextNode(node, i, false, &sense_changed);

            if (node <= next) {
                int next_back = sense_changed ? COMPLEMENT(head_rightmost_base) : head_rightmost_base;
                Kmer next_head_kmer = next->sequence.GetHeadKmer(g__FULLKMER_LENGTH);

                fprint_kmer(head_kmer, g__FULLKMER_LENGTH, file);
                fprintf(file, ":L%d -- ", i);
                fprint_kmer(next_head_kmer, g__FULLKMER_LENGTH, file);
                if (sense_changed) {
                    fprintf(file, ":L%d ", next_back);
                } else {
                    fprintf(file, ":R%d ", next_back);
                }
                fprintf(file, "[weight=%d] ;\n", cnorm(node->left_count[i]));
            }
		}

		// check on the right side
		if (node->right_count[i] != 0 && node->right_color[i] == 0) {
            bool sense_changed;
			SeqNode *next = graph->findNextNode(node, i, true, &sense_changed);

            if (node <= next) {
                int next_back = sense_changed ? COMPLEMENT(tail_leftmost_base) : tail_leftmost_base;
                Kmer next_head_kmer = next->sequence.GetHeadKmer(g__FULLKMER_LENGTH);

                fprint_kmer(head_kmer, g__FULLKMER_LENGTH, file);
                fprintf(file, ":R%d -- ", i);
                fprint_kmer(next_head_kmer, g__FULLKMER_LENGTH, file);
                if (sense_changed) {
                    fprintf(file, ":R%d ", next_back);
                } else {
                    fprintf(file, ":L%d ", next_back);
                }
                fprintf(file, "[weight=%d] ;\n", cnorm(node->right_count[i]));
            }
        }
	}
}

namespace {

template<typename GraphType, typename NodeType>
struct lambda_graphviz {
    lambda_graphviz(GraphType *graph, FILE *file) : graph(graph), file(file) {}
    void operator()(NodeType *node) { emitNode(node, file); emitNodeConnections(graph, node, file); }
  private:
    GraphType *graph;
    FILE *file;
};

} // namespace: anonymous

void emit_scoop_graphviz(SeqGraph *graph, SeqNode *root, intptr_t nodes, char *filename)
{
    FILE *file = fopen(filename, "w");

    fprintf(file, "graph scoopedSequenceGraph {\n");
    fprintf(file, "node [shape=record];\n");
    fprintf(file, "graph [fontsize=8,rankdir=TB];\n");

    lambda_graphviz<SeqGraph, SeqNode> ft(graph, file);

    // scoop
    if (root->claim_tid == 0) {
        std::deque<SeqNode *> worklist;
        std::deque<SeqNode *> nodelist;

        root->claim_tid = 1;
        worklist.push_front(root);
        nodelist.push_front(root);
        ft(root);
        -- nodes;
        while(!worklist.empty() && nodes > 0) {
            SeqNode *node = worklist.back();
            worklist.pop_back();

            // check right side
            counter_t *counters = node->right_count;
            Color *colors = node->right_color;
            for (int i = 0; i < 4; ++ i) {
                if (counters[i] != 0) {
                    if (colors[i] == 0) {
                        SeqNode *next = graph->findNextNode(node, i, GO_RIGHT);
                        assert( !isNodeDead<SeqNode>(next) );
                        if (next->claim_tid == 0) {
                            next->claim_tid = 1;
                            worklist.push_front(next);
                            nodelist.push_front(next);
                            ft(next);
                            -- nodes;
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
                        SeqNode *next = graph->findNextNode(node, i, GO_LEFT);
                        assert( !isNodeDead<SeqNode>(next) );
                        if (next->claim_tid == 0) {
                            next->claim_tid = 1;
                            worklist.push_front(next);
                            nodelist.push_front(next);
                            ft(next);
                            -- nodes;
                        }
                    }
                }
            }
        }

        // reset claim_tid
        while(!nodelist.empty()) {
            SeqNode *node = nodelist.back();
            nodelist.pop_back();
            node->claim_tid = 0;
        }
    }

    // end of graph
    fprintf(file, "}\n");

    fclose(file);
}

void emit_graphviz(KmerGraph *graph, char *filename)
{
    FILE *file = fopen(filename, "w");
    
    fprintf(file, "graph kmerGraph {\n");
    fprintf(file, "node [shape=record];\n");
    fprintf(file, "graph [fontsize=8,rankdir=TB];\n");

    lambda_graphviz<KmerGraph, KmerNode> ft(graph, file);
    kg_for_each(graph, ft);

    // end of graph
    fprintf(file, "}\n");

    fclose(file);
}

void emit_graphviz(SeqGraph *graph, char *filename)
{
    FILE *file = fopen(filename, "w");
    
    fprintf(file, "graph sequenceGraph {\n");
    fprintf(file, "node [shape=record];\n");
    fprintf(file, "graph [fontsize=8,rankdir=TB];\n");

    lambda_graphviz<SeqGraph, SeqNode> ft(graph, file);
    sg_for_each(graph, ft);

    // end of graph
    fprintf(file, "}\n");

    fclose(file);
}
