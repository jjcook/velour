//
// This file is part of Velour and is distributed under the University of
// Illinois Open Source License.  See LICENSE.txt for details.
//

//
// stat_components.cpp
//
//   generates and reports statistics on the components in a graph
//

#include "types.h"

namespace statcomponents {

template<typename G, typename N>
struct NodeAnalysis {
    NodeAnalysis(G *graph, N *node, std::deque<N*>& worklist) :
        is_isolated(false),
        is_single_copy(false),
        is_tip(false),
        is_linear(false),
        is_branching(false),
        is_boundary(false),
        total_arcs(0),
        boundary_arcs(0)
    { runAnalysis(graph, node, worklist); }

    void runAnalysis(G *graph, N* node, std::deque<N*>& worklist);

	bool is_isolated;
	bool is_single_copy;
	bool is_tip;
	bool is_linear;
	bool is_branching;
	bool is_boundary;
	uintptr_t total_arcs;
	uintptr_t boundary_arcs;
};

bool maybeSingleCopy(KmerNode *node) {
    return true;
}

bool maybeSingleCopy(SeqNode *node) {
    return (isNodeWeight_None(node) || isNodeWeight_AllSingleCopy(node));
}

template<typename G, typename N>
inline void NodeAnalysis<G,N>::runAnalysis(G *graph, N *node, std::deque<N*>& worklist)
{
	unsigned left_count, right_count;
    node->getNeighborCounts(&left_count, &right_count);

    unsigned total_count = left_count + right_count;

    is_single_copy = maybeSingleCopy(node); // if skipping this, initialize to true in constructor

	if (total_count == 0) {
		is_isolated = true;
		return;
	}

	if (total_count == 1) {
		is_tip = true;
	} else if (left_count == 1 && right_count == 1) {
		is_linear = true;
	} else {
		is_branching = true;
	}

	total_arcs = left_count + right_count;

	// check right side
	for (int i = 0; i < 4; ++ i) {
		counter_t count = cnorm(node->right_count[i]);
		if (count != 0) {
			if (count > 1) {
				is_single_copy = false;
			}
            if (node->right_color[i] == 0) {
                N *next = graph->findNextNode(node, i, true);
				if (!isNodeMerged<N>(next)) {
					setNodeMerged<N>(next);
					worklist.push_front(next);
				}
			} else {
				is_boundary = true;
				++ boundary_arcs;
			}
		}
	}

	// check left side
	for (int i = 0; i < 4; ++ i) {
		counter_t count = cnorm(node->left_count[i]);
		if (count != 0) {
			if (count > 1) {
				is_single_copy = false;
			}
            if (node->left_color[i] == 0) {
                N *next = graph->findNextNode(node, i, false);
				if (!isNodeMerged<N>(next)) {
					setNodeMerged<N>(next);
					worklist.push_front(next);
				}
			} else {
				is_boundary = true;
				++ boundary_arcs;
			}
		}
	}
}

template<typename G>
struct ComponentsStatistics {
    ComponentsStatistics(G *graph) :
        total_nodes(graph->node_count),
        single_copy_nodes(0),
        isolated_nodes(0),
        //
        total_components(0),
        total_boundary_arcs(0),
        total_internal_arcs(0),
        //
        components(20, 1, 1),
        linear_components(20, 1, 1),
        single_copy_components(20, 1, 1),
        complete_components(20, 1, 1),
        //
        branching_components(21, 1, 0),
        tip_components(21, 1, 0),
        //
        internal_arcs_components(21, 1, 0),
        boundary_arcs_components(21, 1, 0),
        boundary_nodes_components(21, 1, 0)
    {}

    void print(FILE *); //const;

	uintptr_t total_nodes;
	uintptr_t single_copy_nodes;
	uintptr_t isolated_nodes;

	uintptr_t total_components;
	uintptr_t total_boundary_arcs;
	uintptr_t total_internal_arcs;

	Histogram components;                 // parameters: # nodes
	Histogram linear_components;          // parameters: # nodes
	Histogram single_copy_components;     // parameters: # nodes
	Histogram complete_components;        // parameters: # nodes

	Histogram branching_components;       // parameters: # branching nodes -- not # nodes
	Histogram tip_components;             // parameters: # tip nodes -- not # nodes

	Histogram internal_arcs_components;   // parameters: # internal arcs
	Histogram boundary_arcs_components;   // parameters: # boundary arcs
	Histogram boundary_nodes_components;  // parameters: # boundary nodes
};

template<typename G>
inline void ComponentsStatistics<G>::print(FILE *output) //const
{
	fputs("=== GRAPH COMPONENTS ===\n", output);
	fprintf(output, "  total_nodes         = %"PRIuPTR"\n", total_nodes);
	fprintf(output, "  single_copy_nodes   = %"PRIuPTR"\n", single_copy_nodes);
	fprintf(output, "  isolated_nodes      = %"PRIuPTR"\n", isolated_nodes);
	fprintf(output, "  total_components    = %"PRIuPTR"\n", total_components);
	fprintf(output, "  total_boundary_arcs = %"PRIuPTR"\n", total_boundary_arcs);
	fprintf(output, "  total_internal_arcs = %"PRIuPTR"\n", total_internal_arcs);
    fprintf(output, "  *** key:  name, %% of components, %% of nodes, average, max, weighted average -- histogram # nodes\n");
	                components.print("  all_components                 : ", total_components, total_nodes);
	         linear_components.print("  linear_components          : ", total_components, total_nodes);
	    single_copy_components.print("  single_copy_components     : ", total_components, total_nodes);
	       complete_components.print("  complete_components        : ", total_components, total_nodes);

	            tip_components.print("  tip_components             : ", total_components, total_nodes);
	      branching_components.print("  branching_components       : ", total_components, total_nodes);

	  internal_arcs_components.print("  internal_arcs_components   : ", total_components, total_internal_arcs+total_boundary_arcs);
	  boundary_arcs_components.print("  boundary_arcs_components   : ", total_components, total_internal_arcs+total_boundary_arcs);
	 boundary_nodes_components.print("  boundary_nodes_components  : ", total_components, total_nodes);

	fputs("=== END GRAPH COMPONENTS ===\n", output);
}

template<typename G, typename N>
struct lambda_functor {
    lambda_functor(G *graph, ComponentsStatistics<G> &cs) : graph(graph), cs(cs) {}
    void operator()(N *trav) {
        if (!isNodeMerged<N>(trav)) {
            bool is_tip_component = false;
            bool is_linear_component = true;
            bool is_single_copy_component = true;
            bool is_branching_component = false;
            bool is_complete_component = true;

            uintptr_t nodes = 0;
            uintptr_t tip_nodes = 0;
            uintptr_t branching_nodes = 0;
            uintptr_t boundary_nodes = 0;
            uintptr_t boundary_arcs = 0;
            uintptr_t internal_arcs = 0;

            std::deque<N*> worklist;

            setNodeMerged<N>(trav);
            worklist.push_front(trav);

            while (!worklist.empty()) {
                N *node = worklist.back();
                worklist.pop_back();

                NodeAnalysis<G,N> analysis(graph, node, worklist);

                cs.isolated_nodes += analysis.is_isolated;
                cs.single_copy_nodes += analysis.is_single_copy;

                is_tip_component |= analysis.is_tip;
                is_linear_component &= analysis.is_linear | analysis.is_tip; // ?
                is_single_copy_component &= analysis.is_single_copy;
                is_branching_component |= analysis.is_branching;
                is_complete_component &= !analysis.is_boundary;

                nodes += 1;
                tip_nodes += analysis.is_tip;
                branching_nodes += analysis.is_branching;
                boundary_nodes += analysis.is_boundary;
                boundary_arcs += analysis.boundary_arcs;
                internal_arcs += (analysis.total_arcs - analysis.boundary_arcs);
            }
            assert( nodes > 0 );

            ++ cs.total_components;
            cs.total_boundary_arcs += boundary_arcs;

            //assert( (internal_arcs & 0x1) == 0x0 ); // check that count is even before dividing by 2
            internal_arcs /= 2; // since we double count them (FIXME: except self-arcs!)
            cs.total_internal_arcs += internal_arcs;

            cs.components.insert(nodes);
            if (is_tip_component) cs.tip_components.insert(tip_nodes);
            if (is_linear_component) cs.linear_components.insert(nodes);
            if (is_single_copy_component) cs.single_copy_components.insert(nodes);
            if (is_branching_component) cs.branching_components.insert(branching_nodes);
            if (is_complete_component) cs.complete_components.insert(nodes);

            if (internal_arcs > 0) cs.internal_arcs_components.insert(internal_arcs);
            if (!is_complete_component) cs.boundary_arcs_components.insert(boundary_arcs);
            if (!is_complete_component) cs.boundary_nodes_components.insert(boundary_nodes);
        }
    }
  private:
    G *graph;
    ComponentsStatistics<G>& cs;
};

} // namespace: statcomponents

void kg_stat_components(KmerGraph *graph, FILE *output)
{
    statcomponents::ComponentsStatistics<KmerGraph> cs(graph);

    graph->resetFlags();

    statcomponents::lambda_functor<KmerGraph,KmerNode> _ft(graph, cs);
    kg_for_each(graph, _ft);

    cs.print(output);

    graph->resetFlags();
}

void sg_stat_components(SeqGraph *graph, FILE *output)
{
    statcomponents::ComponentsStatistics<SeqGraph> cs(graph);

    graph->resetFlags();

    statcomponents::lambda_functor<SeqGraph,SeqNode> _ft(graph, cs);
    sg_for_each(graph, _ft);

    cs.print(output);

    graph->resetFlags();
}

