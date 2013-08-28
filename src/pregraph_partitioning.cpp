//
// This file is part of Velour and is distributed under the University of
// Illinois Open Source License.  See LICENSE.txt for details.
//

//
// pregraph_partitioning.cpp
//

#include "types.h"

static uint64_t virtual_nodes = 0;
static uint64_t virtual_edges = 0;

//static vid_t current_vid;

//static std::vector< Sequence * > sequences;
static std::vector< vid_t > destinations;

static SeqGraph *graph = NULL;
static FILE * metis_file = NULL;

static void number_nodeside(SeqNode *node, vid_t vid, bool right_side)
{
    //assert( node->left_vid == 0 || node->right_vid == 0 );
    if (right_side) {
        if (node->right_vid == 0) {
            node->right_vid = vid;
            if (node->left_vid != 0 && node->left_vid != vid) { destinations.push_back( abs(node->left_vid) ); }
            if (node->right_side == 0 && node->left_vid == 0) { // ensure tip node sides are both the same vid
                number_nodeside(node, vid, GO_LEFT);
            } else {
                for (unsigned i=0; i < 4; ++i) {
                    if (node->right_count[i] != 0) {
                        bool sense_changed;
                        SeqNode *next = graph->findNextNode(node, i, GO_RIGHT, &sense_changed);
                        assert( next != NULL );
                        number_nodeside(next, vid, !(GO_RIGHT ^ sense_changed));
                    }
                }
            }
            if (node->left_side == 0 && node->left_vid == 0) { // ensure tip node sides are both the same vid
                number_nodeside(node, vid, GO_LEFT);
            }
            // TODO: merge unconcatenated nodes into a single vid
        }
        assert( node->right_vid == vid );
    } else {
        // left side
        if (node->left_vid == 0) {
            node->left_vid = vid;
            if (node->right_vid != 0 && node->right_vid != vid) { destinations.push_back( abs(node->right_vid) ); }
            if (node->left_side == 0 && node->right_vid == 0) { // ensure tip node sides are both the same vid
                number_nodeside(node, vid, GO_RIGHT);
            } else {
                for (unsigned i=0; i < 4; ++i) {
                    if (node->left_count[i] != 0) {
                        bool sense_changed;
                        SeqNode *next = graph->findNextNode(node, i, GO_LEFT, &sense_changed);
                        assert( next != NULL );
                        number_nodeside(next, vid, !(GO_LEFT ^ sense_changed));
                    }
                }
            }
            if (node->right_side == 0 && node->right_vid == 0) { // ensure tip node sides are both the same vid
                number_nodeside(node, vid, GO_RIGHT);
            }
            // TODO: merge unconcatenated nodes into a single vid
        }
        assert( node->left_vid == vid );
    }
    //assert( node->left_vid == vid || node->right_vid == vid );
}

static void memoize_nodeside(SeqNode *node, vid_t vid, bool right_side)
{
    assert( vid > 0 );
    if (right_side) {
        if (node->right_vid == vid) {
            node->right_vid *= -1;
            if (node->left_vid == vid) { // an unconcatenated linear path OR a self-loop
                //sequences.push_back(&node->sequence);
                memoize_nodeside(node, vid, GO_LEFT);
            } else if (abs(node->left_vid) != vid) {
                //sequences.push_back(&node->sequence);
                destinations.push_back( abs(node->left_vid) );
            }
            for (unsigned i=0; i < 4; ++i) {
                if (node->right_count[i] != 0) {
                    bool sense_changed;
                    SeqNode *next = graph->findNextNode(node, i, GO_RIGHT, &sense_changed);
                    assert( next != NULL );
                    memoize_nodeside(next, vid, !(GO_RIGHT ^ sense_changed));
                }
            }
        }
    } else {
        // left side
        if (node->left_vid == vid) {
            node->left_vid *= -1;
            if (node->right_vid == vid) { // an unconcatenated linear path OR a self-loop
                // // sequences.push_back(&node->sequence); // commented to avoid duplicate entry
                memoize_nodeside(node, vid, GO_RIGHT);
            } else if (abs(node->right_vid) != vid) {
                //sequences.push_back(&node->sequence);
                destinations.push_back( abs(node->right_vid) );
            }
            for (unsigned i=0; i < 4; ++i) {
                if (node->left_count[i] != 0) {
                    bool sense_changed;
                    SeqNode *next = graph->findNextNode(node, i, GO_LEFT, &sense_changed);
                    assert( next != NULL );
                    memoize_nodeside(next, vid, !(GO_LEFT ^ sense_changed));
                }
            }
        }
    }
}

static void emit_nodeside(vid_t vid)
{
    //
    // remove duplicates 
    //
    //std::sort( sequences.begin(), sequences.end() );
    //sequences.erase( std::unique( sequences.begin(), sequences.end() ), sequences.end() );
    
    std::sort( destinations.begin(), destinations.end() );
    destinations.erase( std::unique( destinations.begin(), destinations.end() ), destinations.end() );

    //
    // emit sequences and metis file entry
    //
    /*current_vid = vid;
    struct lambda1 {
        static void emit_sequence(Sequence *seq)
        {
            fwrite(&current_vid, sizeof(current_vid), 1, sequence_file);
            seq->Save_BinaryFile(sequence_file);
        }
    };
    std::for_each( sequences.begin(), sequences.end(), lambda1::emit_sequence );*/

    struct lambda2 {
        static void emit_metis_edge(vid_t dest_vid)
        {
            fprintf(metis_file, "%"PRId64" ", dest_vid);
        }
    };
    std::for_each(  destinations.begin(), destinations.end(), lambda2::emit_metis_edge);
    
    if (destinations.size() == 0) {
      // FIXME ?? node/cluster is self contained, e.g. a single node, self-loop single node, etc
      //   assign to random, hopefully okay if random is also a self contained blob
      fprintf(metis_file, "%"PRId64" ", (random() % virtual_nodes)+1);
    }

    fprintf(metis_file, "\n");

    //
    // clean up
    //
    //sequences.clear();
    destinations.clear();
}

void pregraph_partitioning(SeqGraph *sgraph, char *metis_filename)
{
    graph = sgraph;

// VERIFY
// TODO: check that color fields are all zeros!!!
// end VERIFY
   
    // zero the vid fields TODO OPT may not be necessary
    struct lambda_reset {
       static void reset_vids(SeqNode *node)
       {
          node->left_vid = 0; 
          node->right_vid = 0;
       }
    };
    sg_for_each(graph, lambda_reset::reset_vids);

    //
    // compute the virtual hypergraph by numbering the nodes
    //
    struct lambda1 {
        static void number_node(SeqNode *node) 
        {
            if (node->right_vid == 0) {
                vid_t vid = ++ virtual_nodes;
                number_nodeside(node, vid, GO_RIGHT);
   
                // eliminate duplicates 
                std::sort( destinations.begin(), destinations.end() );
                destinations.erase( std::unique( destinations.begin(), destinations.end() ), destinations.end() );
                virtual_edges += destinations.size();
                destinations.clear();
            }
            if (node->left_vid == 0) {
                vid_t vid = ++ virtual_nodes;
                number_nodeside(node, vid, GO_LEFT);

                // eliminate duplicates 
                std::sort( destinations.begin(), destinations.end() );
                destinations.erase( std::unique( destinations.begin(), destinations.end() ), destinations.end() );
                virtual_edges += destinations.size();
                destinations.clear();
            }
        }
    };
    sg_for_each(graph, lambda1::number_node);

// VERIFY
// TODO: check that node numbers are all non-zero
// end VERIFY

    //
    // create file for metis input hypergraph
    //
    metis_file = fopen(metis_filename, "w");
    if (metis_file == NULL) {
        printf("Could not create file %s, exiting...\n", metis_filename);
        exit(EXIT_FAILURE);
    }
    fprintf(metis_file, "%"PRIu64" %"PRIu64"\n", virtual_nodes, virtual_edges); // metis file header

    //
    // generate and emit the metis hypergraph entries
    //
    struct lambda2 {
        static void emit_node(SeqNode *node)
        {
            assert( node->right_vid != 0 );
            assert( node->left_vid != 0 );
            if (node->right_vid > 0) {
                vid_t vid = node->right_vid;
                memoize_nodeside(node, vid, GO_RIGHT);
                emit_nodeside(vid);
            }
            if (node->left_vid > 0) {
                vid_t vid = node->left_vid;
                memoize_nodeside(node, vid, GO_LEFT);
                emit_nodeside(vid);
            }
        }
    };
    sg_for_each(graph, lambda2::emit_node);

    //
    // clean up
    //
    fclose(metis_file);
    metis_file = NULL;

    graph = NULL;
}
