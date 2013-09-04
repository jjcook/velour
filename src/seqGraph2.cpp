//
// This file is part of Velour and is distributed under the University of
// Illinois Open Source License.  See LICENSE.txt for details.
//

//
// seqGraph.cpp
//

#include "types.h"
#include "seqGraph.h"

SeqGraph::SeqGraph(uintptr_t _buckets) :
    buckets(round_up_power_of_2(_buckets)), hashmask(buckets-1),
    head_table(static_cast<SeqNode**>(calloc(buckets, sizeof(SeqNode *)))),
    tail_table(static_cast<SeqNode**>(calloc(buckets, sizeof(SeqNode *)))), node_count(0)
#ifdef VERIFY
#  ifdef VELOUR_TBB
    , iterator_lock() // defaults to zero
#  else
    , iterator_lock(0)
#  endif
#endif
    {}

SeqGraph::~SeqGraph() { free(head_table); free(tail_table); }

namespace {
    struct lambda_list_node {
        lambda_list_node() {}
        void operator()(SeqNode *node) {
            printf("node: %p  claim_tid: %u\n", node, node->claim_tid);
        }
    };
}

void gdb_list_seqgraph(SeqGraph *sgraph)
{
    lambda_list_node l;
    sg_for_each(sgraph, l);
    fflush(stdout);
}

void SeqGraph::insertNodeAndUpdateColors(SeqNode *node)
{
    // first, insert the node
    insertNode(node);

    // update (clear) color on edges (and back color edge)
	// XXX NOTE: optimized this w.r.t the redistribution / recombination scheme, e.g. don't check if color is larger than X
	Nucleotide head_rightmost_base = node->sequence.GetHeadKmerRightmostBase(g__FULLKMER_LENGTH);
	Nucleotide tail_leftmost_base = node->sequence.GetTailKmerLeftmostBase(g__FULLKMER_LENGTH);
	for (int i = 0 ; i < 4 ; ++ i) {
	   	// check on the left side
		if (node->left_color[i] > 0 && node->left_color[i] <= g__PARTITION_INDEX) { // NOTE: checks for update w.r.t. partition index
            assert( node->left_count[i] != 0 );
            bool sense_changed;
			SeqNode *next = this->findNextNode(node, i, GO_LEFT, &sense_changed);
			if (next != NULL) {
                int next_back = sense_changed ? COMPLEMENT(head_rightmost_base) : head_rightmost_base;
                counter_t *next_count = sense_changed ? next->left_count : next->right_count;
                Color *next_color = sense_changed ? next->left_color : next->right_color;
                assert( cnorm(next_count[next_back]) == cnorm(node->left_count[i]) );
                assert( next_color[next_back] != 0 );
                node->left_color[i] = 0;
                next_color[next_back] = 0;
			}
		}

		// check on the right side
        if (node->right_color[i] > 0 && node->right_color[i] <= g__PARTITION_INDEX) { // NOTE: checks for update w.r.t. partition index
            assert( node->right_count[i] != 0 );
            bool sense_changed;
			SeqNode *next = this->findNextNode(node, i, GO_RIGHT, &sense_changed);
            if (next != NULL) {
                int next_back = sense_changed ? COMPLEMENT(tail_leftmost_base) : tail_leftmost_base;
                counter_t *next_count = sense_changed ? next->right_count : next->left_count;
                Color *next_color = sense_changed ? next->right_color : next->left_color;
                assert( cnorm(next_count[next_back]) == cnorm(node->right_count[i]) );
                assert( next_color[next_back] != 0 );
                node->right_color[i] = 0;
                next_color[next_back] = 0;
            }
        }
	}
}

// find the previous SeqNode with the matching head kmer, or return NULL
static inline SeqNode * findPreviousNodeHead(Kmer head_kmer, SeqNode **head_bucket)
{
	SeqNode *prev = NULL;

	for (SeqNode *trav = *head_bucket ; trav != NULL ; prev = trav, trav = trav->head_next) {
        Kmer trav_head_kmer = trav->sequence.GetHeadKmer(g__FULLKMER_LENGTH);
		if (trav_head_kmer == head_kmer) { // not canonical kmers
			return prev;
		}
	}
	return NULL;
}

static inline SeqNode * findPreviousNodeTail(Kmer tail_kmer, SeqNode **tail_bucket)
{
	SeqNode *prev = NULL;

	for (SeqNode *trav = *tail_bucket ; trav != NULL ; prev = trav, trav = trav->tail_next) {
        Kmer trav_tail_kmer = trav->sequence.GetTailKmer(g__FULLKMER_LENGTH);
		if (trav_tail_kmer == tail_kmer) { // not canonical kmers
			return prev;
		}
	}
	return NULL;
}

void SeqGraph::removeNodeHead(SeqNode *node)
{
	Kmer head_kmer = node->sequence.GetHeadKmer(g__FULLKMER_LENGTH);
	Kmer head_canon_kmer = canonicalKmer(head_kmer, g__FULLKMER_LENGTH);
    hash_index_t hash_index = hash(head_canon_kmer);
    SeqNode **head_bucket = &head_table[hash_index];

    assert(*head_bucket != NULL);
    if (*head_bucket == node) {
		*head_bucket = node->head_next;
    } else {
        SeqNode *prev_node = findPreviousNodeHead(head_kmer, head_bucket);  // not the canonical kmer
        assert( prev_node != NULL && "missing or duplicate head kmer?" );
		prev_node->head_next = node->head_next;
    }
    node->head_next = NULL;
}

void SeqGraph::removeNodeTail(SeqNode *node)
{
	Kmer tail_kmer = node->sequence.GetTailKmer(g__FULLKMER_LENGTH);
	Kmer tail_canon_kmer = canonicalKmer(tail_kmer, g__FULLKMER_LENGTH);
    hash_index_t hash_index = hash(tail_canon_kmer);
    SeqNode **tail_bucket = &tail_table[hash_index];

    assert(*tail_bucket != NULL);
    if (*tail_bucket == node) {
		*tail_bucket = node->tail_next;
    } else {
        SeqNode *prev_node = findPreviousNodeTail(tail_kmer, tail_bucket);  // not the canonical kmer
        assert( prev_node != NULL && "missing or duplicate tail kmer?" );
		prev_node->tail_next = node->tail_next;
    }
    node->tail_next = NULL;
}

#ifdef VELOUR_TBB
void SeqGraph::atomic_removeNodeHead(SeqNode *node)
{
	Kmer head_kmer = node->sequence.GetHeadKmer(g__FULLKMER_LENGTH);
	Kmer head_canon_kmer = canonicalKmer(head_kmer, g__FULLKMER_LENGTH);
    hash_index_t hash_index = hash(head_canon_kmer);
    SeqNode **head_bucket = &head_table[hash_index];

    assert(*head_bucket != NULL);
    while(1) {
        SeqNode *old_ptr;
        while ((old_ptr = ATOMIC_LOAD(*head_bucket)) == reinterpret_cast<SeqNode*>(1)) {} // spin waiting for unlock
        if (ATOMIZE(*head_bucket).compare_and_swap(reinterpret_cast<SeqNode*>(1), old_ptr) == old_ptr) { // lock bucket with poison value
            if (old_ptr == node) { // simple case: removing the first node, so just update the bucket pointer
                ATOMIC_STORE(*head_bucket) = node->head_next; // unlock
                break;
            } else {
                SeqNode *prev = NULL;
                for (SeqNode *trav = old_ptr ; trav != NULL ; prev = trav, trav = trav->head_next) {
                    Kmer trav_head_kmer = trav->sequence.GetHeadKmer(g__FULLKMER_LENGTH);
                    if (trav_head_kmer == head_kmer) { // not canonical kmers
                        prev->head_next = node->head_next;
                        break;
                    }
                }
                assert( prev != NULL && "missing or duplicate head kmer?" );
                ATOMIC_STORE(*head_bucket) = old_ptr; // unlock
                break;
            }
		}
	}
    node->head_next = NULL;
}

void SeqGraph::atomic_removeNodeTail(SeqNode *node)
{
	Kmer tail_kmer = node->sequence.GetTailKmer(g__FULLKMER_LENGTH);
	Kmer tail_canon_kmer = canonicalKmer(tail_kmer, g__FULLKMER_LENGTH);
    hash_index_t hash_index = hash(tail_canon_kmer);
    SeqNode **tail_bucket = &tail_table[hash_index];

    assert(*tail_bucket != NULL);
    while(1) {
        SeqNode *old_ptr;
        while ((old_ptr = ATOMIC_LOAD(*tail_bucket)) == reinterpret_cast<SeqNode*>(1)) {} // spin waiting for unlock
        if (ATOMIZE(*tail_bucket).compare_and_swap(reinterpret_cast<SeqNode*>(1), old_ptr) == old_ptr) { // lock bucket with poison value
            if (old_ptr == node) { // simple case: removing the first node, so just update the bucket pointer
                ATOMIC_STORE(*tail_bucket) = node->tail_next; // unlock
                break;
            } else {
                SeqNode *prev = NULL;
                for (SeqNode *trav = old_ptr ; trav != NULL ; prev = trav, trav = trav->tail_next) {
                    Kmer trav_tail_kmer = trav->sequence.GetTailKmer(g__FULLKMER_LENGTH);
                    if (trav_tail_kmer == tail_kmer) { // not canonical kmers
                        prev->tail_next = node->tail_next;
                        break;
                    }
                }
                assert( prev != NULL && "missing or duplicate tail kmer?" );
                ATOMIC_STORE(*tail_bucket) = old_ptr; // unlock
                break;
            }
		}
	}
    node->tail_next = NULL;
}
#endif // VELOUR_TBB

// find the SeqNode with the matching head/tail kmer, or return NULL
// NOTE: can't just compare canonical kmers since need to detect sense changed of single kmer nodes
SeqNode* SeqGraph::findNode(const Kmer kmer, bool moving_right, bool* sense_changed) const
{
    Kmer rc_kmer = reverseComplement(kmer, g__FULLKMER_LENGTH);
    Kmer canon_kmer = canonicalKmer(kmer, g__FULLKMER_LENGTH); // OPT: don't do RC twice
    hash_index_t hash_index = hash(canon_kmer);

	for (SeqNode *trav = this->head_table[hash_index] ; trav != NULL ; trav = trav->head_next) {
        if (isNodeDead<SeqNode>(trav)) continue;
		Kmer head_kmer = trav->sequence.GetHeadKmer(g__FULLKMER_LENGTH);
		if (moving_right) {
            if (head_kmer == kmer) {
                if (sense_changed != NULL) *sense_changed = 0;
                return trav;
            }
		} else {
            if (head_kmer == rc_kmer) {
                if (sense_changed != NULL) *sense_changed = 1;
                return trav;
            }
        }
	}
	
	for (SeqNode *trav = this->tail_table[hash_index] ; trav != NULL ; trav = trav->tail_next) {
        if (isNodeDead<SeqNode>(trav)) continue;
		Kmer tail_kmer = trav->sequence.GetTailKmer(g__FULLKMER_LENGTH);
		if (moving_right) {
            if (tail_kmer == rc_kmer) {
                if (sense_changed != NULL) *sense_changed = 1;
                return trav;
            }
		} else {
            if (tail_kmer == kmer) {
                if (sense_changed != NULL) *sense_changed = 0;
                return trav;
            }
		}
	}

	return NULL;
}

#ifdef VELOUR_TBB
SeqNode* SeqGraph::atomic_findNode(const Kmer kmer, bool moving_right, bool* sense_changed) const
{
    Kmer rc_kmer = reverseComplement(kmer, g__FULLKMER_LENGTH);
    Kmer canon_kmer = canonicalKmer(kmer, g__FULLKMER_LENGTH); // OPT: don't do RC twice
    hash_index_t hash_index = hash(canon_kmer);

    SeqNode **head_bucket = &this->head_table[hash_index];
    SeqNode *retval = NULL;
    bool found = false;
    if (ATOMIC_LOAD(*head_bucket) != NULL) {
        while(1) {
            SeqNode *old_ptr;
            while ((old_ptr = ATOMIC_LOAD(*head_bucket)) == reinterpret_cast<SeqNode*>(1)) {} // spin waiting for unlock
            if (ATOMIZE(*head_bucket).compare_and_swap(reinterpret_cast<SeqNode*>(1), old_ptr) == old_ptr) { // lock bucket with poison value
                for (SeqNode *trav = old_ptr ; trav != NULL ; trav = trav->head_next) {
                    if (isNodeDead<SeqNode>(trav)) continue;
                    Kmer head_kmer = trav->sequence.GetHeadKmer(g__FULLKMER_LENGTH);
                    if (moving_right) {
                        if (head_kmer == kmer) {
                            if (sense_changed != NULL) *sense_changed = 0;
                            retval = trav;
                            found = true;
                            break;
                        }
                    } else {
                        if (head_kmer == rc_kmer) {
                            if (sense_changed != NULL) *sense_changed = 1;
                            retval = trav;
                            found = true;
                            break;
                        }
                    }
                }

                // last, restore old pointer to unlock bucket
                ATOMIC_STORE(*head_bucket) = old_ptr; // unlock
                break;
            }
        }
    }
    if (found) { return retval; }
    
    SeqNode **tail_bucket = &this->tail_table[hash_index];
    if (ATOMIC_LOAD(*tail_bucket) != NULL) {
        while(1) {
            SeqNode *old_ptr;
            while ((old_ptr = ATOMIC_LOAD(*tail_bucket)) == reinterpret_cast<SeqNode*>(1)) {} // spin waiting for unlock
            if (ATOMIZE(*tail_bucket).compare_and_swap(reinterpret_cast<SeqNode*>(1), old_ptr) == old_ptr) { // lock bucket with poison value
                for (SeqNode *trav = old_ptr ; trav != NULL ; trav = trav->tail_next) {
                    if (isNodeDead<SeqNode>(trav)) continue;
                    Kmer tail_kmer = trav->sequence.GetTailKmer(g__FULLKMER_LENGTH);
                    if (moving_right) {
                        if (tail_kmer == rc_kmer) {
                            if (sense_changed != NULL) *sense_changed = 1;
                            retval = trav;
                            found = true;
                            break;
                        }
                    } else {
                        if (tail_kmer == kmer) {
                            if (sense_changed != NULL) *sense_changed = 0;
                            retval = trav;
                            found = true;
                            break;
                        }
                    }
                }

                // last, restore old pointer to unlock bucket
                ATOMIC_STORE(*tail_bucket) = old_ptr; // unlock
                break;
            }
        }
    }
    if (found) {
        assert(retval != NULL);
        return retval;
    } else {
        return NULL;
    }
}
#endif // VELOUR_TBB
     
static Kmer computeNextKmer(const SeqNode *node, Nucleotide next_base, bool moving_right)
{
    int double_kmer_length = g__FULLKMER_LENGTH << 1;
#ifdef LARGE_KMERS
    Kmer mask;
    mask.createMask(double_kmer_length);
#else
    Kmer mask = (Kmer(1) << double_kmer_length) - 1;
#endif

    Kmer next_kmer; 
    if (moving_right) {
        Kmer kmer = node->sequence.GetTailKmer(g__FULLKMER_LENGTH);
        next_kmer = KMER_APPEND(kmer, next_base, double_kmer_length);
    } else {
        Kmer kmer = node->sequence.GetHeadKmer(g__FULLKMER_LENGTH);
        next_kmer = KMER_PREPEND(kmer, next_base, double_kmer_length, mask);
    }
    return next_kmer;
}

SeqNode * SeqGraph::findNextNode(const SeqNode *node, Nucleotide next_base, bool moving_right, bool *sense_changed_ptr) const
{
    Kmer next_kmer = computeNextKmer(node, next_base, moving_right);
    SeqNode *next = this->findNode(next_kmer, moving_right, sense_changed_ptr);

#ifdef VERIFY
    /* TODO if (!g__LOADING_NODES) {
        if (next != NULL) {
            assert( (moving_right ? node->right_color : node->left_color)[next_base] == 0 );
        } else {
            assert( g__PSEUDO_NODES_PRESENT );
            assert( (moving_right ? node->right_color : node->left_color)[next_base] != 0 );
            //assert( (moving_right ? node->right_color : node->left_color)[next_base] != currentPartitionIndex );
        }
    } */
#endif // VERIFY
    return next;
}

#ifdef VELOUR_TBB
SeqNode * SeqGraph::atomic_findNextNode(const SeqNode *node, Nucleotide next_base, bool moving_right, bool *sense_changed_ptr) const
{
    Kmer next_kmer = computeNextKmer(node, next_base, moving_right);
    SeqNode *next = this->atomic_findNode(next_kmer, moving_right, sense_changed_ptr);
    // VERIFY?
    return next;
}
#endif // VELOUR_TBB


void SeqGraph::bulkMoveAllNodes(SeqGraph *target)
{
#ifdef VELOUR_TBB
    std::deque<SeqNode*> node_list;
#endif
    for (uintptr_t i = 0 ; i < this->buckets ; ++ i) {
        SeqNode *node;
        while ((node = this->head_table[i]) != NULL) {
            this->head_table[i] = node->head_next;
            node->head_next = NULL;
            this->removeNodeTail(node);
            -- this->node_count;
            assert( !isNodeDead<SeqNode>(node) );
            resetNodeMergingFlags<SeqNode>(node);
            // resetNodeSliceFlags(node); node->slice_color = 0; XXX enable these???
#ifdef VELOUR_TBB
            // first, we set the slice color to an impossible value so that other threads don't
            //   attempt to find neighbor nodes before we've inserted them all
            ATOMIC_STORE(node->claim_tid) = ParallelComponent<SeqGraph,SeqNode>::RELEASING; // TODO OPT: don't need this fancy jazz when forming initial resident seqgraph
            node_list.push_back( node );
            target->atomic_insertNode(node);
#else
            node->claim_tid = 0;
            target->insertNode(node);
#endif
        }
    }
    assert( this->node_count == 0 );
#ifdef VELOUR_TBB
    // all nodes present in the shared graph, reset the slice color on each node
    for(std::deque<SeqNode*>::iterator it = node_list.begin(); it != node_list.end(); ++it) {
        SeqNode *node = *it;
        ATOMIC_STORE(node->claim_tid) = ParallelComponent<SeqGraph,SeqNode>::UNCLAIMED;
    }
#endif
}

