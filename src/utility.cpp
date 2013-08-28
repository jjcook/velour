//
// This file is part of Velour and is distributed under the University of
// Illinois Open Source License.  See LICENSE.txt for details.
//

//
// utility.cpp
//
// the junk drawer
//

#include "types.h"


void report_unrecoverable_error(void)
{
    fprintf(stderr, "Unrecoverable error.  Exiting...");
}

//********************************************************************************
//***************************  Verification Routines  ****************************
//********************************************************************************

#ifdef SMALL_NODES
void verify_node(KmerNode * node, KmerGraph *hashtable, unsigned kmer_length)
{
    int double_kmer_length = kmer_length << 1;
    Kmer mask = (((Kmer)1) << double_kmer_length) - 1;  // a mask with 2*kmer_length bits
    Kmer kmer = node->kmer;
    Kmer rc_kmer = reverseComplement(kmer, kmer_length);
    char rightmost_base = KMER_GET_TAIL_BASE(kmer, kmer_length);
    char leftmost_base = KMER_GET_HEAD_BASE(kmer, kmer_length);
    KmerNode *node2;

    for (int i = 0 ; i < 4 ; ++ i) {
        // check on the left side
        int count = node->left_count[i];
        int color = node->left_color[i];
        assert( color == 0 || count != 0 );  // count must be non-zero if color is non-zero
        if (color == 0) {
            if (count > 0) {
                Kmer kmer2 = KMER_PREPEND(kmer, i, double_kmer_length, mask);
                node2 = hashtable->findNode(canonicalKmer(kmer2, kmer_length));
                assert( node2 != NULL );
                assert(cnorm(node2->right_count[static_cast<int>(rightmost_base)]) == cnorm(count));
                assert(node2->right_color[static_cast<int>(rightmost_base)] == 0);
            }	else if (count < 0) {
                Kmer kmer2 = KMER_APPEND(rc_kmer, i ^ 0x3, double_kmer_length);
                node2 = hashtable->findNode(canonicalKmer(kmer2, kmer_length));
                assert( node2 != NULL );
                assert(cnorm(node2->left_count[static_cast<int>( COMPLEMENT(rightmost_base) )]) == cnorm(count));
                assert(node2->left_color[static_cast<int>( COMPLEMENT(rightmost_base) )] == 0);
            }
        }

        // check on the right side
        count = node->right_count[i];
        color = node->right_color[i];
        assert( color == 0 || count != 0 );  // count must be non-zero if color is non-zero
        if (color == 0) {
            if (count > 0) {
                Kmer kmer2 = KMER_APPEND(kmer, i, double_kmer_length);
                node2 = hashtable->findNode(canonicalKmer(kmer2, kmer_length));
                assert( node2 != NULL );
                assert(cnorm(node2->left_count[static_cast<int>(leftmost_base)]) == cnorm(count));
                assert(node2->left_color[static_cast<int>(leftmost_base)] == 0);
            } else if (count < 0) {
                Kmer kmer2 = KMER_PREPEND(rc_kmer, i ^ 0x3, double_kmer_length, mask);
                node2 = hashtable->findNode(canonicalKmer(kmer2, kmer_length));
                assert( node2 != NULL );
                assert(cnorm(node2->right_count[static_cast<int>( COMPLEMENT(leftmost_base) )]) == cnorm(count));
                assert(node2->right_color[static_cast<int>( COMPLEMENT(leftmost_base) )] == 0);
            }
        }
    }
}
#else
void
verify_node_orig(kg_node_t * node, unsigned kmer_length) {
    assert( false && "TODO FIX! REVERSED KMER ENDIANNESS" );
  int double_kmer_length = kmer_length << 1;
  Kmer mask = (((Kmer)1) << double_kmer_length) - 1;  // a mask with 2*kmer_length bits
  Kmer kmer = node->kmer;
  Kmer rc_kmer = reverseComplement(kmer, kmer_length);
  char leftmost_base = (kmer >> (double_kmer_length - 2)) & 0x3;
  char rightmost_base = kmer & 0x3;

  for (int i = 0 ; i < 4 ; ++ i) {
    // check on the left side
    kg_node_t * node2 = node->left[i];
    int count = node->left_count[i];

    if (node2) {
      assert (count != 0);
      if (count > 0) {
	Kmer kmer2 = KMER_PREPEND(kmer, i, double_kmer_length, mask);
	assert(kmer2 == node2->kmer);
	assert(node2->right[(int)rightmost_base] == node);
	assert(node2->right_count[(int)rightmost_base] == count);
      } else {
	Kmer kmer2 = KMER_APPEND(rc_kmer, i ^ 0x3, double_kmer_length, mask);
	assert(kmer2 == node2->kmer);
	assert(node2->left[rightmost_base ^ 0x3] == node);
	assert(node2->left_count[rightmost_base ^ 0x3] == count);
      }
    } else {
      assert (count == 0);
    }


    // check on the right side
    node2 = node->right[i];
    count = node->right_count[i];

    if (node2) {
      assert (count != 0);
      if (count > 0) {
	Kmer kmer2 = KMER_APPEND(kmer, i, double_kmer_length, mask);
	assert(kmer2 == node2->kmer);
	assert(node2->left[(int)leftmost_base] == node);
	assert(node2->left_count[(int)leftmost_base] == count);
      } else {
	Kmer kmer2 = KMER_PREPEND(rc_kmer, i ^ 0x3, double_kmer_length, mask);
	assert(kmer2 == node2->kmer);
	assert(node2->right[leftmost_base ^ 0x3] == node);
	assert(node2->right_count[leftmost_base ^ 0x3] == count);
      }
    } else {
      assert (count == 0);
    }
  }
}

void
verify_node(kg_node_t * node, unsigned kmer_length) {
    assert( false && "TODO FIX! REVERSED KMER ENDIANNESS" );
  int double_kmer_length = kmer_length << 1;
  Kmer mask = (((Kmer)1) << double_kmer_length) - 1;  // a mask with 2*kmer_length bits
  Kmer kmer = node->kmer;
  Kmer rc_kmer = reverseComplement(kmer, kmer_length); //  + node->length - 1);
  char leftmost_base = (kmer >> (double_kmer_length - 2)) & 0x3;
  char rightmost_base = kmer & 0x3; // (kmer >> ((node->length-1)<<1)) & 0x3;
  Kmer lkmer = kmer; // >> ((node->length - 1) << 1);        // left kmer
  Kmer rrc_kmer = rc_kmer; // >> ((node->length - 1) << 1);  // right reverse-complement kmer

  for (int i = 0 ; i < 4 ; ++ i) {
    // check on the left side
    kg_node_t * node2 = node->left[i];
    int count = node->left_count[i];

    if (node2) {
      assert (count != 0);
      if (count > 0) {
	Kmer kmer2 = KMER_PREPEND(lkmer, i, double_kmer_length, mask);
	assert(kmer2 == (node2->kmer & mask));
	assert(node2->right[(int)rightmost_base] == node);
	assert(node2->right_count[(int)rightmost_base] == count);
      } else {
	Kmer kmer2 = KMER_APPEND(rc_kmer, i ^ 0x3, double_kmer_length, mask);
	assert(kmer2 == node2->kmer); //(node2->kmer >> ((node2->length-1) << 1)));
	assert(node2->left[rightmost_base ^ 0x3] == node);
	assert(node2->left_count[rightmost_base ^ 0x3] == count);
      }
    } else {
      assert (count == 0);
    }


    // check on the right side
    node2 = node->right[i];
    count = node->right_count[i];

    if (node2) {
      assert (count != 0);
      if (count > 0) {
	Kmer kmer2 = KMER_APPEND(kmer, i, double_kmer_length, mask);
	assert(kmer2 == node2->kmer); // (node2->kmer >> ((node2->length-1) << 1)));
	assert(node2->left[(int)leftmost_base] == node);
	assert(node2->left_count[(int)leftmost_base] == count);
      } else {
	Kmer kmer2 = KMER_PREPEND(rrc_kmer, i ^ 0x3, double_kmer_length, mask);
	assert(kmer2 == (node2->kmer & mask));
	assert(node2->right[leftmost_base ^ 0x3] == node);
	assert(node2->right_count[leftmost_base ^ 0x3] == count);
      }
    } else {
      assert (count == 0);
    }
  }
}
#endif

//********************************************************************************
//******************************  Graph Operations  ******************************
//********************************************************************************

int 
valid_single_successor(counter_t *counters) {
  int valid_dir = NO_SUCCESSORS;  // keep track of first valid direction found
  for (int i = 0 ; i < 4 ; ++ i) {
    if (counters[i] != 0) {
      if (valid_dir != NO_SUCCESSORS) {  // diverges at this point
		  return MULTIPLE_SUCCESSORS;
      }
      valid_dir = i;
    }
  }
  return valid_dir;
}  

/*
int
get_node_unique(sg_node_t *node, int count) {
  if (node == NULL) {
    return 0;
  } 
  int unique = node->unique;
  return (count < 0) ? -unique : unique;
} */ 

counter_t 
get_counter_sum(counter_t *counters) {
  counter_t count = 0;
  for (int i = 0 ; i < 4 ; ++ i) {
      counter_t temp = abs(counters[i]);
      if (temp == CLIP_SINGLECOPY_COUNTER_VALUE) temp = 1;   // XXX: inefficient ?
      count += temp;
  }
  return count;
}

counter_t 
get_counter_max(counter_t *counters) {
  counter_t max_value = 0;
  for (int i = 0 ; i < 4 ; ++ i) { 
	 counter_t temp = abs(counters[i]);  // absolute value
     if (temp == CLIP_SINGLECOPY_COUNTER_VALUE) temp = 1;  // XXX: inefficient ?
	 max_value = max(max_value, temp);
  }
  return max_value;
}

// given a pre-node "one", find out what nucleotide connects "one" to "two" 
// on a given side (specified by "left_not_right")
#ifndef SMALL_NODES
int
find_connection_index(kg_node_t *one, kg_node_t *two, bool left_not_right) {
  kg_node_t **pointers = left_not_right ? one->left : one->right;
  for (int i = 0 ; i < 4 ; ++ i) {
    if (pointers[i] == two) {
      return i;
    }
  }
  assert(0);  // doesn't connect
  return -1;
}

void
add_neighbor(kg_node_t *node, unsigned index, kg_node_t *rnode, counter_t count, bool left_not_right) {
  if (left_not_right) {
	 node->left[index] = rnode;
	 node->left_count[index] = count;
	 node->connections |= (1 << (index+LEFT_CONNECT_OFFSET));
  } else {
	 node->right[index] = rnode;
	 node->right_count[index] = count;
	 node->connections |= (1 << index);
  }
}

void
add_right_neighbor(kg_node_t *node, unsigned index, kg_node_t *rnode, counter_t count) {
  node->right[index] = rnode;
  node->right_count[index] = count;
  node->connections |= (1 << index);
}

void
remove_right_neighbor(kg_node_t *node, unsigned index) {
  node->right[index] = NULL;
  node->right_count[index] = 0;
  node->connections &= ~(1 << index);
}

void
remove_left_neighbor(kg_node_t *node, unsigned index) {
  node->left[index] = NULL;
  node->left_count[index] = 0;
  node->connections &= ~(1 << (index + LEFT_CONNECT_OFFSET));
}
#endif

